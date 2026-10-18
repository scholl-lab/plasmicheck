from __future__ import annotations

import base64
import bisect
import gzip
import hashlib
import json
import logging
import math
import os
import random
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd

from plasmicheck.config import get_config
from plasmicheck.resources import get_resource_path
from plasmicheck.version import __version__ as VERSION

from .plotting.colors import ASSIGNMENT_COLORS
from .utils import (
    categorize_value,
    encode_image_to_base64,
    ensure_plotly_assets,
    truncate_comma_separated,
    verdict_bounds,
)

_cfg = get_config()
DEFAULT_THRESHOLD: float = _cfg["default_threshold"]
UNCLEAR_RANGE: dict[str, float] = _cfg["unclear_range"]
PLOT_SAMPLE_REPORT: dict[str, Any] = _cfg["plot_sample_report"]
TEMPLATE_DIR: str = _cfg["paths"]["template_dir"]
LOGO_PATH: str = _cfg["paths"]["logo_path"]
DOWNSAMPLE_LIMIT: int = PLOT_SAMPLE_REPORT.get("max_rows_for_interactive_plots", 5000)

IGV_VERSION = "3.0.2"
IGV_ASSET_RELPATH = "static/js/igv-3.0.2.min.js.gz"
IGV_GZIP_SHA256 = "0d8b512654b2ef588009453c403c8f1329dce88eedca90ba9e60888af6b2f79f"
IGV_SHA256 = "ab1aa79c514ee3a0d66a0ffc788b6d37803910e62cf6d114d9b2909d96b5e790"

# Setup logging for the libraries used in this script
logging.getLogger("jinja2").setLevel(logging.ERROR)
logging.getLogger("fontTools").setLevel(logging.ERROR)


def load_data(reads_assignment_file: str, summary_file: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    import pandas as pd

    logging.info(f"Loading data from {reads_assignment_file} and {summary_file}")
    reads_df = pd.read_csv(
        reads_assignment_file,
        sep="\t",
        usecols=lambda c: c in {"AssignedTo", "PlasmidScore", "HumanScore"},
    )
    summary_df = pd.read_csv(summary_file, sep="\t")
    return reads_df, summary_df


def downsample_data(df: pd.DataFrame, limit: int) -> tuple[pd.DataFrame, bool]:
    if len(df) > limit:
        logging.info(f"Downsampling data from {len(df)} to {limit} points.")
        df = df.sample(n=limit, random_state=1)  # Random downsample to the limit
        downsampled = True
    else:
        downsampled = False
    return df, downsampled


_CATEGORY_ORDER: list[str] = ["Plasmid", "Human", "Tied", "Backbone_Only", "Ambiguous"]

_PLOTLY_FONT = '-apple-system, BlinkMacSystemFont, "Inter", "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif'

_PLOTLY_LAYOUT: dict[str, Any] = {
    "autosize": True,
    "font": {"family": _PLOTLY_FONT, "size": 13, "color": "#0f172a"},
    "paper_bgcolor": "white",
    "plot_bgcolor": "white",
    "margin": {"l": 64, "r": 24, "t": 36, "b": 52, "pad": 4},
    "xaxis": {
        "showgrid": True,
        "gridwidth": 1,
        "gridcolor": "#f1f5f9",
        "showline": True,
        "linewidth": 1,
        "linecolor": "#cbd5e1",
        "zeroline": False,
        "tickfont": {"size": 12, "color": "#475569"},
        "title_font": {"size": 13, "color": "#334155"},
    },
    "yaxis": {
        "showgrid": True,
        "gridwidth": 1,
        "gridcolor": "#f1f5f9",
        "showline": True,
        "linewidth": 1,
        "linecolor": "#cbd5e1",
        "zeroline": False,
        "tickfont": {"size": 12, "color": "#475569"},
        "title_font": {"size": 13, "color": "#334155"},
    },
    "legend": {
        "orientation": "h",
        "yanchor": "bottom",
        "y": 1.02,
        "xanchor": "right",
        "x": 1.0,
        "bgcolor": "rgba(255,255,255,0.92)",
        "bordercolor": "#e2e8f0",
        "borderwidth": 1,
        "font": {"size": 12, "color": "#334155"},
        "title_text": "Assigned To  ",
    },
    "hoverlabel": {
        "bgcolor": "#0f172a",
        "bordercolor": "#334155",
        "font": {"size": 12, "family": _PLOTLY_FONT, "color": "#f8fafc"},
    },
}

_PLOTLY_CONFIG: dict[str, Any] = {
    "showSendToCloud": False,
    "displayModeBar": "hover",
    "displaylogo": False,
    "modeBarButtonsToRemove": ["pan2d", "lasso2d", "select2d", "autoScale2d"],
    "responsive": True,
    "toImageButtonOptions": {
        "format": "png",
        "filename": "plasmicheck_chart",
        "height": 800,
        "width": 1200,
        "scale": 2,
    },
}


def generate_plots(
    reads_df: pd.DataFrame,
    output_folder: str,
    static_report: bool = False,
    plot_backend: str = "matplotlib",
    reports_dir: str | None = None,
    total_reads: int | None = None,
    threshold: float = DEFAULT_THRESHOLD,
    unclear: dict[str, float] = UNCLEAR_RANGE,
) -> tuple[str, str | None, str, str | None]:
    """Generate interactive and optional static plots.

    Returns (boxplot_div_html, boxplot_png_path, scatter_div_html, scatter_png_path).
    The div HTML strings are for embedding in the main report template.
    Standalone viewable plots are written to reports_dir/plots/.
    """
    import plotly.express as px

    logging.info("Generating plots")
    width = PLOT_SAMPLE_REPORT.get("figsize", {}).get("width", 1000)
    height = PLOT_SAMPLE_REPORT.get("figsize", {}).get("height", 400)
    if not isinstance(width, int):
        width = 1000
    if not isinstance(height, int):
        height = 400

    # Standalone plots go to reports/plots/ subdirectory
    effective_reports_dir = reports_dir if reports_dir else output_folder
    standalone_plots_dir = os.path.join(effective_reports_dir, "plots")
    os.makedirs(standalone_plots_dir, exist_ok=True)

    # Consistent category ordering across all plots
    cat_order = {"AssignedTo": _CATEGORY_ORDER}

    # ── Box plot ──────────────────────────────────────────
    fig_box = px.box(
        reads_df,
        x="AssignedTo",
        y="PlasmidScore",
        points="all",
        color="AssignedTo",
        color_discrete_map=ASSIGNMENT_COLORS,
        category_orders=cat_order,
        labels={
            "PlasmidScore": PLOT_SAMPLE_REPORT["box_plot_y_label"],
            "AssignedTo": "Read Assignment Category",
        },
        height=height,
    )
    fig_box.update_layout(**_PLOTLY_LAYOUT, showlegend=False, hovermode="closest")
    fig_box.update_xaxes(showgrid=False)
    fig_box.add_hline(y=0, line_dash="dot", line_color="#94a3b8", line_width=1)
    fig_box.update_traces(
        marker={"size": 3.5, "opacity": 0.45, "line_width": 0},
        line={"width": 1.6},
    )

    # Get div HTML for embedding in main report (no plotly.js, no full HTML)
    boxplot_div_html: str = fig_box.to_html(
        include_plotlyjs=False,
        full_html=False,
        config=_PLOTLY_CONFIG,
    )

    # Write standalone HTML for direct viewing (CDN plotly)
    boxplot_standalone = os.path.join(
        standalone_plots_dir,
        PLOT_SAMPLE_REPORT["output_box_plot_filename"].replace(".png", ".html"),
    )
    fig_box.write_html(
        boxplot_standalone,
        include_plotlyjs="cdn",
        full_html=True,
        config=_PLOTLY_CONFIG,
    )

    # ── Scatter plot ──────────────────────────────────────
    fig_scatter = px.scatter(
        reads_df,
        x="PlasmidScore",
        y="HumanScore",
        color="AssignedTo",
        color_discrete_map=ASSIGNMENT_COLORS,
        category_orders=cat_order,
        labels={
            "PlasmidScore": PLOT_SAMPLE_REPORT["scatter_plot_x_label"],
            "HumanScore": PLOT_SAMPLE_REPORT["scatter_plot_y_label"],
        },
        height=height,
        render_mode="webgl",
    )
    fig_scatter.update_layout(**_PLOTLY_LAYOUT, showlegend=True, hovermode="closest")
    fig_scatter.update_traces(marker={"size": 4.5, "opacity": 0.55, "line_width": 0})

    # Diagonal reference line (equal scores y = x across full data range)
    import plotly.graph_objects as go

    p_max = float(reads_df["PlasmidScore"].max()) if not reads_df.empty else 1.0
    h_max = float(reads_df["HumanScore"].max()) if not reads_df.empty else 1.0
    p_min = float(reads_df["PlasmidScore"].min()) if not reads_df.empty else 0.0
    h_min = float(reads_df["HumanScore"].min()) if not reads_df.empty else 0.0
    import math

    score_max = max(
        1.0,
        0.0 if math.isnan(p_max) else p_max,
        0.0 if math.isnan(h_max) else h_max,
    )
    score_min = min(
        0.0,
        0.0 if math.isnan(p_min) else p_min,
        0.0 if math.isnan(h_min) else h_min,
    )
    fig_scatter.add_trace(
        go.Scattergl(
            x=[score_min, score_max],
            y=[score_min, score_max],
            mode="lines",
            name="Equal Score (y = x)",
            line={"color": "#64748b", "width": 1.5, "dash": "dash"},
            showlegend=True,
            hoverinfo="skip",
        )
    )

    # Get div HTML for embedding in main report
    scatter_div_html: str = fig_scatter.to_html(
        include_plotlyjs=False,
        full_html=False,
        config=_PLOTLY_CONFIG,
    )

    # Write standalone HTML for direct viewing (CDN plotly)
    scatter_standalone = os.path.join(
        standalone_plots_dir,
        PLOT_SAMPLE_REPORT["output_scatter_plot_filename"].replace(".png", ".html"),
    )
    fig_scatter.write_html(
        scatter_standalone,
        include_plotlyjs="cdn",
        full_html=True,
        config=_PLOTLY_CONFIG,
    )

    # Only generate PNGs when static_report=True
    boxplot_filename_png = None
    scatter_filename_png = None
    if static_report:
        boxplot_filename_png = os.path.join(
            standalone_plots_dir, PLOT_SAMPLE_REPORT["output_box_plot_filename"]
        )
        scatter_filename_png = os.path.join(
            standalone_plots_dir, PLOT_SAMPLE_REPORT["output_scatter_plot_filename"]
        )

        if plot_backend == "matplotlib":
            from .plotting.matplotlib_backend import (
                generate_boxplot_matplotlib,
                generate_scatter_matplotlib,
            )

            # Generate plots using matplotlib
            generate_boxplot_matplotlib(
                reads_df,
                boxplot_filename_png,
                total_reads=total_reads,
                threshold=threshold,
                unclear_range=unclear,
            )
            generate_scatter_matplotlib(reads_df, scatter_filename_png, total_reads=total_reads)
        else:
            # Default plotly/kaleido backend
            import kaleido  # type: ignore[import-untyped]

            kaleido.start_sync_server()

            fig_box.write_image(boxplot_filename_png, width=width, height=height)
            fig_scatter.write_image(scatter_filename_png, width=width, height=height)

    logging.info("Plots generated.")

    return (
        boxplot_div_html,
        boxplot_filename_png,
        scatter_div_html,
        scatter_filename_png,
    )


_ensure_plotly_assets = ensure_plotly_assets


def extract_verdict_from_summary(summary_df: pd.DataFrame) -> str:
    verdict_row: pd.DataFrame = summary_df[summary_df["Category"] == "Verdict"]
    if not verdict_row.empty:
        return str(verdict_row["Count"].values[0])
    else:
        return "Verdict not found in summary file"


def _derive_verdict_info(verdict: str) -> tuple[str, str]:
    """Return (verdict_short, verdict_class) from the full verdict string."""
    if "not contaminated" in verdict.lower():
        return "Clean", "clean"
    elif "unclear" in verdict.lower():
        return "Unclear", "unclear"
    else:
        return "Contaminated", "contaminated"


def _compute_gauge_position(
    ratio: float,
    unclear_lower: float,
    unclear_upper: float,
) -> float:
    """Map ratio to 0-100 gauge position across three equal-width visual zones.

    Zones: [0, unclear_lower] -> 0-33.33%, [unclear_lower, unclear_upper] -> 33.33-66.67%,
    [unclear_upper, 10*unclear_upper] -> 66.67-95%, clamped to [3.0, 95.0] for badge padding.
    """
    if ratio <= 0:
        return 3.0
    if unclear_lower <= 0:
        unclear_lower = 0.001
    if unclear_upper <= unclear_lower:
        if ratio <= unclear_upper:
            return max(3.0, min(50.0, (ratio / max(unclear_upper, 0.001)) * 50.0))
        upper_max = max(10.0 * unclear_upper, unclear_upper + 1.0)
        return min(95.0, 50.0 + ((ratio - unclear_upper) / (upper_max - unclear_upper)) * 45.0)

    if ratio < unclear_lower:
        return max(3.0, (ratio / unclear_lower) * (100.0 / 3.0))
    if ratio <= unclear_upper:
        return (100.0 / 3.0) + ((ratio - unclear_lower) / (unclear_upper - unclear_lower)) * (
            100.0 / 3.0
        )
    upper_max = 10.0 * unclear_upper
    if ratio <= upper_max:
        return (200.0 / 3.0) + ((ratio - unclear_upper) / (upper_max - unclear_upper)) * (
            95.0 - (200.0 / 3.0)
        )
    return 95.0


def _parse_mismatches(raw: str) -> tuple[int, int]:
    """Parse mismatches_near_insert from its raw string representation.

    Returns (with_mismatches, without_mismatches).
    """
    import ast

    if raw in ("N/A", ""):
        return 0, 0
    try:
        d = ast.literal_eval(raw) if isinstance(raw, str) else raw
        return (
            int(d.get("with_mismatches_or_clipping", 0)),
            int(d.get("without_mismatches_or_clipping", 0)),
        )
    except (ValueError, SyntaxError, AttributeError):
        return 0, 0


def _truncate_list_filter(value: str, max_items: int = 5) -> str:
    """Truncate comma-separated list for display. Full value in title attribute."""
    if not value or value in ("N/A", "NA"):
        return value
    return truncate_comma_separated(value, max_items)


def _load_verified_igv_b64() -> str:
    """Load vendored igv.min.js.gz, verify compressed & decompressed SHA-256, return base64."""
    igv_path = get_resource_path(IGV_ASSET_RELPATH)
    gz_bytes = Path(igv_path).read_bytes()
    actual_gz_sha = hashlib.sha256(gz_bytes).hexdigest()
    if actual_gz_sha != IGV_GZIP_SHA256:
        raise RuntimeError(
            f"Vendored {IGV_ASSET_RELPATH} compressed SHA-256 mismatch: "
            f"expected {IGV_GZIP_SHA256}, got {actual_gz_sha}"
        )
    raw_bytes = gzip.decompress(gz_bytes)
    actual_sha = hashlib.sha256(raw_bytes).hexdigest()
    if actual_sha != IGV_SHA256:
        raise RuntimeError(
            f"Vendored {IGV_ASSET_RELPATH} decompressed SHA-256 mismatch: "
            f"expected {IGV_SHA256}, got {actual_sha}"
        )
    return base64.b64encode(gz_bytes).decode("ascii")


def _encode_bedgraph_uri(seq_name: str, depth_arr: Any) -> str:
    """Encode a per-base depth array into a compact RLE bedGraph base64 data URI."""
    import numpy as np

    n = len(depth_arr)
    if n == 0:
        return ""
    step = 5 if n > 2000 else 1
    lines: list[str] = []
    run_start = 0
    first_chunk = depth_arr[0 : min(step, n)]
    run_val = round(float(np.mean(first_chunk)))
    for pos in range(step, n, step):
        chunk = depth_arr[pos : min(pos + step, n)]
        val = round(float(np.mean(chunk)))
        if val != run_val:
            lines.append(f"{seq_name}\t{run_start}\t{pos}\t{run_val}")
            run_start = pos
            run_val = val
    lines.append(f"{seq_name}\t{run_start}\t{n}\t{run_val}")
    raw = ("\n".join(lines) + "\n").encode("utf-8")
    return "data:application/octet-stream;base64," + base64.b64encode(raw).decode("ascii")


def _build_igv_context(
    output_folder: str,
    reports_dir: str | None = None,
    anonymize: bool = False,
    max_reads: int = 3200,
    max_raw_bam_bytes: int = 262_144,
) -> dict[str, Any]:
    """Build self-contained igv.js alignment inspector context from pipeline artifacts."""
    report_dir = Path(reports_dir if reports_dir else output_folder)
    base_dir = Path(output_folder)
    session_xml_path = base_dir / "igv_session.xml"
    session_rel = (
        os.path.relpath(session_xml_path, report_dir) if session_xml_path.is_file() else ""
    )

    base_ctx: dict[str, Any] = {
        "igv_embedded": False,
        "igv_version": IGV_VERSION,
        "igv_sha256": IGV_SHA256,
        "igv_absence_reason": "artifacts_absent",
        "igv_absence_detail": (
            "Alignment BAM and plasmid reference FASTA are not present in this report directory "
            "(run with --keep-intermediate to retain desktop igv_session.xml and raw BAMs)."
        ),
        "igv_loci": [],
        "igv_session_relpath": session_rel,
        "igv_downsampled": False,
        "igv_read_count": 0,
        "igv_total_primary_reads": 0,
        "igv_splice_junction_count": 0,
        "igv_human_genomic_locus": "",
        "igv_payload": "",
    }

    if anonymize:
        base_ctx["igv_absence_reason"] = "anonymized"
        base_ctx["igv_absence_detail"] = (
            "Interactive read-level alignment browser is withheld because --anonymize was enabled "
            "(raw read sequences, QNAMEs, and intermediate BAMs are redacted)."
        )
        return base_ctx

    bam_candidates = [
        base_dir / "alignments" / "plasmid_alignment.bam",
        base_dir / "plasmid_alignment.bam",
    ]
    bam_path = next((p for p in bam_candidates if p.is_file() and p.stat().st_size > 0), None)

    ref_dir = base_dir / "references"
    fasta_candidates = sorted(ref_dir.glob("*.fasta")) if ref_dir.is_dir() else []
    if not fasta_candidates:
        fasta_candidates = sorted(base_dir.glob("*.fasta"))
    fasta_path = next(
        (
            p
            for p in fasta_candidates
            if not p.name.endswith("spliced_reference.fasta") and p.stat().st_size > 0
        ),
        None,
    )

    if bam_path is None or fasta_path is None:
        return base_ctx

    try:
        import numpy as np
        import pysam

        # Verify BGZF magic header before opening with pysam so mock text files fall back cleanly
        if bam_path.read_bytes()[:4] != b"\x1f\x8b\x08\x04":
            return base_ctx

        fasta_bytes = fasta_path.read_bytes()
        seq_name = "plasmid"
        seq_len = 0
        for line in fasta_bytes.decode("utf-8", errors="replace").splitlines():
            line_s = line.strip()
            if line_s.startswith(">"):
                seq_name = line_s[1:].split()[0]
            elif line_s:
                seq_len += len(line_s)
        if seq_len <= 0:
            return base_ctx

        cdna_candidates = [
            base_dir / "data" / "cDNA_positions.txt",
            base_dir / "cDNA_positions.txt",
        ]
        cdna_path = next((p for p in cdna_candidates if p.is_file()), None)
        insert_start: int | None = None
        insert_end: int | None = None
        exon_lengths: list[int] = []
        intron_lengths: list[int] = []
        strand = "+"
        chrom = ""
        genomic_start = ""
        genomic_end = ""
        if cdna_path is not None:
            for line in cdna_path.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith("cDNA start position in plasmid:"):
                    val = line.split(":", 1)[1].strip()
                    if val not in ("N/A", "NA", ""):
                        insert_start = int(val)
                elif line.startswith("cDNA end position in plasmid:"):
                    val = line.split(":", 1)[1].strip()
                    if val not in ("N/A", "NA", ""):
                        insert_end = int(val)
                elif line.startswith("ExonLengths:"):
                    val = line.split(":", 1)[1].strip()
                    if val not in ("N/A", "NA", ""):
                        exon_lengths = [int(x) for x in val.split(",") if x.strip().isdigit()]
                elif line.startswith("IntronLengths:"):
                    val = line.split(":", 1)[1].strip()
                    if val not in ("N/A", "NA", ""):
                        intron_lengths = [int(x) for x in val.split(",") if x.strip().isdigit()]
                elif line.startswith("Strand:"):
                    val = line.split(":", 1)[1].strip()
                    if val in ("+", "-"):
                        strand = val
                elif line.startswith("Chromosome:"):
                    val = line.split(":", 1)[1].strip()
                    if val not in ("N/A", "NA", ""):
                        chrom = val
                elif line.startswith("GenomicStart:"):
                    val = line.split(":", 1)[1].strip()
                    if val not in ("N/A", "NA", ""):
                        genomic_start = val
                elif line.startswith("GenomicEnd:"):
                    val = line.split(":", 1)[1].strip()
                    if val not in ("N/A", "NA", ""):
                        genomic_end = val

        human_genomic_locus = ""
        if chrom and genomic_start.isdigit() and genomic_end.isdigit():
            human_genomic_locus = f"{chrom}:{int(genomic_start):,}-{int(genomic_end):,} ({strand})"

        # Orient exon and intron lengths along plasmid coordinate axis (5' -> 3' on plasmid)
        exon_lengths_plasmid = exon_lengths[::-1] if strand == "-" else list(exon_lengths)
        intron_lengths_plasmid = intron_lengths[::-1] if strand == "-" else list(intron_lengths)

        has_insert = (
            insert_start is not None and insert_end is not None and insert_end > insert_start
        )
        junc_positions: list[int] = []
        junc_intron_bp: list[int | None] = []
        bed_lines: list[str] = []

        if has_insert and insert_start is not None and insert_end is not None:
            bed_lines.append(
                f"{seq_name}\t{insert_start}\t{insert_end + 1}\tcDNA_Insert\t1000\t+\t"
                f"{insert_start}\t{insert_end + 1}\t37,99,235"
            )
            b5_lo = max(0, insert_start - 10)
            b5_hi = min(seq_len, insert_start + 10)
            b3_lo = max(0, insert_end - 10)
            b3_hi = min(seq_len, insert_end + 10)
            bed_lines.append(
                f"{seq_name}\t{b5_lo}\t{b5_hi}\t5p_Vector_Junction\t1000\t+\t{b5_lo}\t{b5_hi}\t124,58,237"
            )
            bed_lines.append(
                f"{seq_name}\t{b3_lo}\t{b3_hi}\t3p_Vector_Junction\t1000\t+\t{b3_lo}\t{b3_hi}\t124,58,237"
            )
            if len(exon_lengths_plasmid) > 1 and sum(exon_lengths_plasmid) <= (
                insert_end - insert_start + 5
            ):
                cursor = insert_start
                for idx, ex_len in enumerate(exon_lengths_plasmid, start=1):
                    ex_end = min(seq_len, cursor + ex_len)
                    ex_color = "13,148,136" if idx % 2 == 1 else "15,118,110"
                    bed_lines.append(
                        f"{seq_name}\t{cursor}\t{ex_end}\tExon_{idx}\t1000\t+\t"
                        f"{cursor}\t{ex_end}\t{ex_color}"
                    )
                    if idx < len(exon_lengths_plasmid):
                        j_pos = ex_end
                        intron_bp = (
                            intron_lengths_plasmid[idx - 1]
                            if (idx - 1) < len(intron_lengths_plasmid)
                            else None
                        )
                        junc_positions.append(j_pos)
                        junc_intron_bp.append(intron_bp)
                        j_lo = max(0, j_pos - 2)
                        j_hi = min(seq_len, j_pos + 2)
                        j_label = (
                            f"J{idx}_Exon{idx}_{idx + 1}_Intron_{intron_bp}bp"
                            if intron_bp is not None
                            else f"J{idx}_Exon{idx}_{idx + 1}"
                        )
                        bed_lines.append(
                            f"{seq_name}\t{j_lo}\t{j_hi}\t{j_label}\t1000\t+\t"
                            f"{j_lo}\t{j_hi}\t220,38,38"
                        )
                    cursor = ex_end

        n_junc = len(junc_positions)
        n_diag_loci = n_junc + (2 if has_insert else 0)
        if max_reads <= 20:
            locus_quota = max(1, max_reads // max(1, n_diag_loci))
        else:
            locus_quota = max(24, min(120, (max_reads * 3 // 4) // max(1, n_diag_loci)))

        rng = random.Random(42)

        def _res_add(
            buf: list[tuple[Any, str]], item: tuple[Any, str], seen: int, cap: int
        ) -> None:
            if cap <= 0:
                return
            if len(buf) < cap:
                buf.append(item)
            else:
                j = rng.randint(0, seen - 1)
                if j < cap:
                    buf[j] = item

        b5_vi_count = 0
        b3_vi_count = 0
        j_ee_counts = [0] * n_junc
        j_ei_counts = [0] * n_junc
        backbone_count = 0
        shared_count = 0

        b5_bucket: list[tuple[Any, str]] = []
        b3_bucket: list[tuple[Any, str]] = []
        j_ee_buckets: list[list[tuple[Any, str]]] = [[] for _ in range(n_junc)]
        j_ei_buckets: list[list[tuple[Any, str]]] = [[] for _ in range(n_junc)]
        backbone_bucket: list[tuple[Any, str]] = []
        shared_bucket: list[tuple[Any, str]] = []

        total_primary_reads = 0
        diff_depth = np.zeros(seq_len + 1, dtype=np.int64)

        with pysam.AlignmentFile(str(bam_path), "rb") as in_bam:
            for r in in_bam.fetch(until_eof=True):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate:
                    continue
                total_primary_reads += 1
                rstart = r.reference_start
                rend = r.reference_end if r.reference_end is not None else rstart
                if 0 <= rstart < seq_len:
                    diff_depth[rstart] += 1
                    diff_depth[max(rstart, min(seq_len, rend))] -= 1

                matched_locus = False
                if has_insert and insert_start is not None and insert_end is not None:
                    cigar = r.cigartuples or []
                    nm_val = int(r.get_tag("NM")) if r.has_tag("NM") else 0
                    has_long_clip = any(op in (4, 5) and ln > 6 for op, ln in cigar)

                    # 1. Vector-Insert boundaries (5' at insert_start, 3' at insert_end + 1)
                    if not has_long_clip and nm_val <= 2:
                        b5 = insert_start
                        if (
                            0 < b5 < seq_len
                            and rstart + 10 <= b5 <= rend - 10
                            and r.get_overlap(b5 - 10, b5 + 10) == 20
                        ):
                            b5_vi_count += 1
                            _res_add(b5_bucket, (r, "Vector_Insert"), b5_vi_count, locus_quota)
                            matched_locus = True
                        b3 = insert_end + 1
                        if (
                            0 < b3 < seq_len
                            and rstart + 10 <= b3 <= rend - 10
                            and r.get_overlap(b3 - 10, b3 + 10) == 20
                        ):
                            b3_vi_count += 1
                            if not matched_locus:
                                _res_add(b3_bucket, (r, "Vector_Insert"), b3_vi_count, locus_quota)
                                matched_locus = True

                    # 2. Exon-Exon splice junction spanning reads (Plasmid_ExonExon)
                    if n_junc > 0 and not has_long_clip and nm_val <= 2 and (rend - rstart) >= 30:
                        idx = bisect.bisect_left(junc_positions, rstart + 15)
                        while idx < n_junc and junc_positions[idx] <= rend - 15:
                            jp = junc_positions[idx]
                            if r.get_overlap(jp - 15, jp + 15) == 30:
                                j_ee_counts[idx] += 1
                                if not matched_locus:
                                    _res_add(
                                        j_ee_buckets[idx],
                                        (r, "Plasmid_ExonExon"),
                                        j_ee_counts[idx],
                                        locus_quota,
                                    )
                                    matched_locus = True
                            idx += 1

                    # 3. Human Exon-Intron soft-clipped reads at splice junctions (Human_ExonIntron)
                    if n_junc > 0 and cigar:
                        if cigar[0][0] == 4 and cigar[0][1] >= 15:
                            idx = bisect.bisect_left(junc_positions, rstart - 2)
                            while idx < n_junc and junc_positions[idx] <= rstart + 2:
                                jp = junc_positions[idx]
                                if (r.get_overlap(jp, min(seq_len, jp + 15)) or 0) >= 15:
                                    j_ei_counts[idx] += 1
                                    if not matched_locus:
                                        _res_add(
                                            j_ei_buckets[idx],
                                            (r, "Human_ExonIntron"),
                                            j_ei_counts[idx],
                                            locus_quota,
                                        )
                                        matched_locus = True
                                idx += 1
                        if cigar[-1][0] == 4 and cigar[-1][1] >= 15:
                            idx = bisect.bisect_left(junc_positions, rend - 2)
                            while idx < n_junc and junc_positions[idx] <= rend + 2:
                                jp = junc_positions[idx]
                                if (r.get_overlap(max(0, jp - 15), jp) or 0) >= 15:
                                    j_ei_counts[idx] += 1
                                    if not matched_locus:
                                        _res_add(
                                            j_ei_buckets[idx],
                                            (r, "Human_ExonIntron"),
                                            j_ei_counts[idx],
                                            locus_quota,
                                        )
                                        matched_locus = True
                                idx += 1

                if not matched_locus:
                    if (
                        has_insert
                        and insert_start is not None
                        and insert_end is not None
                        and (rend <= insert_start or rstart > insert_end)
                    ):
                        backbone_count += 1
                        _res_add(backbone_bucket, (r, "Backbone"), backbone_count, max_reads)
                    else:
                        shared_count += 1
                        _res_add(shared_bucket, (r, "Shared_Exon"), shared_count, max_reads)

            # Assemble stratified diagnostic read sample
            selected_pairs: list[tuple[Any, str]] = []
            b5_sel = b5_bucket[:locus_quota]
            b3_sel = b3_bucket[:locus_quota]
            selected_pairs.extend(b5_sel)
            selected_pairs.extend(b3_sel)

            j_embedded_counts = [0] * n_junc
            for idx in range(n_junc):
                ee_list = j_ee_buckets[idx]
                ei_list = j_ei_buckets[idx]
                half = max(1, locus_quota // 2)
                take_ee = min(len(ee_list), half)
                take_ei = min(len(ei_list), locus_quota - take_ee)
                if take_ee + take_ei < locus_quota:
                    if len(ee_list) > take_ee:
                        take_ee = min(len(ee_list), locus_quota - take_ei)
                    elif len(ei_list) > take_ei:
                        take_ei = min(len(ei_list), locus_quota - take_ee)
                j_sel = ee_list[:take_ee] + ei_list[:take_ei]
                j_embedded_counts[idx] = len(j_sel)
                selected_pairs.extend(j_sel)

            if len(selected_pairs) > max_reads:
                step = max(1, len(selected_pairs) // max_reads)
                selected_pairs = selected_pairs[::step][:max_reads]
            elif len(selected_pairs) < max_reads:
                rem = max_reads - len(selected_pairs)
                half_rem = rem // 2
                take_bb = min(len(backbone_bucket), half_rem)
                take_sh = min(len(shared_bucket), rem - take_bb)
                if take_bb + take_sh < rem and len(backbone_bucket) > take_bb:
                    take_bb = min(len(backbone_bucket), rem - take_sh)
                selected_pairs.extend(backbone_bucket[:take_bb])
                selected_pairs.extend(shared_bucket[:take_sh])

            selected_pairs.sort(key=lambda item: (item[0].reference_id, item[0].reference_start))
            mapped_reads = len(selected_pairs)
            downsampled = total_primary_reads > mapped_reads

            with tempfile.NamedTemporaryFile(suffix=".bam", delete=False) as tmp_bam:
                tmp_bam_name = tmp_bam.name
            try:
                with pysam.AlignmentFile(tmp_bam_name, "wb", header=in_bam.header) as out_bam:
                    for r, pc_cat in selected_pairs:
                        r.set_tag("PC", pc_cat, value_type="Z")
                        out_bam.write(r)
                bam_bytes = Path(tmp_bam_name).read_bytes()
            finally:
                if os.path.exists(tmp_bam_name):
                    os.remove(tmp_bam_name)

        # Build 100% full-BAM depth array (prefer C count_coverage when indexed, fallback to diff_depth)
        try:
            from .coverage_metrics import get_depth_array

            def _primary_cb(read: Any) -> bool:
                return not (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_duplicate
                )

            depth_arr = get_depth_array(
                str(bam_path), seq_name, 0, seq_len, read_callback=_primary_cb
            )
            if len(depth_arr) != seq_len:
                depth_arr = np.cumsum(diff_depth[:-1])
        except Exception:
            depth_arr = np.cumsum(diff_depth[:-1])

        bedgraph_uri = _encode_bedgraph_uri(seq_name, depth_arr)

        # Build rich diagnostic loci metadata with full-BAM and embedded read counts
        loci: list[dict[str, Any]] = [
            {
                "id": "full",
                "short_label": "Full Construct",
                "label": "Full Construct",
                "type": "Plasmid Reference",
                "locus": f"{seq_name}:1-{seq_len}",
                "note": f"{seq_len:,} bp circular construct",
                "ee_count": None,
                "ei_count": None,
                "embedded_count": mapped_reads,
                "is_quick_jump": True,
            }
        ]

        if has_insert and insert_start is not None and insert_end is not None:
            pad_lo = max(1, insert_start - 150)
            pad_hi = min(seq_len, insert_end + 150)
            insert_embedded = sum(
                1
                for r, _ in selected_pairs
                if (r.reference_end or r.reference_start) > insert_start
                and r.reference_start <= insert_end
            )
            loci.append(
                {
                    "id": "insert",
                    "short_label": "cDNA Insert",
                    "label": "cDNA Insert",
                    "type": "Insert Region",
                    "locus": f"{seq_name}:{pad_lo}-{pad_hi}",
                    "note": (
                        f"{insert_end - insert_start + 1:,} bp cDNA insert ({len(exon_lengths_plasmid)} exons)"
                        if exon_lengths_plasmid
                        else f"{insert_end - insert_start + 1:,} bp insert (±150 bp flank)"
                    ),
                    "ee_count": sum(j_ee_counts) if n_junc > 0 else None,
                    "ei_count": sum(j_ei_counts) if n_junc > 0 else None,
                    "embedded_count": insert_embedded,
                    "is_quick_jump": True,
                }
            )
            loci.append(
                {
                    "id": "b5",
                    "short_label": "5' Boundary",
                    "label": "5' Vector-Insert Boundary",
                    "type": "Vector Junction",
                    "locus": f"{seq_name}:{max(1, insert_start - 120)}-{min(seq_len, insert_start + 120)}",
                    "note": f"Coordinate {insert_start:,} bp (Vector → cDNA)",
                    "ee_count": b5_vi_count,
                    "ei_count": None,
                    "embedded_count": len(b5_sel),
                    "is_quick_jump": True,
                }
            )
            loci.append(
                {
                    "id": "b3",
                    "short_label": "3' Boundary",
                    "label": "3' Vector-Insert Boundary",
                    "type": "Vector Junction",
                    "locus": f"{seq_name}:{max(1, insert_end - 120)}-{min(seq_len, insert_end + 120)}",
                    "note": f"Coordinate {insert_end:,} bp (cDNA → Vector)",
                    "ee_count": b3_vi_count,
                    "ei_count": None,
                    "embedded_count": len(b3_sel),
                    "is_quick_jump": True,
                }
            )
            for idx, j_pos in enumerate(junc_positions, start=1):
                intron_bp = junc_intron_bp[idx - 1]
                note_str = (
                    f"Coord {j_pos:,} bp · Spliced human intron: {intron_bp:,} bp"
                    if intron_bp is not None
                    else f"Coordinate {j_pos:,} bp"
                )
                loci.append(
                    {
                        "id": f"j{idx}",
                        "short_label": f"J{idx}",
                        "label": f"Splice Junction {idx} (Exon {idx}/{idx + 1})",
                        "type": "Exon-Exon Junction",
                        "locus": f"{seq_name}:{max(1, j_pos - 100)}-{min(seq_len, j_pos + 100)}",
                        "note": note_str,
                        "ee_count": j_ee_counts[idx - 1],
                        "ei_count": j_ei_counts[idx - 1],
                        "embedded_count": j_embedded_counts[idx - 1],
                        "is_quick_jump": False,
                    }
                )

        igv_gz_b64 = _load_verified_igv_b64()
        fasta_uri = "data:application/octet-stream;base64," + base64.b64encode(fasta_bytes).decode(
            "ascii"
        )
        bam_uri = "data:application/octet-stream;base64," + base64.b64encode(bam_bytes).decode(
            "ascii"
        )
        bed_uri = (
            "data:application/octet-stream;base64,"
            + base64.b64encode(("\n".join(bed_lines) + "\n").encode("utf-8")).decode("ascii")
            if bed_lines
            else ""
        )

        tracks: list[dict[str, Any]] = []
        if bed_uri:
            tracks.append(
                {
                    "name": "Plasmid Construct & Human Splice Sites",
                    "type": "annotation",
                    "format": "bed",
                    "url": bed_uri,
                    "indexed": False,
                    "displayMode": "EXPANDED",
                    "height": 74,
                }
            )
        if bedgraph_uri:
            tracks.append(
                {
                    "name": f"Full-BAM Depth (100% of {total_primary_reads:,} reads)",
                    "type": "wig",
                    "format": "bedgraph",
                    "url": bedgraph_uri,
                    "indexed": False,
                    "color": "#0f172a",
                    "height": 72,
                    "min": 0,
                }
            )
        pileup_label = (
            f"Diagnostic Read Pileup ({mapped_reads:,} stratified of {total_primary_reads:,} reads)"
            if downsampled
            else f"Diagnostic Read Pileup (all {mapped_reads:,} primary reads)"
        )
        tracks.append(
            {
                "name": pileup_label,
                "type": "alignment",
                "format": "bam",
                "url": bam_uri,
                "indexed": False,
                "showCoverage": False,
                "showSoftClips": True,
                "colorBy": "tag:PC",
                "tagColorTable": {
                    "Plasmid_ExonExon": "#dc2626",
                    "Vector_Insert": "#7c3aed",
                    "Human_ExonIntron": "#0284c7",
                    "Backbone": "#d97706",
                    "Shared_Exon": "#64748b",
                    "*": "#94a3b8",
                },
                "samplingDepth": 500,
                "height": 320,
            }
        )

        browser_options = {
            "loadDefaultGenomes": False,
            "reference": {
                "id": seq_name,
                "name": seq_name,
                "fastaURL": fasta_uri,
                "indexed": False,
            },
            "locus": loci[1]["locus"] if len(loci) > 1 else loci[0]["locus"],
            "tracks": tracks,
        }
        payload_js = (
            f"const IGV_GZ_B64 = {json.dumps(igv_gz_b64)};\n"
            f"const IGV_BROWSER_OPTIONS = {json.dumps(browser_options)};"
        )

        return {
            "igv_embedded": True,
            "igv_version": IGV_VERSION,
            "igv_sha256": IGV_SHA256,
            "igv_absence_reason": "",
            "igv_absence_detail": "",
            "igv_loci": loci,
            "igv_session_relpath": session_rel,
            "igv_downsampled": downsampled,
            "igv_read_count": mapped_reads,
            "igv_total_primary_reads": total_primary_reads,
            "igv_splice_junction_count": n_junc,
            "igv_human_genomic_locus": human_genomic_locus,
            "igv_payload": payload_js,
        }
    except Exception as exc:
        logging.warning("Skipping embedded IGV inspector: %s", exc)
        return base_ctx


def generate_report(
    summary_df: pd.DataFrame,
    output_folder: str,
    verdict: str,
    ratio: float,
    threshold: float,
    unclear: dict[str, float],
    command_line: str,
    human_fasta: str = "None",
    plasmid_gb: str = "None",
    sequencing_file: str = "None",
    boxplot_html: str = "",
    boxplot_filename_png: str | None = None,
    scatter_html: str = "",
    scatter_filename_png: str | None = None,
    downsampled: bool = False,
    static_report: bool = False,
    plotly_mode: str = "directory",
    plotly_version: str = "",
    plotly_js_path: str = "",
    plotly_js_inline: str = "",
    plasmid_count: int = 0,
    human_count: int = 0,
    tied_count: int = 0,
    backbone_only_count: int = 0,
    ambiguous_count: int = 0,
    total_reads: int = 0,
    plasmid_pct: str = "0.0",
    human_pct: str = "0.0",
    tied_pct: str = "0.0",
    backbone_only_pct: str = "0.0",
    ambiguous_pct: str = "0.0",
    coverage_outside_insert: str = "N/A",
    mismatches_near_insert: str = "N/A",
    verdict_short: str = "",
    verdict_class: str = "",
    gauge_position: float = 50.0,
    gauge_max_label: str = "10+",
    contam_cutoff: float | None = None,
    mismatches_with: int = 0,
    mismatches_without: int = 0,
    sample_name: str = "",
    plasmid_name: str = "",
    coverage_metrics: dict[str, str] | None = None,
    coverage_fallback: bool = False,
    reports_dir: str | None = None,
    junction_metrics: dict[str, str] | None = None,
    anonymize: bool = False,
    igv_context: dict[str, Any] | None = None,
) -> None:
    from jinja2 import Environment, FileSystemLoader

    if anonymize:

        def _shown(p: str) -> str:
            return os.path.basename(p.rstrip(os.sep)) if os.sep in p else p

        human_fasta = _shown(human_fasta)
        plasmid_gb = _shown(plasmid_gb)
        sequencing_file = _shown(sequencing_file)
        if command_line:
            command_line = " ".join(_shown(a) for a in command_line.split())

    resolved_igv = igv_context or _build_igv_context(
        output_folder, reports_dir=reports_dir, anonymize=anonymize
    )

    logging.info("Generating report")
    env = Environment(
        loader=FileSystemLoader(str(get_resource_path(TEMPLATE_DIR))), autoescape=True
    )
    env.filters["truncate_list"] = _truncate_list_filter
    template = env.get_template("report_template.html")

    logo_base64 = encode_image_to_base64(str(get_resource_path(LOGO_PATH)))

    # Prepare downsample message if data was downsampled
    downsample_message = (
        f"Plots show a random subset of {DOWNSAMPLE_LIMIT:,} reads for rendering performance."
        " The verdict and all statistics are computed from the complete dataset."
        if downsampled
        else ""
    )

    resolved_cutoff = (
        contam_cutoff if contam_cutoff is not None else unclear.get("upper_bound", threshold)
    )
    report_output_dir = reports_dir if reports_dir is not None else output_folder
    out_abs = os.path.abspath(report_output_dir)
    summary_report_relpath = ""
    for rel_cand in (
        "../../../summary/summary_report_interactive.html",
        "../../../summary/reports/summary_interactive.html",
        "../../summary/summary_report_interactive.html",
        "../../summary/reports/summary_interactive.html",
    ):
        if os.path.exists(os.path.normpath(os.path.join(out_abs, rel_cand))):
            summary_report_relpath = rel_cand
            break

    # Render the interactive HTML report
    html_content_interactive = template.render(
        box_plot=boxplot_html,
        scatter_plot=scatter_html,
        verdict=verdict,
        verdict_short=verdict_short,
        verdict_class=verdict_class,
        ratio=f"{ratio:.3f}",
        threshold=threshold,
        unclear_range=unclear,
        contam_cutoff=resolved_cutoff,
        version=VERSION,
        logo_base64=logo_base64,
        human_fasta=human_fasta,
        plasmid_gb=plasmid_gb,
        sequencing_file=sequencing_file,
        run_date=datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
        command_line=command_line,
        interactive=True,
        downsample_message=downsample_message,
        plotly_mode=plotly_mode,
        plotly_version=plotly_version,
        plotly_js_path=plotly_js_path,
        plotly_js_inline=plotly_js_inline,
        plasmid_count=plasmid_count,
        human_count=human_count,
        tied_count=tied_count,
        backbone_only_count=backbone_only_count,
        ambiguous_count=ambiguous_count,
        total_reads=total_reads,
        plasmid_pct=plasmid_pct,
        human_pct=human_pct,
        tied_pct=tied_pct,
        backbone_only_pct=backbone_only_pct,
        ambiguous_pct=ambiguous_pct,
        assignment_colors=ASSIGNMENT_COLORS,
        coverage_outside_insert=coverage_outside_insert,
        mismatches_with=mismatches_with,
        mismatches_without=mismatches_without,
        gauge_position=gauge_position,
        gauge_max_label=gauge_max_label,
        sample_name=sample_name,
        plasmid_name=plasmid_name,
        coverage_metrics=coverage_metrics or {},
        coverage_fallback=coverage_fallback,
        junction_metrics=junction_metrics or {},
        summary_report_relpath=summary_report_relpath,
        **resolved_igv,
    )

    # Save interactive HTML report
    html_report_interactive = os.path.join(report_output_dir, "report_interactive.html")
    with open(html_report_interactive, "w") as f:
        f.write(html_content_interactive)

    # Only generate non-interactive HTML report when static_report=True
    if static_report and boxplot_filename_png and scatter_filename_png:
        # Convert the static PNG files to Base64
        boxplot_png_base64 = encode_image_to_base64(boxplot_filename_png)
        scatter_png_base64 = encode_image_to_base64(scatter_filename_png)
        static_igv = {
            **resolved_igv,
            "igv_embedded": False,
            "igv_payload": "",
            "igv_absence_reason": "static_report",
            "igv_absence_detail": (
                "Interactive IGV alignment browser is only active in report_interactive.html."
            ),
        }

        # Render the non-interactive HTML report
        html_content_non_interactive = template.render(
            box_plot=boxplot_png_base64,
            scatter_plot=scatter_png_base64,
            verdict=verdict,
            verdict_short=verdict_short,
            verdict_class=verdict_class,
            ratio=f"{ratio:.3f}",
            threshold=threshold,
            unclear_range=unclear,
            contam_cutoff=resolved_cutoff,
            version=VERSION,
            logo_base64=logo_base64,
            human_fasta=human_fasta,
            plasmid_gb=plasmid_gb,
            sequencing_file=sequencing_file,
            run_date=datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M:%S"),
            command_line=command_line,
            interactive=False,
            downsample_message=downsample_message,
            plotly_mode=plotly_mode,
            plotly_version=plotly_version,
            plotly_js_path=plotly_js_path,
            plotly_js_inline=plotly_js_inline,
            plasmid_count=plasmid_count,
            human_count=human_count,
            tied_count=tied_count,
            backbone_only_count=backbone_only_count,
            ambiguous_count=ambiguous_count,
            total_reads=total_reads,
            plasmid_pct=plasmid_pct,
            human_pct=human_pct,
            tied_pct=tied_pct,
            backbone_only_pct=backbone_only_pct,
            ambiguous_pct=ambiguous_pct,
            assignment_colors=ASSIGNMENT_COLORS,
            coverage_outside_insert=coverage_outside_insert,
            mismatches_with=mismatches_with,
            mismatches_without=mismatches_without,
            gauge_position=gauge_position,
            gauge_max_label=gauge_max_label,
            sample_name=sample_name,
            plasmid_name=plasmid_name,
            coverage_metrics=coverage_metrics or {},
            coverage_fallback=coverage_fallback,
            junction_metrics=junction_metrics or {},
            **static_igv,
        )

        # Save non-interactive HTML report
        html_report_non_interactive = os.path.join(report_output_dir, "report_non_interactive.html")
        with open(html_report_non_interactive, "w") as f:
            f.write(html_content_non_interactive)


def main(
    reads_assignment_file: str,
    summary_file: str,
    output_folder: str,
    threshold: float = DEFAULT_THRESHOLD,
    unclear: dict[str, float] = UNCLEAR_RANGE,
    human_fasta: str = "None",
    plasmid_gb: str = "None",
    sequencing_file: str = "None",
    command_line: str = "",
    static_report: bool = False,
    plotly_mode: str = "directory",
    plot_backend: str = "matplotlib",
    reports_dir: str | None = None,
    anonymize: bool = False,
    sample_name: str | None = None,
) -> None:
    if anonymize:

        def _shown(p: str) -> str:
            return os.path.basename(p.rstrip(os.sep)) if os.sep in p else p

        human_fasta = _shown(human_fasta)
        plasmid_gb = _shown(plasmid_gb)
        sequencing_file = _shown(sequencing_file)
        if command_line:
            command_line = " ".join(_shown(a) for a in command_line.split())

    reads_df, summary_df = load_data(reads_assignment_file, summary_file)

    # Extract effective threshold from summary_df (Threshold row if present, else fallback)
    effective_threshold = threshold
    thresh_rows = summary_df[summary_df["Category"] == "Threshold"]
    if not thresh_rows.empty:
        try:
            val_col = "Count" if "Count" in thresh_rows.columns else "Value"
            effective_threshold = float(thresh_rows[val_col].iloc[0])
        except (ValueError, TypeError):
            effective_threshold = threshold

    unclear_lower, unclear_upper = verdict_bounds(effective_threshold, unclear)
    effective_unclear = {"lower_bound": unclear_lower, "upper_bound": unclear_upper}

    # Extract category counts from summary_df before downsampling
    def _get_count(df: pd.DataFrame, category: str) -> int:
        rows = df[df["Category"] == category]
        return int(rows["Count"].values[0]) if not rows.empty else 0

    plasmid_count = _get_count(summary_df, "Plasmid")
    human_count = _get_count(summary_df, "Human")
    tied_count = _get_count(summary_df, "Tied")
    backbone_only_count = _get_count(summary_df, "Backbone_Only")
    ambiguous_count = _get_count(summary_df, "Ambiguous")
    total_reads = plasmid_count + human_count + tied_count + backbone_only_count + ambiguous_count
    if total_reads == 0 and not reads_df.empty:
        total_reads = len(reads_df)

    # Downsample data if necessary
    reads_df, downsampled = downsample_data(reads_df, DOWNSAMPLE_LIMIT)

    # Generate plots (returns div HTML strings for embedding, writes standalone files)
    (
        boxplot_div_html,
        boxplot_filename_png,
        scatter_div_html,
        scatter_filename_png,
    ) = generate_plots(
        reads_df,
        output_folder,
        static_report=static_report,
        plot_backend=plot_backend,
        reports_dir=reports_dir,
        total_reads=total_reads,
        threshold=effective_threshold,
        unclear=effective_unclear,
    )

    # Calculate percentages
    def _pct(count: int, total: int) -> str:
        return f"{count / total * 100:.1f}" if total > 0 else "0.0"

    # Extract additional metrics
    cov_rows = summary_df[summary_df["Category"] == "CoverageOutsideINSERT"]
    coverage_outside_insert = str(cov_rows["Count"].values[0]) if not cov_rows.empty else "N/A"
    mis_rows = summary_df[summary_df["Category"] == "MismatchesNearINSERT"]
    mismatches_near_insert = str(mis_rows["Count"].values[0]) if not mis_rows.empty else "N/A"

    # Extract coverage metrics (Phase 9)
    def _get_metric(df: pd.DataFrame, category: str, default: str = "0.00") -> str:
        rows = df[df["Category"] == category]
        if rows.empty:
            return default
        value = rows["Count"].values[0]
        try:
            if value != value:  # NaN check (NaN != NaN is True)
                return default
        except (TypeError, ValueError):
            pass
        return str(value)

    coverage_metrics = {
        "mean_depth_insert": _get_metric(summary_df, "MeanDepthInsert"),
        "median_depth_insert": _get_metric(summary_df, "MedianDepthInsert"),
        "breadth_insert": _get_metric(summary_df, "BreadthInsert"),
        "breadth_insert_5x": _get_metric(summary_df, "BreadthInsert_5x"),
        "cv_insert": _get_metric(summary_df, "CoverageCV_Insert"),
        "mean_depth_backbone": _get_metric(summary_df, "MeanDepthBackbone"),
        "median_depth_backbone": _get_metric(summary_df, "MedianDepthBackbone"),
        "breadth_backbone": _get_metric(summary_df, "BreadthBackbone"),
        "breadth_backbone_5x": _get_metric(summary_df, "BreadthBackbone_5x"),
        "cv_backbone": _get_metric(summary_df, "CoverageCV_Backbone"),
    }

    # Extract junction metrics (Phase 15 & #121/#122)
    est_contam_raw = _get_metric(summary_df, "EstimatedContaminationFraction", default="N/A")
    est_contam_pct = "N/A"
    if est_contam_raw not in ("N/A", "NA", ""):
        try:
            est_contam_pct = f"{float(est_contam_raw) * 100:.1f}%"
        except (ValueError, TypeError):
            est_contam_pct = "N/A"

    junction_metrics = {
        "splice_junction_count": _get_metric(summary_df, "SpliceJunctionCount", default="N/A"),
        "exon_count": _get_metric(summary_df, "ExonCount", default="N/A"),
        "splice_junction_reads": _get_metric(summary_df, "SpliceJunctionReads", default="N/A"),
        "splice_junction_events": _get_metric(summary_df, "SpliceJunctionEvents", default="N/A"),
        "splice_junctions_covered": _get_metric(
            summary_df, "SpliceJunctionsCovered", default="N/A"
        ),
        "exon_intron_boundary_reads": _get_metric(
            summary_df, "ExonIntronBoundaryReads", default="N/A"
        ),
        "estimated_contamination_fraction": est_contam_raw,
        "estimated_contamination_pct": est_contam_pct,
        "vector_insert_junction_reads": _get_metric(
            summary_df, "VectorInsertJunctionReads", default="N/A"
        ),
        "vector_insert_bridge_pairs": _get_metric(
            summary_df, "VectorInsertBridgePairs", default="N/A"
        ),
        "exon_lengths": _get_metric(summary_df, "ExonLengths", default="N/A"),
        "intron_lengths": _get_metric(summary_df, "IntronLengths", default="N/A"),
    }

    # Detect coverage fallback mode
    coverage_fallback_rows = summary_df[summary_df["Category"] == "CoverageFallback"]
    coverage_fallback = (
        not coverage_fallback_rows.empty
        and str(coverage_fallback_rows["Count"].values[0]).lower() == "true"
    )

    # Determine if backbone filtering was unavailable
    ratio_str = _get_metric(summary_df, "Ratio", default="")
    has_stored_ratio = ratio_str != ""
    if has_stored_ratio:
        try:
            ratio = float(ratio_str)
        except (ValueError, TypeError):
            ratio = plasmid_count / human_count if human_count != 0 else float("inf")
    else:
        ratio = plasmid_count / human_count if human_count != 0 else float("inf")

    # Extract the verdict from the summary file or recompute from valid stored ratio
    verdict = extract_verdict_from_summary(summary_df)
    if has_stored_ratio and math.isfinite(ratio):
        verdict_cat = categorize_value(ratio, effective_threshold, effective_unclear)
        verdict = {
            "Contaminated": "Sample is contaminated with plasmid DNA",
            "Unclear": "Sample contamination status is unclear",
            "Clean": "Sample is not contaminated with plasmid DNA",
        }.get(verdict_cat, verdict)

    # Derive new template variables
    verdict_short, verdict_class = _derive_verdict_info(verdict)
    gauge_position = _compute_gauge_position(ratio, unclear_lower, unclear_upper)
    gauge_max_label = (
        f"{int(10.0 * unclear_upper)}+"
        if (10.0 * unclear_upper).is_integer()
        else f"{10.0 * unclear_upper:.1f}+"
    )
    mismatches_with, mismatches_without = _parse_mismatches(mismatches_near_insert)

    # Derive sample and plasmid display names from file paths unless explicitly provided
    sample_name = (
        sample_name
        if sample_name
        else (
            os.path.splitext(os.path.basename(sequencing_file))[0]
            if sequencing_file not in ("None", "")
            else ""
        )
    )
    plasmid_name = (
        os.path.splitext(os.path.basename(plasmid_gb))[0] if plasmid_gb not in ("None", "") else ""
    )

    # Place plotly assets in reports/assets/ so reports dir is self-contained.
    # Fall back to output_folder/assets/ for backward compatibility.
    report_html_dir = reports_dir if reports_dir else output_folder

    plotly_js_path = ""
    plotly_version = ""
    plotly_js_inline = ""

    if plotly_mode in ("directory", "cdn"):
        import plotly

        plotly_version = plotly.__version__
    if plotly_mode == "directory":
        _ensure_plotly_assets(report_html_dir)
        plotly_js_path = os.path.join("assets", "plotly.min.js")
    elif plotly_mode == "embedded":
        import plotly as _plotly

        plotly_js_source = os.path.join(
            os.path.dirname(_plotly.__file__), "package_data", "plotly.min.js"
        )
        with open(plotly_js_source) as f:
            plotly_js_inline = f"<script>{f.read()}</script>"

    # Generate report
    generate_report(
        summary_df,
        output_folder,
        verdict,
        ratio,
        effective_threshold,
        effective_unclear,
        command_line,
        human_fasta,
        plasmid_gb,
        sequencing_file,
        boxplot_html=boxplot_div_html,
        boxplot_filename_png=boxplot_filename_png,
        scatter_html=scatter_div_html,
        scatter_filename_png=scatter_filename_png,
        downsampled=downsampled,
        static_report=static_report,
        plotly_mode=plotly_mode,
        plotly_version=plotly_version,
        plotly_js_path=plotly_js_path,
        plotly_js_inline=plotly_js_inline,
        plasmid_count=plasmid_count,
        human_count=human_count,
        tied_count=tied_count,
        backbone_only_count=backbone_only_count,
        ambiguous_count=ambiguous_count,
        total_reads=total_reads,
        plasmid_pct=_pct(plasmid_count, total_reads),
        human_pct=_pct(human_count, total_reads),
        tied_pct=_pct(tied_count, total_reads),
        backbone_only_pct=_pct(backbone_only_count, total_reads),
        ambiguous_pct=_pct(ambiguous_count, total_reads),
        coverage_outside_insert=coverage_outside_insert,
        mismatches_near_insert=mismatches_near_insert,
        verdict_short=verdict_short,
        verdict_class=verdict_class,
        gauge_position=gauge_position,
        gauge_max_label=gauge_max_label,
        contam_cutoff=unclear_upper,
        mismatches_with=mismatches_with,
        mismatches_without=mismatches_without,
        sample_name=sample_name,
        plasmid_name=plasmid_name,
        coverage_metrics=coverage_metrics,
        coverage_fallback=coverage_fallback,
        reports_dir=reports_dir,
        junction_metrics=junction_metrics,
        anonymize=anonymize,
    )
