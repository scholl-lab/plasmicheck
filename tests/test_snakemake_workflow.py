"""Unit and acceptance tests for snakemake/plasmicheck.smk and cohort configuration (#129)."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

REPO_ROOT = Path(__file__).resolve().parent.parent
SNAKEMAKE_DIR = REPO_ROOT / "snakemake"
SMK_PATH = SNAKEMAKE_DIR / "plasmicheck.smk"
EXAMPLE_CFG_PATH = SNAKEMAKE_DIR / "config.example.yaml"
SCHEMA_PATH = SNAKEMAKE_DIR / "config.schema.yaml"
RUNNER_PATH = SNAKEMAKE_DIR / "run_plasmicheck.sh"
README_PATH = SNAKEMAKE_DIR / "README.md"


def _load_smk_namespace(
    config: dict[str, Any] | None = None,
    cwd: Path | None = None,
) -> dict[str, Any]:
    """Execute the Python preamble of snakemake/plasmicheck.smk before 'rule all:'."""
    source = SMK_PATH.read_text(encoding="utf-8")
    python_preamble = re.split(r"^rule\s+all\s*:", source, maxsplit=1, flags=re.MULTILINE)[0]

    def _expand(pattern: str, **wildcards: list[str]) -> list[str]:
        results = [pattern]
        for key, values in wildcards.items():
            placeholder = "{" + key + "}"
            results = [r.replace(placeholder, str(v)) for r in results for v in values]
        return results

    ns: dict[str, Any] = {
        "__file__": str(SMK_PATH),
        "__name__": "plasmicheck_smk",
        "expand": _expand,
    }
    if config is not None:
        ns["config"] = config

    prev_cwd = Path.cwd()
    try:
        if cwd is not None:
            os.chdir(cwd)
        exec(compile(python_preamble, str(SMK_PATH), "exec"), ns)
    finally:
        if cwd is not None:
            os.chdir(prev_cwd)
    return ns


@pytest.mark.unit
class TestConfigFilesAndSchema:
    def test_files_exist(self) -> None:
        assert SMK_PATH.is_file()
        assert EXAMPLE_CFG_PATH.is_file()
        assert SCHEMA_PATH.is_file()
        assert RUNNER_PATH.is_file()
        assert README_PATH.is_file()

    def test_example_config_matches_schema(self, tmp_path: Path) -> None:
        bams_txt = tmp_path / "data" / "bams.txt"
        plasmids_txt = tmp_path / "data" / "plasmids.txt"
        bams_txt.parent.mkdir(parents=True, exist_ok=True)
        bams_txt.write_text("data/sample1.bam\ndata/sample2.bam\n", encoding="utf-8")
        plasmids_txt.write_text("data/p1.gb\ndata/p2.gb\n", encoding="utf-8")

        cfg = yaml.safe_load(EXAMPLE_CFG_PATH.read_text(encoding="utf-8"))
        schema = yaml.safe_load(SCHEMA_PATH.read_text(encoding="utf-8"))

        assert isinstance(cfg, dict)
        assert schema["$schema"] == "https://json-schema.org/draft/2020-12/schema"
        assert schema["type"] == "object"

        schema_props = set(schema["properties"].keys())
        assert set(cfg.keys()).issubset(schema_props)

        rule_res_schema = schema["$defs"]["rule_resources"]["properties"]
        assert rule_res_schema["runtime"]["maximum"] == 120

        ns = _load_smk_namespace(config=cfg, cwd=tmp_path)
        assert ns["SAMPLES"] == ["sample1", "sample2"]
        assert ns["get_rule_resource"](cfg, "scan_sample", "threads") == 4
        assert ns["get_rule_resource"](cfg, "scan_sample", "mem_mb") == 8000
        assert ns["get_rule_resource"](cfg, "scan_sample", "runtime") == 60

    def test_validate_config_rejects_invalid_values(self, tmp_path: Path) -> None:
        bams = tmp_path / "bams.txt"
        bams.write_text("a.bam\n", encoding="utf-8")
        ns = _load_smk_namespace(config={"samples": {"s1": "a.bam"}})
        validate_config = ns["validate_config"]

        with pytest.raises(ValueError, match="dictionary"):
            validate_config("not_a_dict")

        with pytest.raises(ValueError, match="'samples' must not be empty"):
            validate_config({"samples": {}})

        with pytest.raises(ValueError, match="'samples' must be a dict"):
            validate_config({"samples": 123})

        with pytest.raises(ValueError, match="'constructs' must not be empty"):
            validate_config({"constructs": []})

        with pytest.raises(ValueError, match="'constructs' must be a list"):
            validate_config({"constructs": 123})

        with pytest.raises(ValueError, match="'threshold' must be a positive number"):
            validate_config({"threshold": 0})

        with pytest.raises(ValueError, match="'threshold' must be a positive number"):
            validate_config({"threshold": True})

        with pytest.raises(ValueError, match="'padding' must be a non-negative integer"):
            validate_config({"padding": -1})

        with pytest.raises(ValueError, match="Unknown rule in 'resources'"):
            validate_config({"resources": {"bad_rule": {"threads": 2}}})

        with pytest.raises(ValueError, match="Unknown resource key"):
            validate_config({"resources": {"scan_sample": {"cpus": 2}}})

        with pytest.raises(ValueError, match="exceeds the 120 min"):
            validate_config({"resources": {"scan_sample": {"runtime": 180}}})


@pytest.mark.unit
class TestResolveSamples:
    def test_legacy_fallback_reads_data_bams_txt(self, tmp_path: Path) -> None:
        data_dir = tmp_path / "data"
        data_dir.mkdir()
        (data_dir / "bams.txt").write_text(
            "/path/to/APA1-N.bam\n/path/to/APA1-T.bam\n", encoding="utf-8"
        )
        (data_dir / "plasmids.txt").write_text("data/p1.gb\n", encoding="utf-8")

        ns = _load_smk_namespace(config=None, cwd=tmp_path)
        assert ns["SAMPLE_BAMS"] == {
            "APA1-N": "/path/to/APA1-N.bam",
            "APA1-T": "/path/to/APA1-T.bam",
        }
        assert ns["PLASMIDS_FILE"] == "data/plasmids.txt"
        assert ns["CONSTRUCT_FILES"] == ["data/p1.gb"]

    def test_samples_from_dict_and_list(self) -> None:
        ns = _load_smk_namespace(config={"samples": {"S1": "/bams/s1.bam", "S2": "/bams/s2.bam"}})
        resolve_samples = ns["resolve_samples"]
        assert resolve_samples({"samples": {"S1": "/bams/s1.bam"}}) == {"S1": "/bams/s1.bam"}
        assert resolve_samples({"samples": ["/bams/alpha.bam", "/bams/beta.bam"]}) == {
            "alpha": "/bams/alpha.bam",
            "beta": "/bams/beta.bam",
        }

    def test_samples_from_tsv_with_header_and_headerless(self, tmp_path: Path) -> None:
        ns = _load_smk_namespace(config={"samples": {"S1": "s1.bam"}})
        resolve_samples = ns["resolve_samples"]

        tsv_header = tmp_path / "samples_header.tsv"
        tsv_header.write_text(
            "sample\tbam\nAPA1-N\t/bams/APA1-N.bam\nAPA1-T\t/bams/APA1-T.bam\n",
            encoding="utf-8",
        )
        assert resolve_samples({"samples": str(tsv_header)}) == {
            "APA1-N": "/bams/APA1-N.bam",
            "APA1-T": "/bams/APA1-T.bam",
        }

        tsv_no_header = tmp_path / "samples_no_header.tsv"
        tsv_no_header.write_text(
            "APA2-N\t/bams/APA2-N.bam\nAPA2-T\t/bams/APA2-T.bam\n",
            encoding="utf-8",
        )
        assert resolve_samples({"samples": str(tsv_no_header)}) == {
            "APA2-N": "/bams/APA2-N.bam",
            "APA2-T": "/bams/APA2-T.bam",
        }

    def test_samples_from_apa_sample_map_tsv(self, tmp_path: Path) -> None:
        ns = _load_smk_namespace(config={"samples": {"S1": "s1.bam"}})
        resolve_samples = ns["resolve_samples"]

        sample_map = tmp_path / "sample_map.tsv"
        sample_map.write_text(
            "library_id\toriginal_id\tcorrected_id\trole\n"
            "lib1\tAPA12-T\tAPA13-T\ttumour\n"
            "lib2\tAPA13-T\tAPA12-T\ttumour\n",
            encoding="utf-8",
        )
        resolved = resolve_samples(
            {
                "samples": str(sample_map),
                "bam_folder": "/results/align/bqsr",
                "bam_ext": ".merged.dedup.bqsr.bam",
            }
        )
        assert resolved == {
            "APA13-T": "/results/align/bqsr/APA13-T.merged.dedup.bqsr.bam",
            "APA12-T": "/results/align/bqsr/APA12-T.merged.dedup.bqsr.bam",
        }

    def test_duplicate_samples_raise(self, tmp_path: Path) -> None:
        ns = _load_smk_namespace(config={"samples": {"S1": "s1.bam"}})
        resolve_samples = ns["resolve_samples"]

        with pytest.raises(ValueError, match="Duplicate sample key"):
            resolve_samples({"samples": ["/dir1/dup.bam", "/dir2/dup.bam"]})

        bad_tsv = tmp_path / "dup.tsv"
        bad_tsv.write_text("sample\tbam\nS1\t/a.bam\nS1\t/b.bam\n", encoding="utf-8")
        with pytest.raises(ValueError, match="Duplicate sample key"):
            resolve_samples({"samples": str(bad_tsv)})


@pytest.mark.unit
class TestResolvePlasmidsAndFlags:
    def test_resolve_plasmids_from_constructs_list(self, tmp_path: Path) -> None:
        out_dir = tmp_path / "out"
        constructs = ["data/pKCNJ5_WT.gb", "data/pCACNA1H.gb", "data/pCLCN2.gb"]
        ns = _load_smk_namespace(
            config={
                "samples": {"S1": "s1.bam"},
                "constructs": constructs,
                "output_dir": str(out_dir),
            }
        )
        expected_pf = out_dir / ".plasmids.txt"
        assert Path(ns["PLASMIDS_FILE"]) == expected_pf
        assert expected_pf.read_text(encoding="utf-8").splitlines() == constructs
        assert ns["CONSTRUCT_FILES"] == constructs

    def test_build_pipeline_and_summary_flags(self) -> None:
        ns = _load_smk_namespace(config={"samples": {"S1": "s1.bam"}})
        build_pipeline_flags = ns["build_pipeline_flags"]
        build_summary_flags = ns["build_summary_flags"]
        get_mem_from_threads = ns["get_mem_from_threads"]

        assert build_pipeline_flags({}) == ""
        assert build_summary_flags({}) == ""
        assert get_mem_from_threads(None, 8) == 17600

        flags_default_unmapped = build_pipeline_flags({"regions_bed": "construct_loci.bed"})
        assert "--regions construct_loci.bed" in flags_default_unmapped
        assert "--include-unmapped" in flags_default_unmapped

        flags_no_unmapped = build_pipeline_flags(
            {"regions_bed": "construct_loci.bed", "include_unmapped": False}
        )
        assert "--regions construct_loci.bed" in flags_no_unmapped
        assert "--include-unmapped" not in flags_no_unmapped

        full_flags = build_pipeline_flags(
            {
                "cache_dir": "/tmp/pc_cache",
                "regions_bed": "loci.bed",
                "include_unmapped": True,
                "target_bed": "exome_targets.bed",
                "gtf": "MANE.gtf.gz",
                "vcf": "variants.vcf.gz",
                "threshold": 0.75,
                "padding": 500,
                "static_report": True,
                "archive_output": True,
                "overwrite": True,
            }
        )
        for expected in (
            "--cache-dir /tmp/pc_cache",
            "--regions loci.bed",
            "--include-unmapped",
            "--target-bed exome_targets.bed",
            "--gtf MANE.gtf.gz",
            "--vcf variants.vcf.gz",
            "--threshold 0.75",
            "--padding 500",
            "--static-report",
            "--archive_output",
            "--overwrite",
        ):
            assert expected in full_flags

        summary_flags = build_summary_flags(
            {
                "threshold": 0.75,
                "static_report": True,
                "substring_to_remove": ".merged.dedup.bqsr",
            }
        )
        assert "--threshold 0.75" in summary_flags
        assert "--static-report" in summary_flags
        assert "--substring_to_remove .merged.dedup.bqsr" in summary_flags


@pytest.mark.unit
class TestCohortDag97x21:
    """Verify 97 exome BAMs x 21 constructs produces 97 per-sample scan jobs (no cross-product) and <2h walltimes."""

    def test_97_bams_by_21_constructs_dag_and_resources(self, tmp_path: Path) -> None:
        samples_tsv = tmp_path / "samples_97.tsv"
        rows = ["sample\tbam"] + [
            f"APA{i:02d}\t/data/bqsr/APA{i:02d}.merged.dedup.bqsr.bam" for i in range(1, 98)
        ]
        samples_tsv.write_text("\n".join(rows) + "\n", encoding="utf-8")

        construct_names = [
            "pcDNA5_FRT_TO_KCNJ5_WT.gb",
            "pcDNA5_FRT_TO_KCNJ5_G151R.gb",
            "pcDNA5_FRT_TO_KCNJ5_L168R.gb",
            "pCMV6_CACNA1H.gb",
            "pcDNA5_CACNA1H_S1073C.gb",
            "pcDNA5_CLCN2.gb",
            *[f"construct_{idx:02d}.gb" for idx in range(7, 22)],
        ]
        assert len(construct_names) == 21
        plasmids_txt = tmp_path / "plasmids_21.txt"
        plasmids_txt.write_text(
            "\n".join(f"data/plasmids/{c}" for c in construct_names) + "\n",
            encoding="utf-8",
        )

        cfg: dict[str, Any] = {
            "human_fasta": "reference/GRCh38.fa",
            "samples": str(samples_tsv),
            "plasmids_file": str(plasmids_txt),
            "output_dir": str(tmp_path / "output"),
            "summary_dir": str(tmp_path / "output" / "summary"),
            "log_dir": str(tmp_path / "logs"),
            "cache_dir": str(tmp_path / "cache"),
            "regions_bed": "config/construct_regions.bed",
            "include_unmapped": True,
            "target_bed": "config/targets_union_pad100.bed",
            "gtf": "config/MANE.GRCh38.v1.4.gtf.gz",
            "vcf": "results/hotspot_rescue.tsv",
            "resources": {
                "prepare_constructs": {"threads": 4, "mem_mb": 16000, "runtime": 60},
                "scan_sample": {"threads": 4, "mem_mb": 8000, "runtime": 60},
                "summary_reports": {"threads": 2, "mem_mb": 4000, "runtime": 30},
            },
        }

        ns = _load_smk_namespace(config=cfg)

        # 97 samples and 21 constructs resolved without N_samples x N_plasmids cross-product
        assert len(ns["SAMPLES"]) == 97
        assert len(ns["CONSTRUCT_FILES"]) == 21
        assert "jobs" not in ns

        done_targets = ns["expand"](f"{ns['LOG_DIR']}/{{sample}}.done", sample=ns["SAMPLES"])
        assert len(done_targets) == 97  # 1 job per sample BAM, NOT 97 * 21 = 2,037

        # Verify rule definitions in snakemake/plasmicheck.smk
        smk_text = SMK_PATH.read_text(encoding="utf-8")
        defined_rules = re.findall(r"^rule\s+(\w+)\s*:", smk_text, flags=re.MULTILINE)
        assert defined_rules == ["all", "prepare_constructs", "scan_sample", "summary_reports"]
        assert "itertools.product" not in smk_text
        assert "-pf {input.plasmids_file:q} -sf1 {input.bam:q}" in smk_text or (
            "-pf {input.plasmids_file:q}" in smk_text and "-sf1 {input.bam:q}" in smk_text
        )

        # Verify all per-rule walltimes are < 2 hours (120 min)
        for rule_name in ("prepare_constructs", "scan_sample", "summary_reports"):
            runtime = ns["get_rule_resource"](cfg, rule_name, "runtime")
            threads = ns["get_rule_resource"](cfg, rule_name, "threads")
            mem_mb = ns["get_rule_resource"](cfg, rule_name, "mem_mb")
            assert 0 < runtime <= 120, f"{rule_name} runtime {runtime} exceeds 120 min"
            assert threads >= 1
            assert mem_mb >= 1000
