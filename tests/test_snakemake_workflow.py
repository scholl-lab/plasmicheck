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


def _validate_against_schema(
    value: Any, node: dict[str, Any], root: dict[str, Any] | None = None
) -> None:
    """Validate a config value against the JSON Schema Draft 2020-12 rules in config.schema.yaml."""
    if root is None:
        root = node
    if "$ref" in node:
        ref_name = str(node["$ref"]).removeprefix("#/$defs/")
        node = root["$defs"][ref_name]
    if "oneOf" in node:
        for branch in node["oneOf"]:
            try:
                _validate_against_schema(value, branch, root)
                return
            except ValueError:
                continue
        raise ValueError(f"Value {value!r} did not match any oneOf branch in {node!r}")

    allowed_types = node.get("type")
    if allowed_types is not None:
        types = [allowed_types] if isinstance(allowed_types, str) else list(allowed_types)
        type_ok = False
        for t in types:
            if (
                (t == "null" and value is None)
                or (t == "boolean" and isinstance(value, bool))
                or (t == "integer" and isinstance(value, int) and not isinstance(value, bool))
                or (
                    t == "number"
                    and isinstance(value, (int, float))
                    and not isinstance(value, bool)
                )
                or (t == "string" and isinstance(value, str))
                or (t == "array" and isinstance(value, list))
                or (t == "object" and isinstance(value, dict))
            ):
                type_ok = True
        if not type_ok:
            raise ValueError(f"Expected type {types} for {value!r}")

    if isinstance(value, str) and "minLength" in node and len(value) < int(node["minLength"]):
        raise ValueError(f"String {value!r} shorter than minLength {node['minLength']}")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in node and value < node["minimum"]:
            raise ValueError(f"Value {value} < minimum {node['minimum']}")
        if "exclusiveMinimum" in node and value <= node["exclusiveMinimum"]:
            raise ValueError(f"Value {value} <= exclusiveMinimum {node['exclusiveMinimum']}")
        if "maximum" in node and value > node["maximum"]:
            raise ValueError(f"Value {value} > maximum {node['maximum']}")
    if isinstance(value, list):
        if "minItems" in node and len(value) < int(node["minItems"]):
            raise ValueError(f"Array shorter than minItems {node['minItems']}")
        if "items" in node:
            for item in value:
                _validate_against_schema(item, node["items"], root)
    if isinstance(value, dict):
        if "minProperties" in node and len(value) < int(node["minProperties"]):
            raise ValueError(f"Object has fewer than minProperties {node['minProperties']}")
        props = node.get("properties", {})
        addl = node.get("additionalProperties", True)
        for k, v in value.items():
            if k in props:
                _validate_against_schema(v, props[k], root)
            elif addl is False:
                raise ValueError(f"Additional property not allowed: {k!r}")
            elif isinstance(addl, dict):
                _validate_against_schema(v, addl, root)


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
        _validate_against_schema(cfg, schema)

        rule_res_schema = schema["$defs"]["rule_resources"]["properties"]
        assert rule_res_schema["runtime"]["maximum"] == 120

        ns = _load_smk_namespace(config=cfg, cwd=tmp_path)
        assert ns["SAMPLES"] == ["sample1", "sample2"]
        assert ns["get_rule_resource"](cfg, "prepare_constructs", "runtime") == 30
        assert ns["get_rule_resource"](cfg, "scan_sample", "threads") == 4
        assert ns["get_rule_resource"](cfg, "scan_sample", "mem_mb") == 4000
        assert ns["get_rule_resource"](cfg, "scan_sample", "runtime") == 20

    def test_v042_keys_and_aliases_validate_against_schema(self) -> None:
        schema = yaml.safe_load(SCHEMA_PATH.read_text(encoding="utf-8"))
        props = schema["properties"]

        for key in ("target_bed", "gtf_file", "gtf", "mane_gtf", "vcf_file", "vcf"):
            assert key in props, f"Missing {key!r} in config.schema.yaml"
            assert props[key]["type"] == ["string", "null"]
            assert props[key]["default"] is None
            assert isinstance(props[key].get("description"), str) and props[key]["description"]

        assert props["shift_bases"]["type"] == ["integer", "null"]
        assert props["shift_bases"]["default"] == 500
        assert props["generate_shifted"]["type"] == "boolean"
        assert props["generate_shifted"]["default"] is False

        canonical_cfg = {
            "human_fasta": "reference/GRCh38.fa",
            "samples": {"S1": "data/s1.bam"},
            "constructs": ["data/p1.gb"],
            "target_bed": "config/exome_targets.bed",
            "gtf_file": "reference/MANE.GRCh38.v1.3.refseq_genomic.gtf.gz",
            "vcf_file": "results/variants/cohort.vcf.gz",
            "shift_bases": 500,
            "generate_shifted": False,
        }
        _validate_against_schema(canonical_cfg, schema)

        for alias_key, alias_val in (
            ("gtf", "reference/genes.gtf"),
            ("mane_gtf", "reference/mane.gtf.gz"),
            ("vcf", "results/variants/calls.bcf"),
        ):
            alias_cfg = {
                "samples": {"S1": "data/s1.bam"},
                "constructs": ["data/p1.gb"],
                alias_key: alias_val,
            }
            _validate_against_schema(alias_cfg, schema)

        with pytest.raises(ValueError, match="Additional property not allowed"):
            _validate_against_schema({"samples": {"S1": "s1.bam"}, "undeclared_key": "x"}, schema)

        with pytest.raises(ValueError, match="Expected type"):
            _validate_against_schema({"samples": {"S1": "s1.bam"}, "target_bed": 123}, schema)

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

        with pytest.raises(ValueError, match="'shift_bases' must be a non-negative integer"):
            validate_config({"shift_bases": -1})

        with pytest.raises(ValueError, match="'shift_bases' must be a non-negative integer"):
            validate_config({"shift_bases": True})

        with pytest.raises(ValueError, match="'generate_shifted' must be a boolean"):
            validate_config({"generate_shifted": "yes"})

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
        build_prepare_kwargs = ns["build_prepare_kwargs"]

        assert build_pipeline_flags({}) == ""
        assert build_summary_flags({}) == ""
        assert build_prepare_kwargs(
            {
                "extra_args": "--gtf=/path/mane.gtf --padding=1500 --shift_bases=250 --cache-dir=/tmp/c -g -w"
            }
        ) == {
            "cache_dir": "/tmp/c",
            "padding": 1500,
            "shift_bases": 250,
            "generate_shifted": True,
            "overwrite": True,
            "gtf_file": "/path/mane.gtf",
        }

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
                "gtf_file": "mane.gtf.gz",
                "vcf_file": "cohort.vcf.gz",
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
            "--gtf mane.gtf.gz",
            "--vcf cohort.vcf.gz",
            "--threshold 0.75",
            "--padding 500",
            "--static-report",
            "--archive_output",
            "--overwrite",
        ):
            assert expected in full_flags

        assert "--gtf alias.gtf" in build_pipeline_flags({"gtf": "alias.gtf"})
        assert "--gtf mane_alias.gtf" in build_pipeline_flags({"mane_gtf": "mane_alias.gtf"})
        assert "--vcf alias.vcf.gz" in build_pipeline_flags({"vcf": "alias.vcf.gz"})

        summary_flags = build_summary_flags(
            {
                "threshold": 0.75,
                "static_report": True,
                "substring_to_remove": ".merged.dedup.bqsr",
                "vcf_file": "cohort.vcf.gz",
            }
        )
        assert "--threshold 0.75" in summary_flags
        assert "--static-report" in summary_flags
        assert "--substring_to_remove .merged.dedup.bqsr" in summary_flags
        assert "--vcf cohort.vcf.gz" in summary_flags
        assert "--vcf alias.vcf.gz" in build_summary_flags({"vcf": "alias.vcf.gz"})

    def test_generated_pipeline_flags_execute_in_cli(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Verify flags emitted by build_pipeline_flags and build_summary_flags are accepted and forwarded by plasmicheck.cli.main."""
        import shlex
        import sys
        from unittest.mock import patch

        from plasmicheck import cli

        ns = _load_smk_namespace(config={"samples": {"S1": "s1.bam"}})
        bed_file = tmp_path / "construct_loci.bed"
        bed_file.write_text("chr1\t100\t500\n", encoding="utf-8")
        target_bed = tmp_path / "exome_targets.bed"
        target_bed.write_text("chr1\t0\t1000\n", encoding="utf-8")
        gtf_file = tmp_path / "mane.gtf.gz"
        vcf_file = tmp_path / "cohort.vcf.gz"

        flags_str = ns["build_pipeline_flags"](
            {
                "cache_dir": str(tmp_path / "cache"),
                "regions_bed": str(bed_file),
                "include_unmapped": True,
                "target_bed": str(target_bed),
                "gtf_file": str(gtf_file),
                "vcf_file": str(vcf_file),
                "threshold": 0.75,
                "padding": 500,
                "static_report": True,
                "archive_output": True,
                "overwrite": True,
            }
        )
        argv = [
            "plasmicheck",
            "pipeline",
            "-hf",
            "human.fa",
            "-pf",
            "plasmids.txt",
            "-sf1",
            "sample.bam",
            "--sample-name",
            "APA01",
            "-o",
            str(tmp_path / "out"),
            *shlex.split(flags_str),
        ]
        monkeypatch.setattr(sys, "argv", argv)
        with (
            patch.object(cli, "setup_logging"),
            patch.object(cli, "check_requirements"),
            patch("plasmicheck.scripts.run_pipeline.run_pipeline") as mock_run_pipeline,
        ):
            cli.main()

        mock_run_pipeline.assert_called_once()
        kwargs = mock_run_pipeline.call_args.kwargs
        assert kwargs.get("regions") == str(bed_file) or kwargs.get("regions_bed") == str(bed_file)
        assert kwargs["include_unmapped"] is True
        assert kwargs["target_bed"] == str(target_bed)
        assert kwargs["gtf_file"] == str(gtf_file)
        assert kwargs["vcf_file"] == str(vcf_file)
        assert kwargs["sample_name"] == "APA01"
        assert kwargs["cache_dir"] == str(tmp_path / "cache")

        summary_flags_str = ns["build_summary_flags"]({"vcf_file": str(vcf_file)})
        summary_argv = [
            "plasmicheck",
            "summary_reports",
            "-i",
            str(tmp_path / "manifest.txt"),
            "-o",
            str(tmp_path / "summary"),
            *shlex.split(summary_flags_str),
        ]
        monkeypatch.setattr(sys, "argv", summary_argv)
        with (
            patch.object(cli, "setup_logging"),
            patch.object(cli, "check_requirements"),
            patch("plasmicheck.scripts.generate_summary_reports.main") as mock_summary_reports,
        ):
            cli.main()

        mock_summary_reports.assert_called_once()
        assert mock_summary_reports.call_args.kwargs["vcf_file"] == str(vcf_file)

    def test_bed_restricted_bam_extraction_command(self, tmp_path: Path) -> None:
        """Verify _bam_to_fastq_cmd restricts BAM extraction to BED regions and optional unmapped reads."""
        from plasmicheck.scripts.align_reads import _bam_to_fastq_cmd

        bam = tmp_path / "sample.bam"
        bam.write_bytes(b"BAM\x01")
        (tmp_path / "sample.bam.bai").write_bytes(b"BAI\x01")
        bed = tmp_path / "regions.bed"
        bed.write_text("chr11\t128890000\t128940000\n", encoding="utf-8")

        cmd_with_unmapped = _bam_to_fastq_cmd(
            str(bam), 2, regions_bed=str(bed), include_unmapped=True
        )
        assert "chr11:128890001-128940000" in cmd_with_unmapped
        assert "'*'" in cmd_with_unmapped or '"*"' in cmd_with_unmapped

        cmd_regions_only = _bam_to_fastq_cmd(
            str(bam), 2, regions_bed=str(bed), include_unmapped=False
        )
        assert "chr11:128890001-128940000" in cmd_regions_only
        assert "'*'" not in cmd_regions_only and '"*"' not in cmd_regions_only

        with pytest.raises(FileNotFoundError, match="Regions BED file not found"):
            _bam_to_fastq_cmd(str(bam), 2, regions_bed=str(tmp_path / "missing.bed"))


@pytest.mark.unit
class TestCohortExecutionAndGuards:
    def test_prepare_references_populates_cache_and_first_scan_reuses_it(
        self, tmp_path: Path
    ) -> None:
        """Verify prepare_references registers indexes and spliced artifacts in CacheManager and first scan hits cache."""
        from unittest.mock import patch

        from plasmicheck.cache_manager import CacheManager
        from plasmicheck.scripts.run_pipeline import prepare_references, run_pipeline

        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\n" + "ACGT" * 50 + "\n", encoding="utf-8")
        plasmid_gb = tmp_path / "pKCNJ5.gb"
        plasmid_gb.write_text(
            "LOCUS       pKCNJ5 40 bp DNA circular\nORIGIN\n 1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n//\n",
            encoding="utf-8",
        )
        plasmids_txt = tmp_path / "plasmids.txt"
        plasmids_txt.write_text(f"{plasmid_gb}\n", encoding="utf-8")
        cache_dir = tmp_path / "cache"

        def fake_convert(infile: str, outfile: str, *args: object, **kwargs: object) -> None:
            Path(outfile).write_text(">pKCNJ5\n" + "ACGT" * 10 + "\n", encoding="utf-8")

        def fake_create_indexes(fasta: str, *args: object, **kwargs: object) -> None:
            base = os.path.splitext(fasta)[0]
            Path(f"{base}.mmi").write_bytes(b"MMI_INDEX")
            Path(f"{fasta}.fai").write_text("seq\t40\t6\t40\t41\n", encoding="utf-8")

        def fake_spliced_alignment(
            human_idx: str, plasmid_fa: str, out_bam: str, *args: object, **kwargs: object
        ) -> list[tuple[str, int, int]]:
            Path(out_bam).write_bytes(b"BAM_DATA")
            return [("chr1", 0, 100)]

        def fake_extract_human(
            human_fa: str, regions: object, out_fa: str, *args: object, **kwargs: object
        ) -> None:
            Path(out_fa).write_text(">chr1_slice\n" + "ACGT" * 25 + "\n", encoding="utf-8")

        def fake_extract_cdna(
            plasmid_fa: str, spliced_bam: str, out_txt: str, *args: object, **kwargs: object
        ) -> None:
            assert kwargs.get("padding") == 750
            assert kwargs.get("gtf_file") is None
            Path(out_txt).write_text(
                "cDNA start position in plasmid: 1\ncDNA end position in plasmid: 40\n",
                encoding="utf-8",
            )

        with (
            patch("plasmicheck.scripts.run_pipeline.quality_control"),
            patch(
                "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
                side_effect=fake_convert,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.create_indexes",
                side_effect=fake_create_indexes,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.spliced_alignment",
                side_effect=fake_spliced_alignment,
            ) as mock_spliced_prep,
            patch(
                "plasmicheck.scripts.run_pipeline.extract_human_reference",
                side_effect=fake_extract_human,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions",
                side_effect=fake_extract_cdna,
            ),
        ):
            prepare_references(
                str(human_fasta),
                str(plasmids_txt),
                cache_dir=str(cache_dir),
                padding=750,
                threads=2,
            )
            # Re-running with overwrite=True refreshes cached construct/spliced entries
            prepare_references(
                str(human_fasta),
                str(plasmids_txt),
                cache_dir=str(cache_dir),
                padding=750,
                overwrite=True,
                threads=2,
            )

        assert mock_spliced_prep.call_count == 2
        cache = CacheManager(cache_dir)
        info = cache.info()
        # Expect plasmid index + spliced alignment + spliced reference index = 3 entries (no human.mmi copy)
        assert info["total_entries"] == 3

        # Now run first sample scan and verify cached plasmid/spliced artifacts are reused without recomputation
        seq_fq = tmp_path / "sample.fastq"
        seq_fq.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
        out_dir = tmp_path / "output"

        with (
            patch("plasmicheck.scripts.run_pipeline.quality_control"),
            patch(
                "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
                side_effect=fake_convert,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.create_indexes",
                side_effect=fake_create_indexes,
            ) as mock_idx_scan,
            patch("plasmicheck.scripts.run_pipeline.spliced_alignment") as mock_spliced_scan,
            patch("plasmicheck.scripts.run_pipeline.extract_human_reference") as mock_ext_scan,
            patch(
                "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"
            ) as mock_cdna_scan,
            patch("plasmicheck.scripts.run_pipeline.pysam.index"),
            patch("plasmicheck.scripts.run_pipeline.align_reads"),
            patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            patch("plasmicheck.scripts.run_pipeline.generate_report"),
        ):
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmids_txt),
                output_folder=str(out_dir),
                sequencing_files_r1=str(seq_fq),
                cache_dir=str(cache_dir),
                padding=750,
                md5_level="none",
            )

        mock_spliced_scan.assert_not_called()
        mock_ext_scan.assert_not_called()
        mock_cdna_scan.assert_not_called()
        # create_indexes is only called for the upfront human_fasta check, never for plasmid or spliced FASTA
        assert [c.args[0] for c in mock_idx_scan.call_args_list] == [str(human_fasta)]

    def test_prepare_references_with_gtf_matches_scan_sample_cache_key(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Verify rule prepare_constructs writes the exact spliced cache key that scan_sample CLI (build_pipeline_flags -> cli.main) reads (#133)."""
        import logging
        import shlex
        import sys
        from unittest.mock import patch

        from plasmicheck import cli
        from plasmicheck.cache_manager import CacheManager

        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\n" + "ACGT" * 50 + "\n", encoding="utf-8")
        plasmid_gb = tmp_path / "pKCNJ5.gb"
        plasmid_gb.write_text(
            "LOCUS       pKCNJ5 40 bp DNA circular\nORIGIN\n 1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n//\n",
            encoding="utf-8",
        )
        plasmids_txt = tmp_path / "plasmids.txt"
        plasmids_txt.write_text(f"{plasmid_gb}\n", encoding="utf-8")
        gtf_file = tmp_path / "mane.gtf"
        gtf_file.write_text(
            'chr1\tMANE\texon\t1\t120\t.\t+\t.\tgene_id "G1"; transcript_id "TX1"; tag "MANE_Select";\n',
            encoding="utf-8",
        )
        cache_dir = tmp_path / "cache"

        cfg = {
            "human_fasta": str(human_fasta),
            "samples": {"S1": "s1.bam"},
            "plasmids_file": str(plasmids_txt),
            "cache_dir": str(cache_dir),
            "gtf_file": str(gtf_file),
            "padding": 1500,
            "shift_bases": 250,
            "generate_shifted": True,
        }
        ns = _load_smk_namespace(config=cfg)
        prep_kwargs = ns["PREPARE_KWARGS"]
        assert prep_kwargs == {
            "cache_dir": str(cache_dir),
            "padding": 1500,
            "shift_bases": 250,
            "generate_shifted": True,
            "overwrite": False,
            "gtf_file": str(gtf_file),
        }

        # Extract the exact python3 -c one-liner from rule prepare_constructs in snakemake/plasmicheck.smk
        smk_text = SMK_PATH.read_text(encoding="utf-8")
        m_one_liner = re.search(r"python3 -c '([^']+)'", smk_text)
        assert m_one_liner is not None
        prep_one_liner = m_one_liner.group(1)

        def fake_convert(
            infile: str,
            outfile: str,
            file_type: str,
            shift_bases: int = 500,
            generate_shifted: bool = False,
            *args: object,
            **kwargs: object,
        ) -> None:
            Path(outfile).write_text(
                f">pKCNJ5_{shift_bases}_{int(generate_shifted)}\n" + "ACGT" * 10 + "\n",
                encoding="utf-8",
            )

        def fake_create_indexes(fasta: str, *args: object, **kwargs: object) -> None:
            base = os.path.splitext(fasta)[0]
            Path(f"{base}.mmi").write_bytes(b"MMI_INDEX")
            Path(f"{fasta}.fai").write_text("seq\t40\t6\t40\t41\n", encoding="utf-8")

        def fake_spliced_alignment(
            human_idx: str, plasmid_fa: str, out_bam: str, *args: object, **kwargs: object
        ) -> list[tuple[str, int, int]]:
            Path(out_bam).write_bytes(b"BAM_DATA")
            return [("chr1", 0, 100)]

        def fake_extract_human(
            human_fa: str, regions: object, out_fa: str, *args: object, **kwargs: object
        ) -> None:
            Path(out_fa).write_text(">chr1:0-100\n" + "ACGT" * 25 + "\n", encoding="utf-8")

        def fake_extract_cdna(
            plasmid_fa: str, spliced_bam: str, out_txt: str, *args: object, **kwargs: object
        ) -> None:
            assert kwargs.get("gtf_file") == str(gtf_file)
            assert kwargs.get("padding") == 1500
            Path(out_txt).write_text(
                "cDNA start position in plasmid: 1\ncDNA end position in plasmid: 40\n",
                encoding="utf-8",
            )

        prep_argv = [
            "-c",
            ns["HUMAN_FASTA"],
            ns["PLASMIDS_FILE"],
            ns["CACHE_DIR"],
            str(prep_kwargs["padding"]),
            "1" if prep_kwargs["overwrite"] else "0",
            "2",
            ns["GTF_FILE"] or "",
            str(prep_kwargs["shift_bases"]),
            "1" if prep_kwargs["generate_shifted"] else "0",
        ]
        monkeypatch.setattr(sys, "argv", prep_argv)
        with (
            patch("plasmicheck.scripts.run_pipeline.quality_control"),
            patch(
                "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
                side_effect=fake_convert,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.create_indexes",
                side_effect=fake_create_indexes,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.spliced_alignment",
                side_effect=fake_spliced_alignment,
            ) as mock_spliced_prep,
            patch(
                "plasmicheck.scripts.run_pipeline.extract_human_reference",
                side_effect=fake_extract_human,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions",
                side_effect=fake_extract_cdna,
            ),
        ):
            exec(compile(prep_one_liner, str(SMK_PATH), "exec"), {})

        mock_spliced_prep.assert_called_once()
        spliced_keys_written = [
            k for k, v in CacheManager(cache_dir).manifest.items() if v.get("type") == "spliced"
        ]
        assert len(spliced_keys_written) == 1

        seq_fq = tmp_path / "sample.fastq"
        seq_fq.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")
        out_dir = tmp_path / "output"

        scan_argv = [
            "plasmicheck",
            "pipeline",
            "-hf",
            str(human_fasta),
            "-pf",
            str(plasmids_txt),
            "-sf1",
            str(seq_fq),
            "--sample-name",
            "S1",
            "-o",
            str(out_dir),
            "--threads",
            "2",
            *shlex.split(ns["build_pipeline_flags"](cfg)),
        ]
        monkeypatch.setattr(sys, "argv", scan_argv)

        caplog.set_level(logging.INFO)
        caplog.clear()
        with (
            patch.object(cli, "setup_logging"),
            patch.object(cli, "check_requirements"),
            patch("plasmicheck.scripts.run_pipeline.quality_control"),
            patch(
                "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
                side_effect=fake_convert,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.create_indexes",
                side_effect=fake_create_indexes,
            ),
            patch("plasmicheck.scripts.run_pipeline.spliced_alignment") as mock_spliced_scan,
            patch("plasmicheck.scripts.run_pipeline.extract_human_reference") as mock_ext_scan,
            patch(
                "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"
            ) as mock_cdna_scan,
            patch("plasmicheck.scripts.run_pipeline.pysam.index"),
            patch("plasmicheck.scripts.run_pipeline.align_reads"),
            patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            patch("plasmicheck.scripts.run_pipeline.generate_report"),
            patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
        ):
            cli.main()

        mock_spliced_scan.assert_not_called()
        mock_ext_scan.assert_not_called()
        mock_cdna_scan.assert_not_called()
        assert "Cache HIT for spliced alignment" in caplog.text
        assert "Cache MISS for spliced alignment" not in caplog.text

    def test_multi_construct_bam_single_extract_single_bedcov_zero_genome_loads(
        self, tmp_path: Path
    ) -> None:
        """Verify a multi-construct BAM scan extracts the BAM once, calls bedcov once (even when median is 0.0), and performs zero whole-genome index builds or spliced alignments."""
        from unittest.mock import patch

        from plasmicheck.scripts.run_pipeline import prepare_references, run_pipeline

        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr11\n" + "ACGT" * 50 + "\n", encoding="utf-8")
        constructs: list[Path] = []
        for idx, gene in enumerate(("pKCNJ5", "pCACNA1H", "pCLCN2"), start=1):
            gb = tmp_path / f"{gene}.gb"
            gb.write_text(
                f"LOCUS       {gene} 40 bp DNA circular\nORIGIN\n 1 {'acgt' * (9 + idx)}\n//\n",
                encoding="utf-8",
            )
            constructs.append(gb)
        plasmids_txt = tmp_path / "plasmids.txt"
        plasmids_txt.write_text("\n".join(str(c) for c in constructs) + "\n", encoding="utf-8")

        bam = tmp_path / "sample.bam"
        bam.write_bytes(b"BAM" * 64)
        (tmp_path / "sample.bam.bai").write_bytes(b"BAI")
        regions_bed = tmp_path / "loci.bed"
        regions_bed.write_text("chr11\t100\t500\n", encoding="utf-8")
        target_bed = tmp_path / "targets.bed"
        target_bed.write_text("chr11\t100\t200\n", encoding="utf-8")
        cache_dir = tmp_path / "cache"

        def fake_convert(infile: str, outfile: str, *args: object, **kwargs: object) -> None:
            Path(outfile).write_text(
                f">{Path(infile).stem}\n" + "ACGT" * 10 + "\n", encoding="utf-8"
            )

        def fake_create_indexes(fasta: str, *args: object, **kwargs: object) -> None:
            base = os.path.splitext(fasta)[0]
            Path(f"{base}.mmi").write_bytes(b"MMI_INDEX")
            Path(f"{fasta}.fai").write_text("seq\t40\t6\t40\t41\n", encoding="utf-8")

        def fake_spliced(
            human_idx: str, plasmid_fa: str, out_bam: str, *args: object, **kwargs: object
        ) -> list[tuple[str, int, int]]:
            Path(out_bam).write_bytes(b"BAM_DATA")
            return [("chr11", 100, 200)]

        def fake_extract_human(
            human_fa: str, regions: object, out_fa: str, *args: object, **kwargs: object
        ) -> None:
            Path(out_fa).write_text(">chr11:101-200\n" + "ACGT" * 25 + "\n", encoding="utf-8")

        def fake_extract_cdna(
            plasmid_fa: str, spliced_bam: str, out_txt: str, *args: object, **kwargs: object
        ) -> None:
            Path(out_txt).write_text(
                "cDNA start position in plasmid: 1\ncDNA end position in plasmid: 40\n",
                encoding="utf-8",
            )

        with (
            patch("plasmicheck.scripts.run_pipeline.quality_control"),
            patch(
                "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
                side_effect=fake_convert,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.create_indexes",
                side_effect=fake_create_indexes,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.spliced_alignment",
                side_effect=fake_spliced,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.extract_human_reference",
                side_effect=fake_extract_human,
            ),
            patch(
                "plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions",
                side_effect=fake_extract_cdna,
            ),
        ):
            prepare_references(str(human_fasta), str(plasmids_txt), cache_dir=str(cache_dir))

        collate_cmds: list[str] = []

        def fake_run_cmd(cmd: str, *args: object, **kwargs: object) -> None:
            collate_cmds.append(cmd)
            if ">" in cmd:
                target = cmd.split(">")[-1].strip().strip("'\"")
                Path(target).write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")

        with (
            patch("plasmicheck.scripts.run_pipeline.quality_control"),
            patch(
                "plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta",
                side_effect=fake_convert,
            ),
            patch("plasmicheck.scripts.create_indexes.subprocess.run") as mock_idx_subproc,
            patch("plasmicheck.scripts.run_pipeline.spliced_alignment") as mock_spliced_scan,
            patch("plasmicheck.scripts.utils.run_command", side_effect=fake_run_cmd),
            patch(
                "plasmicheck.scripts.run_pipeline.compute_bed_median_depth", return_value=0.0
            ) as mock_bed_dp,
            patch(
                "plasmicheck.scripts.compare_alignments.compute_bed_median_depth"
            ) as mock_cmp_bed_dp,
            patch("plasmicheck.scripts.run_pipeline.pysam.index"),
            patch("plasmicheck.scripts.run_pipeline.align_reads"),
            patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            patch("plasmicheck.scripts.run_pipeline.generate_report"),
        ):
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmids_txt),
                output_folder=str(tmp_path / "output"),
                sequencing_files_r1=str(bam),
                cache_dir=str(cache_dir),
                regions=str(regions_bed),
                include_unmapped=False,
                target_bed=str(target_bed),
                md5_level="none",
            )

        assert len(collate_cmds) == 1
        assert "'*'" not in collate_cmds[0] and '"*"' not in collate_cmds[0]
        mock_bed_dp.assert_called_once_with(str(bam), str(target_bed))
        mock_cmp_bed_dp.assert_not_called()
        mock_spliced_scan.assert_not_called()
        mock_idx_subproc.assert_not_called()

    def test_distinct_samples_sharing_bam_basename_do_not_collide(self, tmp_path: Path) -> None:
        """Verify samples: {A: /a/sample.bam, B: /b/sample.bam} write isolated outputs and preserve sample IDs."""
        from unittest.mock import patch

        from plasmicheck.scripts.generate_summary_reports import read_compare_outputs
        from plasmicheck.scripts.run_pipeline import run_pipeline

        ns = _load_smk_namespace(config={"samples": {"A": "/a/sample.bam", "B": "/b/sample.bam"}})
        assert ns["SAMPLE_BAMS"] == {"A": "/a/sample.bam", "B": "/b/sample.bam"}

        # Colliding sanitized keys in dict must raise ValueError
        with pytest.raises(ValueError, match="Duplicate sample key"):
            ns["resolve_samples"]({"samples": {"S/1": "/a/sample.bam", "S:1": "/b/sample.bam"}})

        dir_a = tmp_path / "a"
        dir_b = tmp_path / "b"
        dir_a.mkdir()
        dir_b.mkdir()
        bam_a = dir_a / "sample.bam"
        bam_b = dir_b / "sample.bam"
        bam_a.write_text("", encoding="utf-8")
        bam_b.write_text("", encoding="utf-8")

        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        plasmid_gb = tmp_path / "p1.gb"
        plasmid_gb.write_text("LOCUS       p1 1 bp\nORIGIN\n 1 a\n//\n", encoding="utf-8")
        out_dir = tmp_path / "output"

        def write_fake_compare(
            plasmid_bam: str,
            human_bam: str,
            output_basename: str,
            *args: object,
            **kwargs: object,
        ) -> None:
            ratio = "1.5" if f"{os.sep}A{os.sep}" in output_basename else "0.1"
            verdict = (
                "Sample is contaminated with plasmid DNA"
                if ratio == "1.5"
                else "Sample is not contaminated with plasmid DNA"
            )
            Path(f"{output_basename}.summary.tsv").write_text(
                f"Category\tCount\nPlasmid\t10\nHuman\t10\nRatio\t{ratio}\nVerdict\t{verdict}\n",
                encoding="utf-8",
            )
            Path(f"{output_basename}.reads_assignment.tsv").write_text(
                "ReadName\tAssignedTo\tPlasmidScore\tHumanScore\n",
                encoding="utf-8",
            )

        with (
            patch("plasmicheck.scripts.run_pipeline.quality_control"),
            patch("plasmicheck.scripts.run_pipeline.create_indexes"),
            patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            patch("plasmicheck.scripts.run_pipeline.align_reads"),
            patch(
                "plasmicheck.scripts.run_pipeline.compare_alignments",
                side_effect=write_fake_compare,
            ),
            patch("plasmicheck.scripts.run_pipeline.generate_report") as mock_report,
        ):
            for sample_id, bam_path in (("A", bam_a), ("B", bam_b)):
                run_pipeline(
                    human_fasta=str(human_fasta),
                    plasmid_files=str(plasmid_gb),
                    output_folder=str(out_dir),
                    sequencing_files_r1=str(bam_path),
                    sample_name=sample_id,
                    no_cache=True,
                    md5_level="none",
                )

        assert (out_dir / "A" / "p1" / "data" / "comparison_result.summary.tsv").is_file()
        assert (out_dir / "B" / "p1" / "data" / "comparison_result.summary.tsv").is_file()
        assert not (out_dir / "sample").exists()

        reported_sample_names = [c.kwargs.get("sample_name") for c in mock_report.call_args_list]
        assert reported_sample_names == ["A", "B"]

        df = read_compare_outputs(str(out_dir))
        assert set(df["Sample"].unique()) == {"A", "B"}

    def test_partial_construct_failure_blocks_done_marker(self, tmp_path: Path) -> None:
        """Verify one successful plus one failed construct fails verification and prevents .done creation."""
        import subprocess
        from unittest.mock import patch

        from plasmicheck.scripts.run_pipeline import run_pipeline

        ns = _load_smk_namespace(config={"samples": {"S1": "s1.bam"}})
        expected_sample_outputs = ns["expected_sample_outputs"]
        verify_sample_outputs = ns["verify_sample_outputs"]

        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
        p_ok = tmp_path / "p_ok.gb"
        p_fail = tmp_path / "p_fail.gb"
        p_ok.write_text("LOCUS       pok 1 bp\nORIGIN\n 1 a\n//\n", encoding="utf-8")
        p_fail.write_text("LOCUS       pfail 1 bp\nORIGIN\n 1 a\n//\n", encoding="utf-8")
        plasmids_txt = tmp_path / "plasmids.txt"
        plasmids_txt.write_text(f"{p_ok}\n{p_fail}\n", encoding="utf-8")
        bam = tmp_path / "s1.bam"
        bam.write_text("", encoding="utf-8")
        out_dir = tmp_path / "output"

        def compare_one_fails(
            plasmid_bam: str,
            human_bam: str,
            output_basename: str,
            *args: object,
            **kwargs: object,
        ) -> None:
            if "p_fail" in output_basename:
                raise RuntimeError("Simulated construct failure")
            Path(f"{output_basename}.summary.tsv").write_text(
                "Category\tCount\nPlasmid\t1\nHuman\t10\nRatio\t0.1\nVerdict\tClean\n",
                encoding="utf-8",
            )
            Path(f"{output_basename}.reads_assignment.tsv").write_text(
                "ReadName\tAssignedTo\tPlasmidScore\tHumanScore\n",
                encoding="utf-8",
            )

        with (
            patch("plasmicheck.scripts.run_pipeline.quality_control"),
            patch("plasmicheck.scripts.run_pipeline.create_indexes"),
            patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            patch("plasmicheck.scripts.run_pipeline.align_reads"),
            patch(
                "plasmicheck.scripts.run_pipeline.compare_alignments",
                side_effect=compare_one_fails,
            ),
            patch("plasmicheck.scripts.run_pipeline.generate_report"),
        ):
            # run_pipeline logs partial failure and continues
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmids_txt),
                output_folder=str(out_dir),
                sequencing_files_r1=str(bam),
                sample_name="S1",
                no_cache=True,
                md5_level="none",
            )

        constructs = [str(p_ok), str(p_fail)]
        with pytest.raises(RuntimeError, match="p_fail"):
            verify_sample_outputs(str(out_dir), "S1", constructs)

        # Also test the exact shell guard from rule scan_sample
        done_file = tmp_path / "logs" / "S1.done"
        done_file.parent.mkdir(parents=True, exist_ok=True)
        expected_files = expected_sample_outputs(str(out_dir), "S1", constructs)
        quoted_files = " ".join(f"'{f}'" for f in expected_files)
        proc = subprocess.run(
            [
                "bash",
                "-c",
                f'for f in {quoted_files}; do [ -s "$f" ] || exit 1; done; touch "{done_file}"',
            ],
            capture_output=True,
            text=True,
        )
        assert proc.returncode != 0
        assert not done_file.exists()

    def test_summary_manifest_excludes_stale_removed_samples_and_constructs(
        self, tmp_path: Path
    ) -> None:
        """Verify removing a sample and a construct from config excludes their stale TSVs from cohort summaries."""
        from plasmicheck.output_paths import OutputPaths
        from plasmicheck.scripts.generate_summary_reports import read_compare_outputs

        ns = _load_smk_namespace(config={"samples": {"S1": "s1.bam", "S2": "s2.bam"}})
        write_summary_manifest = ns["write_summary_manifest"]

        out_dir = tmp_path / "output"
        for sample in ("S1", "S2"):
            for plasmid in ("p1", "p2"):
                paths = OutputPaths(output_root=out_dir, sample_name=sample, plasmid_name=plasmid)
                paths.ensure_directories()
                Path(f"{paths.comparison_prefix}.summary.tsv").write_text(
                    "Category\tCount\nPlasmid\t5\nHuman\t20\nRatio\t0.25\n"
                    "Verdict\tSample is not contaminated with plasmid DNA\n",
                    encoding="utf-8",
                )

        manifest = out_dir / ".summary_manifest.txt"

        # Initial run with 2 samples x 2 constructs
        write_summary_manifest(
            str(manifest), str(out_dir), ["S1", "S2"], ["data/p1.gb", "data/p2.gb"]
        )
        df_full = read_compare_outputs(str(manifest))
        assert set(df_full["Sample"].unique()) == {"S1", "S2"}
        assert set(df_full["Plasmid"].unique()) == {"p1", "p2"}

        # Remove S2 and p2 from config while leaving their files on disk in out_dir
        write_summary_manifest(str(manifest), str(out_dir), ["S1"], ["data/p1.gb"])
        df_pruned = read_compare_outputs(str(manifest))
        assert set(df_pruned["Sample"].unique()) == {"S1"}
        assert set(df_pruned["Plasmid"].unique()) == {"p1"}


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
        assert (
            len(
                ns["expected_cohort_summary_files"](
                    ns["OUTPUT_DIR"], ns["SAMPLES"], ns["CONSTRUCT_FILES"]
                )
            )
            == 97 * 21
        )

        done_targets = ns["expand"](f"{ns['LOG_DIR']}/{{sample}}.done", sample=ns["SAMPLES"])
        assert len(done_targets) == 97  # 1 job per sample BAM, NOT 97 * 21 = 2,037

        # Verify rule definitions in snakemake/plasmicheck.smk
        smk_text = SMK_PATH.read_text(encoding="utf-8")
        defined_rules = re.findall(r"^rule\s+(\w+)\s*:", smk_text, flags=re.MULTILINE)
        assert defined_rules == ["all", "prepare_constructs", "scan_sample", "summary_reports"]
        assert "itertools.product" not in smk_text
        assert "-pf {input.plasmids_file:q}" in smk_text and "-sf1 {input.bam:q}" in smk_text
        assert "--sample-name {wildcards.sample:q}" in smk_text

        # Verify all per-rule walltimes are < 2 hours (120 min)
        for rule_name in ("prepare_constructs", "scan_sample", "summary_reports"):
            runtime = ns["get_rule_resource"](cfg, rule_name, "runtime")
            threads = ns["get_rule_resource"](cfg, rule_name, "threads")
            mem_mb = ns["get_rule_resource"](cfg, rule_name, "mem_mb")
            assert 0 < runtime <= 120, f"{rule_name} runtime {runtime} exceeds 120 min"
            assert threads >= 1
            assert mem_mb >= 1000
