"""Unit tests for plasmicheck.scripts.run_pipeline."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from unittest import mock

import pytest

from plasmicheck.output_paths import OutputPaths
from plasmicheck.scripts.run_pipeline import (
    CombinationResult,
    PipelinePlan,
    PipelineStep,
    SequencingInput,
    build_plan,
    get_file_list,
    print_plan,
    read_file_list,
    resolve_sequencing_inputs,
    run_pipeline,
    write_igv_session,
)


class TestSequencingInput:
    @pytest.mark.unit
    def test_single_end(self) -> None:
        si = SequencingInput("reads.fastq")
        assert si.file1 == "reads.fastq"
        assert si.file2 is None

    @pytest.mark.unit
    def test_paired_end(self) -> None:
        si = SequencingInput("r1.fastq", "r2.fastq")
        assert si.file1 == "r1.fastq"
        assert si.file2 == "r2.fastq"

    @pytest.mark.unit
    def test_frozen(self) -> None:
        si = SequencingInput("reads.fastq")
        with pytest.raises(AttributeError):
            si.file1 = "other.fastq"  # type: ignore[misc]


class TestResolveSequencingInputs:
    @pytest.mark.unit
    def test_new_style_paired(self, tmp_path: Path) -> None:
        r1 = tmp_path / "r1.fastq"
        r2 = tmp_path / "r2.fastq"
        r1.write_text("@read\nACGT\n+\nIIII\n")
        r2.write_text("@read\nTGCA\n+\nIIII\n")
        result = resolve_sequencing_inputs(sequencing_files_r1=str(r1), sequencing_files_r2=str(r2))
        assert len(result) == 1
        assert result[0].file1 == str(r1)
        assert result[0].file2 == str(r2)

    @pytest.mark.unit
    def test_new_style_file_list(self, tmp_path: Path) -> None:
        r1a = tmp_path / "a_R1.fastq"
        r1b = tmp_path / "b_R1.fastq"
        r2a = tmp_path / "a_R2.fastq"
        r2b = tmp_path / "b_R2.fastq"
        for f in [r1a, r1b, r2a, r2b]:
            f.write_text("@read\nACGT\n+\nIIII\n")

        list_r1 = tmp_path / "r1_files.txt"
        list_r2 = tmp_path / "r2_files.txt"
        list_r1.write_text(f"{r1a}\n{r1b}\n")
        list_r2.write_text(f"{r2a}\n{r2b}\n")

        result = resolve_sequencing_inputs(
            sequencing_files_r1=str(list_r1), sequencing_files_r2=str(list_r2)
        )
        assert len(result) == 2
        assert result[0].file2 == str(r2a)

    @pytest.mark.unit
    def test_new_style_r1_only(self, tmp_path: Path) -> None:
        r1 = tmp_path / "reads.bam"
        r1.write_text("fake")
        result = resolve_sequencing_inputs(sequencing_files_r1=str(r1))
        assert len(result) == 1
        assert result[0].file2 is None

    @pytest.mark.unit
    def test_mismatched_counts_raises(self, tmp_path: Path) -> None:
        r1 = tmp_path / "r1.fastq"
        r1.write_text("@read\nACGT\n+\nIIII\n")
        list_r2 = tmp_path / "r2_files.txt"
        list_r2.write_text(f"{r1}\n{r1}\n")  # 2 files vs 1
        with pytest.raises(ValueError, match="counts must match"):
            resolve_sequencing_inputs(sequencing_files_r1=str(r1), sequencing_files_r2=str(list_r2))

    @pytest.mark.unit
    def test_single_file_r1_only(self, tmp_path: Path) -> None:
        bam = tmp_path / "sample.bam"
        bam.write_text("fake")
        result = resolve_sequencing_inputs(sequencing_files_r1=str(bam))
        assert len(result) == 1
        assert result[0].file2 is None


class TestGetFileList:
    @pytest.mark.unit
    def test_single_file(self, tmp_path: Path) -> None:
        f = tmp_path / "sample.bam"
        f.write_text("fake")
        result = get_file_list(str(f))
        assert result == [str(f)]

    @pytest.mark.unit
    def test_txt_file_list(self, tmp_path: Path) -> None:
        a = tmp_path / "a.fastq"
        b = tmp_path / "b.fastq"
        a.write_text("@read\nACGT\n+\nIIII\n")
        b.write_text("@read\nTGCA\n+\nIIII\n")
        fl = tmp_path / "files.txt"
        fl.write_text(f"{a}\n{b}\n")
        result = get_file_list(str(fl))
        assert len(result) == 2


class TestReadFileList:
    @pytest.mark.unit
    def test_skips_blank_lines(self, tmp_path: Path) -> None:
        fl = tmp_path / "files.txt"
        fl.write_text("/path/a.fastq\n\n  \n/path/b.fastq\n")
        result = read_file_list(str(fl))
        assert result == ["/path/a.fastq", "/path/b.fastq"]


class TestPipelineStep:
    @pytest.mark.unit
    def test_defaults(self) -> None:
        step = PipelineStep(name="test", description="A test step")
        assert step.inputs == []
        assert step.outputs == []
        assert step.skip is False
        assert step.skip_reason == ""

    @pytest.mark.unit
    def test_skip(self) -> None:
        step = PipelineStep(
            name="idx", description="Index", skip=True, skip_reason="Already exists"
        )
        assert step.skip is True


class TestBuildPlan:
    @pytest.mark.unit
    def test_single_combination(self, synthetic_data_dir: Path, tmp_path: Path) -> None:
        plan = build_plan(
            human_fasta=str(synthetic_data_dir / "human_ref.fasta"),
            plasmid_files=str(synthetic_data_dir / "plasmid.gb"),
            output_folder=str(tmp_path / "output"),
            overwrite=False,
            sequencing_files_r1=str(synthetic_data_dir / "contaminated_R1.fastq"),
        )
        assert len(plan.combinations) == 1
        assert plan.total_steps == 7
        assert plan.skipped_steps == 0

    @pytest.mark.unit
    def test_overwrite_warning(self, synthetic_data_dir: Path, tmp_path: Path) -> None:
        plan = build_plan(
            human_fasta=str(synthetic_data_dir / "human_ref.fasta"),
            plasmid_files=str(synthetic_data_dir / "plasmid.gb"),
            output_folder=str(tmp_path / "output"),
            overwrite=True,
            sequencing_files_r1=str(synthetic_data_dir / "contaminated_R1.fastq"),
        )
        assert any("Overwrite" in w for w in plan.warnings)

    @pytest.mark.unit
    def test_plan_properties(self) -> None:
        plan = PipelinePlan(
            human_fasta="ref.fasta",
            plasmid_files=["p.gb"],
            sequencing_inputs=[SequencingInput("r.fastq")],
            output_folder="out",
            overwrite=False,
        )
        plan.combination_steps["combo1"] = [
            PipelineStep(name="a", description="step a"),
            PipelineStep(name="b", description="step b", skip=True, skip_reason="exists"),
        ]
        assert plan.total_steps == 2
        assert plan.skipped_steps == 1


class TestPrintPlan:
    @pytest.mark.unit
    def test_print_plan_outputs_text(
        self, synthetic_data_dir: Path, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        plan = build_plan(
            human_fasta=str(synthetic_data_dir / "human_ref.fasta"),
            plasmid_files=str(synthetic_data_dir / "plasmid.gb"),
            output_folder=str(tmp_path / "output"),
            overwrite=False,
            sequencing_files_r1=str(synthetic_data_dir / "contaminated_R1.fastq"),
        )
        print_plan(plan)
        captured = capsys.readouterr()
        assert "PlasmiCheck Pipeline Dry-Run" in captured.out
        assert "Execution Plan" in captured.out
        assert "Summary:" in captured.out


class TestPipelinePlanBuiltIndexes:
    @pytest.mark.unit
    def test_pipeline_plan_has_built_indexes_field(self) -> None:
        """Verify PipelinePlan has built_indexes field that starts as empty set."""
        plan = PipelinePlan(
            human_fasta="ref.fasta",
            plasmid_files=["p.gb"],
            sequencing_inputs=[SequencingInput("r.fastq")],
            output_folder="out",
            overwrite=False,
        )
        assert isinstance(plan.built_indexes, set)
        assert len(plan.built_indexes) == 0


class TestCombinationResultDataclass:
    @pytest.mark.unit
    def test_combination_result_success(self) -> None:
        """Verify CombinationResult dataclass for successful combination."""
        result = CombinationResult(
            combo_label="plasmid.gb x sample.bam", success=True, duration=42.5, error=None
        )
        assert result.combo_label == "plasmid.gb x sample.bam"
        assert result.success is True
        assert result.duration == 42.5
        assert result.error is None

    @pytest.mark.unit
    def test_combination_result_failure(self) -> None:
        """Verify CombinationResult dataclass for failed combination."""
        error = ValueError("Test error")
        result = CombinationResult(
            combo_label="plasmid.gb x sample.bam", success=False, duration=10.2, error=error
        )
        assert result.combo_label == "plasmid.gb x sample.bam"
        assert result.success is False
        assert result.duration == 10.2
        assert result.error is error


@pytest.mark.unit
def test_pipeline_uses_requested_threshold_for_comparison_and_report(
    synthetic_data_dir: Path, tmp_path: Path
) -> None:
    with (
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments") as compare,
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report") as report,
    ):
        run_pipeline(
            human_fasta=str(synthetic_data_dir / "human_ref.fasta"),
            plasmid_files=str(synthetic_data_dir / "plasmid.gb"),
            sequencing_files_r1=str(synthetic_data_dir / "contaminated_R1.fastq"),
            output_folder=str(tmp_path / "output"),
            threshold=1.25,
            md5_level="none",
            no_cache=True,
        )

    assert compare.call_count == 1
    assert compare.call_args.kwargs["threshold"] == 1.25
    assert report.call_args.kwargs["threshold"] == 1.25


class TestIndexDeduplication:
    @pytest.mark.unit
    def test_index_dedup_human_index_built_once(self, tmp_path: Path) -> None:
        """Verify human index is built once upfront, not per combination."""
        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid1 = tmp_path / "p1.gb"
        plasmid1.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        plasmid2 = tmp_path / "p2.gb"
        plasmid2.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq_file = tmp_path / "sample.bam"
        seq_file.write_text("")  # Empty BAM placeholder

        plasmid_list = tmp_path / "plasmids.txt"
        plasmid_list.write_text(f"{plasmid1}\n{plasmid2}\n")

        output = tmp_path / "output"

        # Mock all pipeline steps
        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes") as mock_create_idx,
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
        ):
            mock_create_idx.return_value = None

            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid_list),
                output_folder=str(output),
                sequencing_files_r1=str(seq_file),
                overwrite=False,
            )

            # Human index should be created exactly once (upfront phase)
            human_index_calls = [
                call for call in mock_create_idx.call_args_list if str(human_fasta) in str(call)
            ]
            assert len(human_index_calls) == 1, (
                f"Expected 1 human index call, got {len(human_index_calls)}"
            )

    @pytest.mark.unit
    def test_index_dedup_plasmid_indexes_still_per_combination(self, tmp_path: Path) -> None:
        """Verify plasmid indexes are still created per-combination (no dedup)."""
        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid1 = tmp_path / "p1.gb"
        plasmid1.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        plasmid2 = tmp_path / "p2.gb"
        plasmid2.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq_file = tmp_path / "sample.bam"
        seq_file.write_text("")

        plasmid_list = tmp_path / "plasmids.txt"
        plasmid_list.write_text(f"{plasmid1}\n{plasmid2}\n")

        output = tmp_path / "output"

        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes") as mock_create_idx,
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
        ):
            mock_create_idx.return_value = None

            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid_list),
                output_folder=str(output),
                sequencing_files_r1=str(seq_file),
                overwrite=False,
            )

            # Plasmid indexes: 2 plasmid FASTAs + 2 spliced FASTAs = 4 total
            # (p1.fasta, p2.fasta, spliced_reference.fasta x2)
            # Plus human index (1) = 5 total
            assert mock_create_idx.call_count >= 4, (
                f"Expected at least 4 plasmid index calls, got {mock_create_idx.call_count}"
            )


class TestBatchResilience:
    @pytest.mark.unit
    def test_batch_resilience_continues_on_failure(self, tmp_path: Path) -> None:
        """Verify pipeline continues processing after one combination fails."""
        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid1 = tmp_path / "p1.gb"
        plasmid1.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        plasmid2 = tmp_path / "p2.gb"
        plasmid2.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq_file = tmp_path / "sample.bam"
        seq_file.write_text("")

        plasmid_list = tmp_path / "plasmids.txt"
        plasmid_list.write_text(f"{plasmid1}\n{plasmid2}\n")

        output = tmp_path / "output"

        # Make compare_alignments fail on first call, succeed on second
        call_count = 0

        def compare_side_effect(*args, **kwargs):  # type: ignore[no-untyped-def]
            nonlocal call_count
            call_count += 1
            if call_count == 1:
                raise ValueError("First combination failed")
            return None

        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch(
                "plasmicheck.scripts.run_pipeline.compare_alignments",
                side_effect=compare_side_effect,
            ),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report") as mock_report,
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
        ):
            # Should NOT raise despite first combination failing
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid_list),
                output_folder=str(output),
                sequencing_files_r1=str(seq_file),
                overwrite=False,
            )

            # Second combination should have succeeded and generated report
            assert mock_report.call_count == 1, "Second combination should have generated report"

    @pytest.mark.unit
    def test_batch_resilience_all_fail_raises(self, tmp_path: Path) -> None:
        """Verify pipeline raises RuntimeError if ALL combinations fail."""
        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid1 = tmp_path / "p1.gb"
        plasmid1.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        plasmid2 = tmp_path / "p2.gb"
        plasmid2.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq_file = tmp_path / "sample.bam"
        seq_file.write_text("")

        plasmid_list = tmp_path / "plasmids.txt"
        plasmid_list.write_text(f"{plasmid1}\n{plasmid2}\n")

        output = tmp_path / "output"

        # Make compare_alignments fail on ALL calls
        def compare_side_effect(*args, **kwargs):  # type: ignore[no-untyped-def]
            raise ValueError("All combinations failed")

        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch(
                "plasmicheck.scripts.run_pipeline.compare_alignments",
                side_effect=compare_side_effect,
            ),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
            pytest.raises(RuntimeError, match="All 2 combinations failed"),
        ):
            # Should raise RuntimeError when all combinations fail
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid_list),
                output_folder=str(output),
                sequencing_files_r1=str(seq_file),
                overwrite=False,
            )

    @pytest.mark.unit
    def test_per_combination_timing_logged(
        self, tmp_path: Path, caplog: pytest.LogCaptureFixture
    ) -> None:
        """Verify per-combination timing is logged at INFO level."""
        import logging

        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid = tmp_path / "p.gb"
        plasmid.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq_file = tmp_path / "sample.bam"
        seq_file.write_text("")

        output = tmp_path / "output"

        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
            caplog.at_level(logging.INFO),
        ):
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid),
                output_folder=str(output),
                sequencing_files_r1=str(seq_file),
                overwrite=False,
            )

            # Check for timing log message (format: "p.gb x sample.bam: NN.Ns")
            timing_logs = [
                record.message
                for record in caplog.records
                if "p.gb x sample.bam:" in record.message and "s" in record.message
            ]
            assert len(timing_logs) >= 1, "Expected at least one timing log message"


# --- Cache integration tests ---


class TestPipelineCacheIntegration:
    @pytest.mark.unit
    def test_pipeline_cache_args_passed(self) -> None:
        """Verify that run_pipeline accepts cache_dir and no_cache parameters."""
        # This is a signature test - just verify the function accepts the parameters
        import inspect

        sig = inspect.signature(run_pipeline)
        assert "cache_dir" in sig.parameters
        assert "no_cache" in sig.parameters

    @pytest.mark.unit
    def test_no_cache_flag_skips_initialization(self, tmp_path: Path) -> None:
        """Verify that no_cache=True prevents CacheManager instantiation."""
        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid = tmp_path / "p.gb"
        plasmid.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq_file = tmp_path / "sample.bam"
        seq_file.write_text("")

        output = tmp_path / "output"

        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
            mock.patch("plasmicheck.scripts.run_pipeline.CacheManager") as mock_cache_mgr,
        ):
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid),
                output_folder=str(output),
                sequencing_files_r1=str(seq_file),
                overwrite=False,
                no_cache=True,
            )

            # CacheManager should NOT have been instantiated
            mock_cache_mgr.assert_not_called()

    @pytest.mark.unit
    def test_cache_dir_override(self, tmp_path: Path) -> None:
        """Verify that cache_dir parameter overrides default cache location."""
        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid = tmp_path / "p.gb"
        plasmid.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq_file = tmp_path / "sample.bam"
        seq_file.write_text("")

        output = tmp_path / "output"
        custom_cache = tmp_path / "custom_cache"

        # Mock CacheManager to return None for cache misses
        mock_cache = mock.MagicMock()
        mock_cache.get_index.return_value = None
        mock_cache.get_spliced.return_value = None

        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
            mock.patch(
                "plasmicheck.scripts.run_pipeline.CacheManager", return_value=mock_cache
            ) as mock_cache_mgr,
        ):
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid),
                output_folder=str(output),
                sequencing_files_r1=str(seq_file),
                overwrite=False,
                cache_dir=str(custom_cache),
            )

            # CacheManager should have been called with custom path
            mock_cache_mgr.assert_called_once()
            call_args = mock_cache_mgr.call_args
            assert str(custom_cache) in str(call_args)

    @pytest.mark.unit
    def test_cache_initialization_failure_continues(self, tmp_path: Path) -> None:
        """Verify pipeline continues gracefully when cache initialization fails."""
        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid = tmp_path / "p.gb"
        plasmid.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq_file = tmp_path / "sample.bam"
        seq_file.write_text("")

        output = tmp_path / "output"

        # Make CacheManager.__init__ raise OSError
        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
            mock.patch(
                "plasmicheck.scripts.run_pipeline.CacheManager",
                side_effect=OSError("Permission denied"),
            ),
        ):
            # Should NOT raise - pipeline continues without cache
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid),
                output_folder=str(output),
                sequencing_files_r1=str(seq_file),
                overwrite=False,
            )

    @pytest.mark.unit
    def test_cache_within_run_dedup(self, tmp_path: Path) -> None:
        """Verify within-run deduplication via cache mechanism.

        When the same plasmid appears in multiple combinations:
        - First combo: cache MISS, build index, store in cache
        - Second combo: cache HIT, copy from cache, skip build
        """
        # Create minimal dummy files
        human_fasta = tmp_path / "human.fasta"
        human_fasta.write_text(">chr1\nACGT\n")

        plasmid = tmp_path / "p.gb"
        plasmid.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

        seq1 = tmp_path / "s1.bam"
        seq2 = tmp_path / "s2.bam"
        seq1.write_text("")
        seq2.write_text("")

        seq_list = tmp_path / "seqs.txt"
        seq_list.write_text(f"{seq1}\n{seq2}\n")

        output = tmp_path / "output"

        # Mock CacheManager to track get_index/put_index calls
        mock_cache = mock.MagicMock()
        mock_cache.get_spliced.return_value = None
        # First call to get_index returns None (MISS), second returns cached data (HIT)
        cached_mmi = tmp_path / "cached.mmi"
        cached_fai = tmp_path / "cached.fai"
        cached_mmi.write_text("fake mmi")
        cached_fai.write_text("fake fai")

        mock_cache.get_index.side_effect = [
            None,  # First combination: MISS
            {"mmi": cached_mmi, "fai": cached_fai},  # Second combination: HIT
        ]

        with (
            mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
            mock.patch("plasmicheck.scripts.run_pipeline.create_indexes") as mock_create_idx,
            mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
            mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
            mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
            mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
            mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
            mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
            mock.patch("plasmicheck.scripts.run_pipeline.write_md5sum"),
            mock.patch("plasmicheck.scripts.run_pipeline.shutil.copy2") as mock_copy,
            mock.patch("plasmicheck.scripts.run_pipeline.CacheManager", return_value=mock_cache),
        ):
            run_pipeline(
                human_fasta=str(human_fasta),
                plasmid_files=str(plasmid),
                output_folder=str(output),
                sequencing_files_r1=str(seq_list),
                overwrite=False,
            )

            # Verify cache.get_index was called twice (once per combination)
            assert mock_cache.get_index.call_count == 2

            # Verify cache.put_index was called once (after first build)
            assert mock_cache.put_index.call_count == 1

            # Verify shutil.copy2 was called on second combination (cache HIT)
            # Should be 2 calls: one for .mmi, one for .fai
            copy_calls = [call for call in mock_copy.call_args_list if "cached" in str(call)]
            assert len(copy_calls) >= 2, "Expected at least 2 copy calls for cache HIT"

            # Count create_indexes calls for plasmid FASTA
            # First combo builds, second combo skips due to cache hit
            # Plus spliced FASTA indexes (2x) and human index (1x) = variable count
            # Key verification: plasmid FASTA indexed only once despite 2 combos
            plasmid_fasta_idx_calls = [
                call for call in mock_create_idx.call_args_list if "p.fasta" in str(call)
            ]
            # Should be 1 call for plasmid FASTA (first combo only)
            assert len(plasmid_fasta_idx_calls) == 1, (
                f"Expected 1 plasmid index build, got {len(plasmid_fasta_idx_calls)}"
            )


@pytest.mark.unit
def test_pipeline_passes_cdna_positions_file_to_compare_alignments(
    synthetic_data_dir: Path, tmp_path: Path
) -> None:
    """Verify run_pipeline passes cdna_positions_file to compare_alignments."""
    custom_cdna = str(tmp_path / "custom_cdna.tsv")
    with (
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments") as compare,
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
    ):
        run_pipeline(
            human_fasta=str(synthetic_data_dir / "human_ref.fasta"),
            plasmid_files=str(synthetic_data_dir / "plasmid.gb"),
            sequencing_files_r1=str(synthetic_data_dir / "contaminated_R1.fastq"),
            output_folder=str(tmp_path / "output"),
            cdna_output=custom_cdna,
            md5_level="none",
            no_cache=True,
        )

    assert compare.call_count == 1
    assert compare.call_args.kwargs["cdna_positions_file"] == custom_cdna


@pytest.mark.unit
def test_write_igv_session_relative_paths(tmp_path: Path) -> None:
    """Verify write_igv_session produces valid XML with relative paths."""
    paths = OutputPaths(tmp_path, "sample1", "plasmid1")
    write_igv_session(paths)

    assert paths.igv_session.exists()
    tree = ET.parse(paths.igv_session)
    root = tree.getroot()
    assert root.tag == "Session"
    assert root.attrib["genome"] == paths.plasmid_fasta.relative_to(paths.base).as_posix()
    assert root.attrib["version"] == "8"

    resources = root.find("Resources")
    assert resources is not None
    resource = resources.find("Resource")
    assert resource is not None
    assert resource.attrib["path"] == paths.plasmid_bam.relative_to(paths.base).as_posix()

    content = paths.igv_session.read_text()
    assert str(tmp_path) not in content
    assert str(paths.base) not in content


@pytest.mark.unit
@pytest.mark.parametrize("keep_intermediate,expected_exists", [(True, True), (False, False)])
def test_pipeline_generates_igv_session_when_keep_intermediate(
    synthetic_data_dir: Path, tmp_path: Path, keep_intermediate: bool, expected_exists: bool
) -> None:
    """Verify IGV session generation is gated on keep_intermediate."""
    output_dir = tmp_path / f"output_{keep_intermediate}"

    def fake_align(index: str, reads: str, out_bam: str, *args: object, **kwargs: object) -> None:
        Path(out_bam).parent.mkdir(parents=True, exist_ok=True)
        Path(out_bam).write_bytes(b"BAM data")

    with (
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads", side_effect=fake_align),
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
    ):
        run_pipeline(
            human_fasta=str(synthetic_data_dir / "human_ref.fasta"),
            plasmid_files=str(synthetic_data_dir / "plasmid.gb"),
            sequencing_files_r1=str(synthetic_data_dir / "contaminated_R1.fastq"),
            output_folder=str(output_dir),
            keep_intermediate=keep_intermediate,
            md5_level="none",
            no_cache=True,
        )

    igv_files = list(output_dir.glob("**/igv_session.xml"))
    if expected_exists:
        assert len(igv_files) == 1
    else:
        assert len(igv_files) == 0


@pytest.mark.unit
def test_pipeline_anonymize_flag(
    synthetic_data_dir: Path, tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """Verify run_pipeline --anonymize forces keep_intermediate=False, logs warning, and redacts paths."""
    output_dir = tmp_path / "output_anon"

    with (
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments") as mock_compare,
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report") as mock_report,
    ):
        run_pipeline(
            human_fasta=str(synthetic_data_dir / "human_ref.fasta"),
            plasmid_files=str(synthetic_data_dir / "plasmid.gb"),
            sequencing_files_r1=str(synthetic_data_dir / "contaminated_R1.fastq"),
            output_folder=str(output_dir),
            keep_intermediate=True,
            anonymize=True,
            md5_level="none",
            no_cache=True,
        )

    # Check warning was logged
    assert any("--anonymize deletes intermediate BAMs" in r.message for r in caplog.records)

    # Check compare_alignments received anonymize=True and 16-byte key
    assert mock_compare.call_count == 1
    assert mock_compare.call_args.kwargs["anonymize"] is True
    key = mock_compare.call_args.kwargs["anonymize_key"]
    assert isinstance(key, bytes)
    assert len(key) == 16

    # Check generate_report received anonymize=True and redacted paths
    assert mock_report.call_count == 1
    kwargs = mock_report.call_args.kwargs
    assert kwargs["anonymize"] is True
    assert kwargs["human_fasta"] == "human_ref.fasta"
    assert kwargs["plasmid_gb"] == "plasmid.gb"
    assert kwargs["sequencing_file"] == "contaminated_R1.fastq"
    assert "/" not in kwargs["human_fasta"]
    assert "/" not in kwargs["plasmid_gb"]
    assert "/" not in kwargs["sequencing_file"]


@pytest.mark.unit
def test_spliced_fasta_index_caching(tmp_path: Path) -> None:
    """Verify spliced FASTA indexes (.mmi and .fai) are queried and stored in CacheManager."""
    human_fasta = tmp_path / "human.fasta"
    human_fasta.write_text(">chr1\nACGT\n")

    plasmid = tmp_path / "p.gb"
    plasmid.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

    seq = tmp_path / "s.fastq"
    seq.write_text("@read\nACGT\n+\nIIII\n")

    output = tmp_path / "output"

    # Part 1: Cache MISS on spliced FASTA index
    mock_cache = mock.MagicMock()
    mock_cache.get_index.return_value = None
    mock_cache.get_spliced.return_value = None

    with (
        mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
        mock.patch("plasmicheck.scripts.run_pipeline.CacheManager", return_value=mock_cache),
    ):
        run_pipeline(
            human_fasta=str(human_fasta),
            plasmid_files=str(plasmid),
            sequencing_files_r1=str(seq),
            output_folder=str(output),
            overwrite=False,
            md5_level="none",
        )

        get_index_paths = [call.args[0] for call in mock_cache.get_index.call_args_list]
        assert any("spliced_reference.fasta" in p for p in get_index_paths)

        put_index_paths = [call.args[0] for call in mock_cache.put_index.call_args_list]
        assert any("spliced_reference.fasta" in p for p in put_index_paths)

    # Part 2: Cache HIT on spliced FASTA index
    cached_mmi = tmp_path / "cached_spliced.mmi"
    cached_fai = tmp_path / "cached_spliced.fai"
    cached_mmi.write_text("mmi content")
    cached_fai.write_text("fai content")

    def cache_get_index(fasta_path: str) -> dict[str, Path] | None:
        if "spliced_reference.fasta" in fasta_path:
            return {"mmi": cached_mmi, "fai": cached_fai}
        return None

    mock_cache_hit = mock.MagicMock()
    mock_cache_hit.get_index.side_effect = cache_get_index
    mock_cache_hit.get_spliced.return_value = None

    with (
        mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes") as mock_create_idx_hit,
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
        mock.patch("plasmicheck.scripts.run_pipeline.CacheManager", return_value=mock_cache_hit),
    ):
        run_pipeline(
            human_fasta=str(human_fasta),
            plasmid_files=str(plasmid),
            sequencing_files_r1=str(seq),
            output_folder=str(tmp_path / "output_hit"),
            overwrite=False,
            md5_level="none",
        )

        spliced_idx_calls = [
            call
            for call in mock_create_idx_hit.call_args_list
            if "spliced_reference.fasta" in str(call)
        ]
        assert len(spliced_idx_calls) == 0


@pytest.mark.unit
def test_sample_scoped_temp_fastq_cleanup(tmp_path: Path) -> None:
    """Verify sample_temp_fastq is scoped per-sample and cleaned up even if error occurs."""
    human_fasta = tmp_path / "human.fasta"
    human_fasta.write_text(">chr1\nACGT\n")

    plasmid1 = tmp_path / "p1.gb"
    plasmid1.write_text("LOCUS       test1 1 bp\nORIGIN\n 1 a\n//\n")
    plasmid2 = tmp_path / "p2.gb"
    plasmid2.write_text("LOCUS       test2 1 bp\nORIGIN\n 1 a\n//\n")

    plasmid_list = tmp_path / "plasmids.txt"
    plasmid_list.write_text(f"{plasmid1}\n{plasmid2}\n")

    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"A" * 128)

    output = tmp_path / "output"

    created_temp_fastqs: list[Path] = []

    def fake_run_cmd(cmd: str, *args: object, **kwargs: object) -> None:
        if ">" in cmd:
            target = cmd.split(">")[-1].strip().strip("'\"")
            fq_path = Path(target)
            fq_path.write_text("@read\nACGT\n+\nIIII\n")
            created_temp_fastqs.append(fq_path)

    with (
        mock.patch("plasmicheck.scripts.utils.run_command", side_effect=fake_run_cmd),
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
        mock.patch(
            "plasmicheck.scripts.run_pipeline.compare_alignments",
            side_effect=RuntimeError("simulated comparison error"),
        ),
        pytest.raises(RuntimeError, match="All 2 combinations failed"),
    ):
        run_pipeline(
            human_fasta=str(human_fasta),
            plasmid_files=str(plasmid_list),
            sequencing_files_r1=str(bam),
            output_folder=str(output),
            no_cache=True,
            md5_level="none",
        )

    assert len(created_temp_fastqs) == 1
    assert not str(created_temp_fastqs[0]).startswith(str(output))
    assert "plasmicheck_prep_" in str(created_temp_fastqs[0])
    assert not created_temp_fastqs[0].exists()


@pytest.mark.unit
def test_cross_sample_plasmid_prep_dedup(tmp_path: Path) -> None:
    """Verify plasmid preparation (convert, index, spliced align, cDNA) is done only once across multiple samples."""
    human_fasta = tmp_path / "human.fasta"
    human_fasta.write_text(">chr1\nACGT\n")

    plasmid = tmp_path / "shared_p.gb"
    plasmid.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

    seq1 = tmp_path / "sample1.fastq"
    seq2 = tmp_path / "sample2.fastq"
    seq1.write_text("@read1\nACGT\n+\nIIII\n")
    seq2.write_text("@read2\nACGT\n+\nIIII\n")

    seq_list = tmp_path / "seqs.txt"
    seq_list.write_text(f"{seq1}\n{seq2}\n")

    output = tmp_path / "output_dedup"

    with (
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta") as mock_conv,
        mock.patch(
            "plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]
        ) as mock_splice,
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference") as mock_extract_ref,
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions") as mock_cdna,
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads") as mock_align,
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments") as mock_compare,
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report") as mock_report,
    ):
        run_pipeline(
            human_fasta=str(human_fasta),
            plasmid_files=str(plasmid),
            sequencing_files_r1=str(seq_list),
            output_folder=str(output),
            no_cache=True,
            md5_level="none",
        )

    assert mock_conv.call_count == 1
    assert mock_splice.call_count == 1
    assert mock_extract_ref.call_count == 1
    assert mock_cdna.call_count == 1

    assert mock_align.call_count == 4
    assert mock_compare.call_count == 2
    assert mock_report.call_count == 2


@pytest.mark.unit
def test_auto_derived_and_explicit_regions_for_indexed_bam(tmp_path: Path) -> None:
    """Verify auto-derived spliced_reference.fasta regions and explicit --regions override on indexed BAM."""
    human_fasta = tmp_path / "human.fasta"
    human_fasta.write_text(">chr11\nACGT\n")

    plasmid = tmp_path / "p1.gb"
    plasmid.write_text("LOCUS       test1 1 bp\nORIGIN\n 1 a\n//\n")

    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"A" * 128)
    (tmp_path / "sample.bam.bai").write_bytes(b"BAI")

    output = tmp_path / "output_regions"
    captured_cmds: list[str] = []

    def fake_extract_ref(_hf: str, _regions: object, out_fasta: str) -> None:
        Path(out_fasta).write_text(">chr11:128890001-128920000\nACGT\n")

    def fake_run_cmd(cmd: str, *args: object, **kwargs: object) -> None:
        captured_cmds.append(cmd)
        if ">" in cmd:
            target = cmd.split(">")[-1].strip().strip("'\"")
            Path(target).write_text("@read\nACGT\n+\nIIII\n")

    with (
        mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
        mock.patch("plasmicheck.scripts.utils.run_command", side_effect=fake_run_cmd),
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch(
            "plasmicheck.scripts.run_pipeline.spliced_alignment",
            return_value=[("chr11", 128890000, 128920000)],
        ),
        mock.patch(
            "plasmicheck.scripts.run_pipeline.extract_human_reference",
            side_effect=fake_extract_ref,
        ),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads") as mock_align,
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
    ):
        run_pipeline(
            human_fasta=str(human_fasta),
            plasmid_files=str(plasmid),
            sequencing_files_r1=str(bam),
            output_folder=str(output),
            include_unmapped=True,
            no_cache=True,
            md5_level="none",
        )

    assert len(captured_cmds) == 1
    assert "samtools view" in captured_cmds[0]
    assert "chr11:128890001-128920000" in captured_cmds[0]
    assert "'*'" in captured_cmds[0]
    # Because sample_temp_fastq was created, align_reads gets regions=None
    for call in mock_align.call_args_list:
        assert call.kwargs["regions"] is None

    # Explicit --regions BED override with include_unmapped=False (default)
    captured_cmds.clear()
    bed = tmp_path / "override.bed"
    bed.write_text("chr3\t53500000\t53900000\n")

    with (
        mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
        mock.patch("plasmicheck.scripts.utils.run_command", side_effect=fake_run_cmd),
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch(
            "plasmicheck.scripts.run_pipeline.spliced_alignment",
            return_value=[("chr11", 128890000, 128920000)],
        ),
        mock.patch(
            "plasmicheck.scripts.run_pipeline.extract_human_reference",
            side_effect=fake_extract_ref,
        ),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions"),
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments"),
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
    ):
        run_pipeline(
            human_fasta=str(human_fasta),
            plasmid_files=str(plasmid),
            sequencing_files_r1=str(bam),
            output_folder=str(output / "override"),
            regions=str(bed),
            no_cache=True,
            md5_level="none",
        )

    assert len(captured_cmds) == 1
    assert "chr3:53500001-53900000" in captured_cmds[0]
    assert "chr11:128890001-128920000" not in captured_cmds[0]
    assert "'*'" not in captured_cmds[0] and '"*"' not in captured_cmds[0]


@pytest.mark.unit
def test_pipeline_forwards_gtf_vcf_and_depth_options(tmp_path: Path) -> None:
    """Verify run_pipeline forwards gtf_file/gtf_md5, vcf_file, target_bed, and depth options."""
    human_fasta = tmp_path / "human.fasta"
    human_fasta.write_text(">chr1\nACGT\n")

    plasmid = tmp_path / "p.gb"
    plasmid.write_text("LOCUS       test 1 bp\nORIGIN\n 1 a\n//\n")

    bam = tmp_path / "sample.bam"
    bam.write_bytes(b"fake bam")

    gtf = tmp_path / "mane.gtf"
    gtf.write_text('chr1\tMANE\texon\t10\t50\t.\t+\t.\ttranscript_id "tx1"; tag "MANE_Select";\n')

    target_bed = tmp_path / "exome.bed"
    target_bed.write_text("chr1\t0\t100\n")

    vcf = tmp_path / "calls.tsv"
    vcf.write_text("sample\tchr\tpos\tref\talt\n")

    mock_cache = mock.MagicMock()
    mock_cache.get_index.return_value = None
    mock_cache.get_spliced.return_value = None

    with (
        mock.patch("plasmicheck.scripts.run_pipeline.quality_control"),
        mock.patch("plasmicheck.scripts.run_pipeline.create_indexes"),
        mock.patch("plasmicheck.scripts.run_pipeline.convert_plasmidfile_to_fasta"),
        mock.patch("plasmicheck.scripts.run_pipeline.spliced_alignment", return_value=[]),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_human_reference"),
        mock.patch("plasmicheck.scripts.run_pipeline.extract_plasmid_cdna_positions") as mock_cdna,
        mock.patch("plasmicheck.scripts.run_pipeline.align_reads"),
        mock.patch(
            "plasmicheck.scripts.run_pipeline.compute_bed_median_depth", return_value=112.5
        ) as mock_bed_dp,
        mock.patch("plasmicheck.scripts.run_pipeline.compare_alignments") as mock_compare,
        mock.patch("plasmicheck.scripts.run_pipeline.generate_report"),
        mock.patch("plasmicheck.scripts.run_pipeline.CacheManager", return_value=mock_cache),
    ):
        run_pipeline(
            human_fasta=str(human_fasta),
            plasmid_files=str(plasmid),
            sequencing_files_r1=str(bam),
            output_folder=str(tmp_path / "out"),
            gtf_file=str(gtf),
            vcf_file=str(vcf),
            target_bed=str(target_bed),
            dp_locus=240.0,
            r_cohort=0.58,
            md5_level="none",
        )

    assert mock_cdna.call_args.kwargs["gtf_file"] == str(gtf)
    assert mock_cache.get_spliced.call_args.kwargs["gtf_md5"] != ""
    assert mock_cache.put_spliced.call_args.kwargs["gtf_md5"] != ""
    mock_bed_dp.assert_called_once_with(str(bam), str(target_bed))

    cmp_kwargs = mock_compare.call_args.kwargs
    assert cmp_kwargs["vcf_file"] == str(vcf)
    assert cmp_kwargs["target_bed"] == str(target_bed)
    assert cmp_kwargs["dp_sample_median"] == 112.5
    assert cmp_kwargs["dp_locus"] == 240.0
    assert cmp_kwargs["r_cohort"] == 0.58
    assert cmp_kwargs["sample_bam"] == str(bam)
