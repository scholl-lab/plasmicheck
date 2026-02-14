"""Unit tests for plasmicheck.scripts.run_pipeline."""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

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
