"""Unit tests for plasmicheck.scripts.run_pipeline."""

from __future__ import annotations

import warnings
from pathlib import Path

import pytest

from plasmicheck.scripts.run_pipeline import (
    PipelinePlan,
    PipelineStep,
    SequencingInput,
    build_plan,
    get_file_list,
    print_plan,
    read_file_list,
    resolve_sequencing_inputs,
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
    def test_legacy_comma_emits_deprecation(self, tmp_path: Path) -> None:
        r1 = tmp_path / "r1.fastq"
        r2 = tmp_path / "r2.fastq"
        r1.write_text("@read\nACGT\n+\nIIII\n")
        r2.write_text("@read\nTGCA\n+\nIIII\n")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = resolve_sequencing_inputs(sequencing_files=f"{r1},{r2}")
            assert len(w) == 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "deprecated" in str(w[0].message).lower()
        assert len(result) == 1
        assert result[0].file1 == str(r1)
        assert result[0].file2 == str(r2)

    @pytest.mark.unit
    def test_legacy_single_file(self, tmp_path: Path) -> None:
        bam = tmp_path / "sample.bam"
        bam.write_text("fake")
        result = resolve_sequencing_inputs(sequencing_files=str(bam))
        assert len(result) == 1
        assert result[0].file2 is None

    @pytest.mark.unit
    def test_neither_provided_raises(self) -> None:
        with pytest.raises(ValueError, match="Either -sf or -sf1"):
            resolve_sequencing_inputs()


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
            sequencing_files=str(synthetic_data_dir / "contaminated_R1.fastq"),
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
            sequencing_files=str(synthetic_data_dir / "contaminated_R1.fastq"),
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
            sequencing_files=str(synthetic_data_dir / "contaminated_R1.fastq"),
        )
        print_plan(plan)
        captured = capsys.readouterr()
        assert "PlasmiCheck Pipeline Dry-Run" in captured.out
        assert "Execution Plan" in captured.out
        assert "Summary:" in captured.out
