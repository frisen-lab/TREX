import logging
import subprocess
from pathlib import Path
from shutil import copytree

from trex.cli import add_file_logging
from trex.cli.qc import make_qc_report
from trex.cli.run10x import run_trex, read_excluded_clone_ids
from trex.cli.smartseq2 import run_smartseq2
from trex.cli.smartseq3 import run_smartseq3
from trex.quality_control import read_cells


def diff(expected, actual, ignore=None, recursive=False):
    args = ["diff", "-u"]
    if recursive:
        args.append("-r")
    if ignore is not None:
        args += [f"-x{name}" for name in ignore]
    subprocess.run(args + [expected, actual]).check_returncode()


def bam_diff(expected_bam, actual_bam, tmp_path):
    expected_sam = tmp_path / "expected.sam"
    actual_sam = tmp_path / "actual.sam"
    with open(expected_sam, "w") as expected, open(actual_sam, "w") as actual:
        subprocess.run(
            ["samtools", "view", "--no-PG", "-h", expected_bam], stdout=expected
        )
        subprocess.run(["samtools", "view", "--no-PG", "-h", actual_bam], stdout=actual)
    diff(expected_sam, actual_sam)


def test_run_trex(tmp_path):
    add_file_logging(tmp_path / "log.txt")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    run_trex(
        tmp_path,
        transcriptome_inputs=["tests/data/"],
        amplicon_inputs=[],
        start=694,
        end=724,
        keep_doublets=True,
        should_write_loom=True,
        should_write_umi_matrix=True,
    )
    diff(
        "tests/expected", tmp_path, ignore=["data.loom", "entries.bam"], recursive=True
    )
    bam_diff("tests/expected/entries.bam", tmp_path / "entries.bam", tmp_path)


def test_run_trex_per_cell(tmp_path):
    # Ensure --per-cell works
    add_file_logging(tmp_path / "log.txt")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    run_trex(
        tmp_path,
        transcriptome_inputs=["tests/data/"],
        amplicon_inputs=[],
        start=694,
        end=724,
        correct_per_cell=True,
    )
    diff("tests/expected_per_cell/log.txt", tmp_path / "log.txt")


def test_run_smartseq2(tmp_path):
    add_file_logging(tmp_path / "log.txt")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    run_smartseq2(
        tmp_path,
        transcriptome_inputs=[Path("tests/data/smartseq2_test.bam")],
        amplicon_inputs=[],
        start=2329,
        end=2359,
        should_write_read_matrix=True,
        readcount_threshold=1,
        chromosome="EGFP-30N",
    )
    diff("tests/expected_smartseq2/", tmp_path, ignore=["entries.bam"], recursive=True)


def test_run_smartseq3(tmp_path):
    add_file_logging(tmp_path / "log.txt")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    run_smartseq3(
        tmp_path,
        transcriptome_inputs=[Path("tests/data/smartseq3_test.bam")],
        amplicon_inputs=[],
        start=2329,
        end=2359,
        chromosome="EGFP-30N",
        should_write_umi_matrix=True,
    )
    diff("tests/expected_smartseq3/", tmp_path, ignore=["entries.bam"], recursive=True)
    assert tmp_path.joinpath("umi_count_matrix.csv").exists()


def test_qc(tmp_path):
    add_file_logging(tmp_path / "log.txt")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    pdf_path = tmp_path / "quality_report.pdf"

    copytree("tests/expected", tmp_path / "trexrun")
    make_qc_report(
        tmp_path / "trexrun",
        pdf_path,
        plot_jaccard=True,
        plot_hamming=True,
    )
    assert pdf_path.exists()


EXCLUDED_CLONE_IDS = {
    "GGTCTCCCTATACCAACAGTATCGTCTCAA",
    "GGGTTCTGGGATATTACGTTGACTTGAGAG",
}


def test_filter_clone_ids(tmp_path):
    add_file_logging(tmp_path / "log.txt")
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    run_trex(
        tmp_path,
        transcriptome_inputs=["tests/data/"],
        amplicon_inputs=[],
        start=694,
        end=724,
        correct_per_cell=True,
        excluded_clone_ids=EXCLUDED_CLONE_IDS,
    )
    cells = read_cells(tmp_path, filtered=False)

    diff("tests/expected_filter_clone_id/log.txt", tmp_path / "log.txt")
    for clone_id in EXCLUDED_CLONE_IDS:
        assert clone_id not in cells["clone_id"].values


def test_read_excluded_clone_ids():
    cloneids = read_excluded_clone_ids(Path("tests/data/exclude_list.csv"))
    assert cloneids == EXCLUDED_CLONE_IDS
