import os
import pickle
from collections import Counter
from pathlib import Path

import pysam


def write_alignment_stats(bam_file: Path, out_file: Path):
    """
    Read through a BAM file and collect some distributions:
      - alignment scores for uniquely-mapped reads
      - number of mismatches per uniquely-mapped read
      - ratio of matches / total for uniquely-mapped reads
      - the above three stats but for multi-mapped reads

    And write them into a pickle

    :param bam_file: aligned BAM file to check (output from STAR)
    :param out_file: file to output distributions (as a pickle)
    """
    aligned_bam = pysam.AlignmentFile(bam_file, mode="rb")

    mp = Counter()
    for a in aligned_bam:
        mp[a.qname] += 1

    aligned_bam.reset()

    unique_score = Counter()
    unique_mismatch = Counter()
    unique_ratio = Counter()
    multi_score = Counter()
    multi_mismatch = Counter()
    multi_ratio = Counter()

    for a in aligned_bam:
        if mp[a.qname] == 1:
            if a.has_tag("AS"):
                unique_score[a.get_tag("AS")] += 1
            if a.has_tag("nM"):
                unique_mismatch[a.get_tag("nM")] += 1
            r = 100 * round(a.get_cigar_stats()[0][0] / a.qlen)
            unique_ratio[r] += 1
        else:
            if a.has_tag("AS"):
                multi_score[a.get_tag("AS")] += 1
            if a.has_tag("nM"):
                multi_mismatch[a.get_tag("nM")] += 1
            r = 100 * round(a.get_cigar_stats()[0][0] / a.qlen)
            multi_ratio[r] += 1

    with out_file.open("wb") as out:
        pickle.dump(
            (
                unique_score,
                unique_mismatch,
                unique_ratio,
                multi_score,
                multi_mismatch,
                multi_ratio,
            ),
            out,
        )


def parse_star_log(log_file: Path) -> dict[str, str]:
    """
    Read the Log.final.out file from STAR and return the parsed data

    :param log_file: Path to STAR's final output log
    :return: Dictionary containing the run stats
    """
    with log_file.open() as fh:
        rows = [line.strip().split("|\t") for line in fh if "|" in line]
        log_data = {r[0].strip(): r[1].strip() for r in rows}

    return log_data


def combine_alignment_stats(
    stats_files: list[Path],
    star_logs: list[Path],
    out_base: Path,
):
    """
    Combines the alignment statistics files from multiple lanes into one aggregate,
    and outputs aggregated versions into text files. Also aggregates the STAR log
    reports into one combined file.

    :param stats_files: List of pickles output by check_alignment_quality
    :param star_logs: List of STAR log output files
    :param out_base: Base path for the summary files to write
    """
    unique_score = Counter()
    unique_mismatch = Counter()
    unique_ratio = Counter()
    multi_score = Counter()
    multi_mismatch = Counter()
    multi_ratio = Counter()

    for stat_file in stats_files:
        with stat_file.open("rb") as fh:
            us, um, ur, ms, mm, mr = pickle.load(fh)
            unique_score += us
            unique_mismatch += um
            unique_ratio += ur
            multi_score += ms
            multi_mismatch += mm
            multi_ratio += mr

        os.remove(stat_file)

    # write out the combined statistics
    for suffix, dist in (
        (".unique.score", unique_score),
        (".unique.mismatch", unique_mismatch),
        (".unique.ratio", unique_ratio),
        (".multi.score", unique_score),
        (".multi.mismatch", multi_mismatch),
        (".multi.ratio", multi_ratio),
    ):
        with out_base.with_suffix(suffix).open("w") as out:
            for k, v in sorted(dist.items()):
                print(f"{k}\t{v}", file=out)

    total_reads = 0
    unique_reads = 0
    multi_reads = 0
    too_many_reads = 0

    # collect some stats from star log files
    for star_log in star_logs:
        log_data = parse_star_log(star_log)
        total_reads += int(log_data["Number of input reads"])
        unique_reads += int(log_data["Uniquely mapped reads number"])
        multi_reads += int(log_data["Number of reads mapped to multiple loci"])
        too_many_reads += int(log_data["Number of reads mapped to too many loci"])

    mismatch1 = unique_mismatch[1] + multi_mismatch[1]
    mismatch2 = unique_mismatch[2] + multi_mismatch[2]
    mismatch3 = unique_mismatch[3] + multi_mismatch[3]

    with out_base.with_suffix(".mapping_rate.txt").open("w") as out:
        print(f"total_reads\t{total_reads}", file=out)
        print(f"unique_aligned_reads\t{unique_reads}", file=out)
        print(f"unique_aligned_ratio\t{100 * unique_reads / total_reads:.3g}", file=out)
        print(f"multi_aligned_reads\t{multi_reads}", file=out)
        print("multi_aligned_ratio\t{100 * multi_reads / total_reads:.3g}", file=out)
        print(f"too_many_aligned_reads\t{too_many_reads}", file=out)
        print(
            f"too_many_aligned_ratio\t{100 * too_many_reads / total_reads:.3g}",
            file=out,
        )
        print(f"mismatch1_rate\t{100 * mismatch1 / total_reads:.3g}", file=out)
        print(f"mismatch2_rate\t{100 * mismatch2 / total_reads:.3g}", file=out)
        print(f"mismatch3_rate\t{100 * mismatch3 / total_reads:.3g}", file=out)
