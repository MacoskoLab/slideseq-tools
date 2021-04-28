import pathlib
import pickle
from collections import Counter

import pysam


def write_alignment_stats(bam_file: pathlib.Path, out_file: pathlib.Path):
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
