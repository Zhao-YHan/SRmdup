#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pysam


def end_clip_len(read: pysam.AlignedSegment, which: str) -> int:
    ct = read.cigartuples or []
    if not ct:
        return 0

    ops = ct if which == "left" else ct[::-1]
    n = 0
    for op, ln in ops:
        if op in (pysam.CSOFT_CLIP, pysam.CHARD_CLIP):  # S or H
            n += ln
        else:
            break
    return n


def real_len(read: pysam.AlignedSegment) -> int:
    ct = read.cigartuples or []
    if not ct:
        return 0

    n = 0
    for op, ln in ct:
        if op in (
            pysam.CMATCH,        # M
            pysam.CINS,          # I
            pysam.CSOFT_CLIP,    # S
            pysam.CEQUAL,        # =
            pysam.CDIFF,         # X
            pysam.CHARD_CLIP     # H
        ):
            n += ln
    return n


def frag_len(read: pysam.AlignedSegment) -> int:
    if read.is_unmapped or read.reference_start is None or read.reference_end is None:
        return 0

    ref_span = read.reference_end - read.reference_start
    left_clip = end_clip_len(read, "left")
    right_clip = end_clip_len(read, "right")
    return ref_span + left_clip + right_clip


def main():
    if len(sys.argv) != 3:
        print("Usage: python reads_length.py input.bam output.tsv")
        sys.exit(1)

    bam_path = sys.argv[1]
    out_path = sys.argv[2]

    n = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam, open(out_path, "w", encoding="utf-8") as out:
        out.write("qname\tchr\tboundary_len\treal_len\tcigar\n")

        for a in bam:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue

            f_len = frag_len(a)
            r_len = real_len(a)

            if f_len <= 0 or r_len <= 0:
                continue

            chrom = bam.get_reference_name(a.reference_id)
            cigar = a.cigarstring or "*"

            out.write(f"{a.query_name}\t{chrom}\t{f_len}\t{r_len}\t{cigar}\n")
            n += 1

    print(f"Done. {n} reads.")


if __name__ == "__main__":
    main()