#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SRmdup: read-level duplicate filtering for single-end low-coverage ancient human DNA.

This script performs read-level duplicate filtering on BAM files using
softclip-corrected outer boundaries. Terminal softclips at both ends
are considered when defining outer coordinates.

The current implementation removes duplicates in four categories:
  1. same start and same end
  2. same start
  3. same end
  4. start and end within 1 bp

Winner selection is currently based on qsum, with real_len used as a
secondary tie-breaker.

Input:
  - BAM file before duplicate filtering

Main outputs:
  - filtered BAM file
  - a table with the count of deleted reads by strand and category

Optional outputs:
  - deleted-read tables
  - kept-read table
"""


import argparse
import sys
import os
import tempfile
import heapq
import struct
import zlib
from collections import defaultdict
from typing import Tuple, List, Callable, Optional
import pysam


# ----------------------------- Basic metrics -----------------------------

def flag01(a: pysam.AlignedSegment) -> int:
    return 16 if a.is_reverse else 0


def qsum(a: pysam.AlignedSegment) -> int:
    q = a.query_qualities
    return int(sum(q)) if q is not None else 0


def real_len(a: pysam.AlignedSegment) -> int:
    return a.query_length  # includes soft clips


def score_tuple(qsum_v: int, rlen_v: int) -> Tuple[int, int]:
    # higher is better
    return (qsum_v, rlen_v)


# ---------------------- Terminal softclips / outer coords ---------------

def terminal_softclips(read: pysam.AlignedSegment) -> Tuple[int, int]:
    """
    Return (left_softclip, right_softclip) in reference-coordinate orientation.

    Supports terminal H+S patterns such as:
      left side : 5H3S92M
      right side: 92M3S5H

    Only terminal softclips are considered; internal S operations are ignored.
    """
    ct = read.cigartuples or []
    if not ct:
        return 0, 0

    left_sc = 0
    i = 0
    while i < len(ct) and ct[i][0] == pysam.CHARD_CLIP:
        i += 1
    if i < len(ct) and ct[i][0] == pysam.CSOFT_CLIP:
        left_sc = ct[i][1]

    right_sc = 0
    j = len(ct) - 1
    while j >= 0 and ct[j][0] == pysam.CHARD_CLIP:
        j -= 1
    if j >= 0 and ct[j][0] == pysam.CSOFT_CLIP:
        right_sc = ct[j][1]

    return left_sc, right_sc


def outer_coords_1based(read: pysam.AlignedSegment):
    """
    Softclip-corrected outer boundaries using BOTH termini.

    start = (reference_start + 1) - left_terminal_softclip
    end   = reference_end + right_terminal_softclip
    """
    if read.is_unmapped:
        return None

    left_sc, right_sc = terminal_softclips(read)
    s1 = read.reference_start + 1
    e1 = read.reference_end

    start = max(1, s1 - left_sc)
    end = e1 + right_sc
    return start, end


def stable_key_outer_1based(read: pysam.AlignedSegment):
    oc = outer_coords_1based(read)
    if oc is None:
        return None
    s1, e1 = oc
    return (read.query_name, flag01(read), read.reference_id, s1, e1 + 1, read.cigarstring or "*")


# --------------------------- 64-bit key hashing --------------------------

def key64_from_fields(qname: str, flag: int, tid: int, s: int, e_plus1: int, cigar: str) -> int:
    b = f"{qname}\t{flag}\t{tid}\t{s}\t{e_plus1}\t{cigar}".encode("utf-8", "ignore")
    h = zlib.crc32(b)
    b2 = struct.pack("<II", h, len(b))
    h2 = zlib.crc32(b2, h)
    return (h << 32) | (h2 & 0xffffffff)


# Row = (qname, tid, flag, s, e, qsum, rlen, mapq, cigar)
Row = Tuple[str, int, int, int, int, int, int, int, str]


def key64_from_row(row: Row) -> int:
    qname, tid, flag, s, e, qsum_v, rlen_v, mapq, cigar = row
    return key64_from_fields(qname, flag, tid, s, e + 1, cigar or "*")


# -------------------------- External sort helpers ------------------------

def parse_line(line: str) -> Row:
    qname, tid, flag, s, e, qsum_v, rlen_v, mapq, cigar = line.rstrip("\n").split("\t")
    return (qname, int(tid), int(flag), int(s), int(e), int(qsum_v), int(rlen_v), int(mapq), cigar)


def row_to_line(r: Row) -> str:
    return "\t".join((
        r[0], str(r[1]), str(r[2]), str(r[3]), str(r[4]),
        str(r[5]), str(r[6]), str(r[7]), r[8]
    )) + "\n"


def chunk_sort_write(lines: List[str], keyfn: Callable[[Row], Tuple], tmpdir: str) -> str:
    rows = [parse_line(ln) for ln in lines]
    rows.sort(key=keyfn)
    fh = tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir, prefix="chunk_", suffix=".tsv")
    with fh as w:
        w.writelines(row_to_line(r) for r in rows)
    return fh.name


def ext_sort(inp_path: str, keyfn: Callable[[Row], Tuple], tmpdir: str,
             chunk_lines: int = 1_000_000, progress_every: int = 2_000_000) -> str:
    """External sort a TSV file of Row records by keyfn. Returns path to merged file."""
    chunk_files: List[str] = []
    total = 0

    with open(inp_path, "r", encoding="utf-8") as f:
        buf: List[str] = []
        for ln in f:
            buf.append(ln)
            total += 1
            if len(buf) >= chunk_lines:
                chunk_files.append(chunk_sort_write(buf, keyfn, tmpdir))
                buf = []
                if total % progress_every == 0:
                    print(f"\r  sorted {total} rows into {len(chunk_files)} chunks...", end="", flush=True)
        if buf:
            chunk_files.append(chunk_sort_write(buf, keyfn, tmpdir))

    print(f"\r  sorting complete: {total} rows, {len(chunk_files)} chunks.            ")

    if len(chunk_files) == 1:
        return chunk_files[0]

    out_path = tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir, prefix="merge_", suffix=".tsv").name
    fps = [open(p, "r", encoding="utf-8") for p in chunk_files]

    class Stream:
        __slots__ = ("fp", "keyfn", "row", "key")

        def __init__(self, fp):
            self.fp = fp
            self.keyfn = keyfn
            self.row = None
            self.key = None

        def read_next(self):
            ln = self.fp.readline()
            if not ln:
                self.row = None
                self.key = None
                return False
            self.row = parse_line(ln)
            self.key = self.keyfn(self.row)
            return True

    streams = []
    for fp in fps:
        st = Stream(fp)
        if st.read_next():
            streams.append(st)

    with open(out_path, "w", encoding="utf-8") as out:
        heap = [(st.key, i) for i, st in enumerate(streams)]
        heapq.heapify(heap)
        emitted = 0
        while heap:
            _, i = heapq.heappop(heap)
            st = streams[i]
            out.write(row_to_line(st.row))
            emitted += 1
            if not st.read_next():
                continue
            heapq.heappush(heap, (st.key, i))
            if emitted % progress_every == 0:
                print(f"\r  merged {emitted} rows...", end="", flush=True)

    print(f"\r  merge complete.                                     ")

    for fp in fps:
        fp.close()
    for p in chunk_files:
        try:
            os.unlink(p)
        except OSError:
            pass

    return out_path


# ----------------------- Dedup step implementations ----------------------

def best_of(rows: List[Row]) -> Row:
    b = rows[0]
    bs = score_tuple(b[5], b[6])  # qsum, rlen
    for r in rows[1:]:
        s = score_tuple(r[5], r[6])
        if s > bs:
            b, bs = r, s
    return b


def stepA_both_equal(inp_path: str, tmpdir: str, del_writer, delkey_writer, stats_flag: dict,
                     progress_every: int) -> str:
    print("  Step A: both_equal (same (s,e)) ...")
    sorted_path = ext_sort(inp_path, keyfn=lambda r: (r[3], r[4]), tmpdir=tmpdir, progress_every=progress_every)
    out_path = tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir, prefix="afterA_", suffix=".tsv").name

    try:
        with open(sorted_path, "r", encoding="utf-8") as f, open(out_path, "w", encoding="utf-8") as out:
            group_rows: List[Row] = []
            last_key = None  # (s, e)
            processed = 0

            for ln in f:
                r = parse_line(ln)
                k = (r[3], r[4])

                if last_key is None:
                    group_rows = [r]
                    last_key = k
                elif k == last_key:
                    group_rows.append(r)
                else:
                    if len(group_rows) > 1:
                        winner = best_of(group_rows)
                        for rr in group_rows:
                            if rr is not winner:
                                del_writer(rr)
                                delkey_writer(rr)
                                stats_flag["both_equal"] += 1
                        out.write(row_to_line(winner))
                    else:
                        out.write(row_to_line(group_rows[0]))

                    group_rows = [r]
                    last_key = k

                processed += 1
                if processed % progress_every == 0:
                    print(f"\r    processed {processed} rows...", end="", flush=True)

            if group_rows:
                if len(group_rows) > 1:
                    winner = best_of(group_rows)
                    for rr in group_rows:
                        if rr is not winner:
                            del_writer(rr)
                            delkey_writer(rr)
                            stats_flag["both_equal"] += 1
                    out.write(row_to_line(winner))
                else:
                    out.write(row_to_line(group_rows[0]))
    finally:
        try:
            os.unlink(sorted_path)
        except OSError:
            pass

    print("\r    Step A done.                                  ")
    return out_path


def stepB_start_equal(inp_path: str, tmpdir: str, del_writer, delkey_writer, stats_flag: dict,
                      global_best_hash: Optional[int], progress_every: int) -> str:
    print("  Step B: start_equal (same s) ...")
    sorted_path = ext_sort(inp_path, keyfn=lambda r: (r[3], r[4]), tmpdir=tmpdir, progress_every=progress_every)
    out_path = tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir, prefix="afterB_", suffix=".tsv").name

    try:
        with open(sorted_path, "r", encoding="utf-8") as f, open(out_path, "w", encoding="utf-8") as out:
            group_rows: List[Row] = []
            last_key = None  # s
            processed = 0

            for ln in f:
                r = parse_line(ln)
                k = r[3]

                if last_key is None:
                    group_rows = [r]
                    last_key = k
                elif k == last_key:
                    group_rows.append(r)
                else:
                    winner = None
                    if global_best_hash is not None:
                        for rr in group_rows:
                            if key64_from_row(rr) == global_best_hash:
                                winner = rr
                                break
                    if winner is None:
                        winner = best_of(group_rows)

                    for rr in group_rows:
                        if rr is not winner:
                            del_writer(rr)
                            delkey_writer(rr)
                            stats_flag["start_equal"] += 1
                    out.write(row_to_line(winner))

                    group_rows = [r]
                    last_key = k

                processed += 1
                if processed % progress_every == 0:
                    print(f"\r    processed {processed} rows...", end="", flush=True)

            if group_rows:
                winner = None
                if global_best_hash is not None:
                    for rr in group_rows:
                        if key64_from_row(rr) == global_best_hash:
                            winner = rr
                            break
                if winner is None:
                    winner = best_of(group_rows)

                for rr in group_rows:
                    if rr is not winner:
                        del_writer(rr)
                        delkey_writer(rr)
                        stats_flag["start_equal"] += 1
                out.write(row_to_line(winner))
    finally:
        try:
            os.unlink(sorted_path)
        except OSError:
            pass

    print("\r    Step B done.                                  ")
    return out_path


def stepC_end_equal(inp_path: str, tmpdir: str, del_writer, delkey_writer, stats_flag: dict,
                    global_best_hash: Optional[int], progress_every: int) -> str:
    print("  Step C: end_equal (same e) ...")
    sorted_path = ext_sort(inp_path, keyfn=lambda r: (r[4], r[3]), tmpdir=tmpdir, progress_every=progress_every)
    out_path = tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir, prefix="afterC_", suffix=".tsv").name

    try:
        with open(sorted_path, "r", encoding="utf-8") as f, open(out_path, "w", encoding="utf-8") as out:
            group_rows: List[Row] = []
            last_key = None  # e
            processed = 0

            for ln in f:
                r = parse_line(ln)
                k = r[4]

                if last_key is None:
                    group_rows = [r]
                    last_key = k
                elif k == last_key:
                    group_rows.append(r)
                else:
                    winner = None
                    if global_best_hash is not None:
                        for rr in group_rows:
                            if key64_from_row(rr) == global_best_hash:
                                winner = rr
                                break
                    if winner is None:
                        winner = best_of(group_rows)

                    for rr in group_rows:
                        if rr is not winner:
                            del_writer(rr)
                            delkey_writer(rr)
                            stats_flag["end_equal"] += 1
                    out.write(row_to_line(winner))

                    group_rows = [r]
                    last_key = k

                processed += 1
                if processed % progress_every == 0:
                    print(f"\r    processed {processed} rows...", end="", flush=True)

            if group_rows:
                winner = None
                if global_best_hash is not None:
                    for rr in group_rows:
                        if key64_from_row(rr) == global_best_hash:
                            winner = rr
                            break
                if winner is None:
                    winner = best_of(group_rows)

                for rr in group_rows:
                    if rr is not winner:
                        del_writer(rr)
                        delkey_writer(rr)
                        stats_flag["end_equal"] += 1
                out.write(row_to_line(winner))
    finally:
        try:
            os.unlink(sorted_path)
        except OSError:
            pass

    print("\r    Step C done.                                  ")
    return out_path


def stepD_diag_neighbors(inp_path: str, tmpdir: str, del_writer, delkey_writer, stats_flag: dict,
                         global_best_hash: Optional[int], progress_every: int) -> str:
    print("  Step D: both_±1 ...")

    def _run_pass(in_path: str, keyfn, pass_name: str) -> str:
        sorted_path = ext_sort(in_path, keyfn=keyfn, tmpdir=tmpdir, progress_every=progress_every)
        out_path = tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir,
                                               prefix=f"afterD_{pass_name}_", suffix=".tsv").name

        try:
            with open(sorted_path, "r", encoding="utf-8") as f, open(out_path, "w", encoding="utf-8") as out:
                current_key = None
                run: List[Row] = []
                last_s = None
                processed = 0

                def flush_run(rrun: List[Row]):
                    if not rrun:
                        return

                    winner = None
                    if global_best_hash is not None:
                        for rr in rrun:
                            if key64_from_row(rr) == global_best_hash:
                                winner = rr
                                break
                    if winner is None:
                        winner = best_of(rrun)

                    for rr in rrun:
                        if rr is not winner:
                            del_writer(rr)
                            delkey_writer(rr)
                            stats_flag["both_±1"] += 1
                    out.write(row_to_line(winner))

                for ln in f:
                    r = parse_line(ln)
                    s = r[3]
                    gkey = keyfn(r)
                    key0 = gkey[:2]

                    if current_key is None:
                        current_key = key0
                        run = [r]
                        last_s = s
                    elif key0 != current_key:
                        flush_run(run)
                        current_key = key0
                        run = [r]
                        last_s = s
                    else:
                        if last_s is not None and s == last_s + 1:
                            run.append(r)
                            last_s = s
                        else:
                            flush_run(run)
                            run = [r]
                            last_s = s

                    processed += 1
                    if processed % progress_every == 0:
                        print(f"\r    [{pass_name}] processed {processed} rows...", end="", flush=True)

                flush_run(run)
        finally:
            try:
                os.unlink(sorted_path)
            except OSError:
                pass

        print(f"\r    [{pass_name}] pass done.                                  ")
        return out_path

    def keyfn_diag(r: Row):
        return (r[1], r[3] - r[4], r[3])

    def keyfn_anti(r: Row):
        return (r[1], r[3] + r[4], r[3])

    out1 = _run_pass(inp_path, keyfn_diag, "diag")
    out2 = _run_pass(out1, keyfn_anti, "anti")

    try:
        os.unlink(out1)
    except OSError:
        pass

    print("  Step D done.                                  ")
    return out2


# ------------------------------- Main flow -------------------------------

def main():
    ap = argparse.ArgumentParser(description="Low-memory BAM dedup with external sorting and progress output")
    ap.add_argument("--bam", required=True, help="Input BAM (coordinate-sorted recommended)")
    ap.add_argument("--out-bam", required=True, help="Output BAM with duplicates removed")
    ap.add_argument("--stats-tsv", required=True, help="Summary stats TSV")
    ap.add_argument("--del-prefix", default=None, help="Write deleted read details to <prefix>.forward.tsv / .reverse.tsv")
    ap.add_argument("--keep-tsv", default=None, help="(Optional) TSV of kept reads written in pass-2")
    ap.add_argument("--tmp-dir", default=None, help="Directory for temporary files (defaults to system temp)")
    ap.add_argument("--read-threads", type=int, default=4, help="BGZF read threads for pysam/htslib")
    ap.add_argument("--write-threads", type=int, default=4, help="BGZF write threads for pysam/htslib")
    ap.add_argument("--progress-every", type=int, default=100000, help="Print progress every N rows")
    ap.add_argument("--accept-unsorted", action="store_true", help="Allow unsorted BAM (still works; slower)")

    args = ap.parse_args()

    tmpdir = args.tmp_dir or tempfile.mkdtemp(prefix="dedup_tmp_")
    os.makedirs(tmpdir, exist_ok=True)
    print(f"[tmp] Using temp dir: {tmpdir}") 

    stats_flag = {0: defaultdict(int), 16: defaultdict(int)}
    input_counts = {0: 0, 16: 0}

    del_writer = {0: (lambda _r, _rn, _cat: None), 16: (lambda _r, _rn, _cat: None)}
    del_fhs = {0: None, 16: None}

    if args.del_prefix:
        for f01 in (0, 16):
            path = f"{args.del_prefix}.{'forward' if f01 == 0 else 'reverse'}.tsv"
            fh = open(path, "w", encoding="utf-8")
            fh.write("\t".join([
                "qname", "rname", "start_outer_1based", "end_outer_1based",
                "flag", "qual_sum", "real_len", "mapq", "cigar", "category"
            ]) + "\n")

            def make_w(fh_local, flag_val):
                def _w(row: Row, rname: str, cat: str):
                    fh_local.write("\t".join(map(str, [
                        row[0], rname, row[3], row[4], flag_val, row[5], row[6], row[7], row[8] or "*", cat
                    ])) + "\n")
                return _w

            del_writer[f01] = make_w(fh, f01)
            del_fhs[f01] = fh

    print("[1/3] Scanning BAM and spilling to per-(tid,flag) buckets ...")

    buckets = defaultdict(lambda: {
        0: tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir, prefix="bk_", suffix="_forward.tsv"),
        16: tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir, prefix="bk_", suffix="_reverse.tsv")
    })
    bucket_best = defaultdict(lambda: {0: (None, (-1, -1)), 16: (None, (-1, -1))})

    pos_counter = 0
    with pysam.AlignmentFile(args.bam, "rb", threads=max(1, args.read_threads)) as bam:
        try:
            so = bam.header.to_dict().get("HD", {}).get("SO", None)
        except Exception:
            so = None
        if (so != "coordinate") and (not args.accept_unsorted):
            print("[WARN] Input BAM not coordinate-sorted (HD:SO:coordinate missing). Consider sorting or pass --accept-unsorted.")

        for a in bam:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue

            oc = outer_coords_1based(a)
            if oc is None:
                continue
            s1, e1 = oc
            if e1 < s1:
                continue

            tid = a.reference_id
            f01 = flag01(a)
            input_counts[f01] += 1

            row: Row = (
                a.query_name, tid, f01, s1, e1,
                qsum(a), real_len(a), a.mapping_quality, a.cigarstring or "*"
            )
            buckets[tid][f01].write(row_to_line(row))

            sc = score_tuple(row[5], row[6])
            best_row, best_sc = bucket_best[tid][f01]
            if sc > best_sc:
                bucket_best[tid][f01] = (row, sc)

            pos_counter += 1
            if pos_counter % args.progress_every == 0:
                print(f"\r  {pos_counter} reads processed...", end="", flush=True)

    print(f"\r  BAM scan complete. Total reads processed: {pos_counter}           ")

    delkey_paths = {}
    print("[2/3] Deduplicating per contig ...")

    for tid, flag_files in buckets.items():
        with pysam.AlignmentFile(args.bam, "rb") as bam_hdr:
            rname = bam_hdr.get_reference_name(tid)

        delkey_path = tempfile.NamedTemporaryFile("w+", delete=False, dir=tmpdir, prefix=f"delkey_tid{tid}_", suffix=".txt").name
        delkey_paths[tid] = delkey_path
        delkey_fh = open(delkey_path, "w", encoding="utf-8")

        for f01 in (0, 16):
            fh = flag_files[f01]
            fh.flush()
            fh.seek(0)
            raw_path = fh.name
            fh.close()

            def del_detail_writer(row: Row, cat: str):
                del_writer[f01](row, rname, cat)

            def delkey_writer(row: Row):
                delkey_fh.write(str(key64_from_row(row)) + "\n")

            best_row, _ = bucket_best[tid][f01]
            gbest_hash = key64_from_row(best_row) if best_row is not None else None

            pathA = stepA_both_equal(
                raw_path, tmpdir,
                del_writer=lambda rr: del_detail_writer(rr, "both_equal"),
                delkey_writer=delkey_writer,
                stats_flag=stats_flag[f01],
                progress_every=args.progress_every
            )
            try:
                os.unlink(raw_path)
            except OSError:
                pass

            pathB = stepB_start_equal(
                pathA, tmpdir,
                del_writer=lambda rr: del_detail_writer(rr, "start_equal"),
                delkey_writer=delkey_writer,
                stats_flag=stats_flag[f01],
                global_best_hash=gbest_hash,
                progress_every=args.progress_every
            )
            try:
                os.unlink(pathA)
            except OSError:
                pass

            pathC = stepC_end_equal(
                pathB, tmpdir,
                del_writer=lambda rr: del_detail_writer(rr, "end_equal"),
                delkey_writer=delkey_writer,
                stats_flag=stats_flag[f01],
                global_best_hash=gbest_hash,
                progress_every=args.progress_every
            )
            try:
                os.unlink(pathB)
            except OSError:
                pass

            pathD = stepD_diag_neighbors(
                pathC, tmpdir,
                del_writer=lambda rr: del_detail_writer(rr, "both_±1"),
                delkey_writer=delkey_writer,
                stats_flag=stats_flag[f01],
                global_best_hash=gbest_hash,
                progress_every=args.progress_every
            )
            try:
                os.unlink(pathC)
            except OSError:
                pass
            try:
                os.unlink(pathD)
            except OSError:
                pass

        delkey_fh.close()

    print("[3/3] Writing kept BAM (streaming, per-contig deletion set) ...")
    keep_fh = open(args.keep_tsv, "w", encoding="utf-8") if args.keep_tsv else None
    if keep_fh:
        keep_fh.write("\t".join([
            "qname", "rname", "start_outer_1based", "end_outer_1based",
            "flag", "qual_sum", "real_len", "mapq", "cigar"
        ]) + "\n")

    total_written = 0

    with pysam.AlignmentFile(args.bam, "rb", threads=max(1, args.read_threads)) as bam2, \
         pysam.AlignmentFile(args.out_bam, "wb", header=bam2.header, threads=max(1, args.write_threads)) as kept:

        tid2name = {i: bam2.get_reference_name(i) for i in range(len(bam2.header.references))}
        current_tid = None
        delset = set()

        def load_delset_for_tid(tid_val: int):
            delset.clear()
            p = delkey_paths.get(tid_val, None)
            if not p or (not os.path.exists(p)):
                return
            with open(p, "r", encoding="utf-8") as f:
                for ln in f:
                    s = ln.strip()
                    if s:
                        delset.add(int(s))

        align_counter = 0
        for a in bam2:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue

            tid = a.reference_id
            if current_tid is None:
                current_tid = tid
                load_delset_for_tid(current_tid)
            elif tid != current_tid:
                current_tid = tid
                load_delset_for_tid(current_tid)

            k = stable_key_outer_1based(a)
            if k is None:
                continue

            h = key64_from_fields(k[0], k[1], k[2], k[3], k[4], k[5])
            if h in delset:
                continue

            kept.write(a)
            total_written += 1

            if keep_fh:
                s1, e1 = outer_coords_1based(a)
                keep_fh.write("\t".join(map(str, [
                    a.query_name, tid2name[tid], s1, e1, flag01(a), qsum(a), real_len(a),
                    a.mapping_quality, a.cigarstring or "*"
                ])) + "\n")

            align_counter += 1
            if align_counter % args.progress_every == 0:
                print(f"\r  wrote {align_counter} kept reads...", end="", flush=True)

    if keep_fh:
        keep_fh.close()
    print(f"\r  kept BAM writing complete. Total kept: {total_written}                 ")

    # Stats
    with open(args.stats_tsv, "w", encoding="utf-8") as sf:
        sf.write("\t".join([
            "scope", "total_input", "total_deleted", "total_kept", "kept_rate%",
            "both_equal", "start_equal", "end_equal", "both_±1"
        ]) + "\n")

        def get_counts(st):
            return (
                st.get("both_equal", 0),
                st.get("start_equal", 0),
                st.get("end_equal", 0),
                st.get("both_±1", 0),
            )

        def row(scope, inf, both_equal, start_equal, end_equal, both_pm1):
            deleted = both_equal + start_equal + end_equal + both_pm1
            keptn = inf - deleted
            kept_rate = (keptn / inf * 100.0) if inf > 0 else 0.0
            return [
                scope, str(inf), str(deleted), str(keptn), f"{kept_rate:.3f}",
                str(both_equal), str(start_equal), str(end_equal), str(both_pm1)
            ]

        f_both, f_start, f_end, f_pm1 = get_counts(stats_flag[0])
        r_both, r_start, r_end, r_pm1 = get_counts(stats_flag[16])

        rev_disp_start = r_end
        rev_disp_end = r_start

        all_both = f_both + r_both
        all_start = f_start + r_end
        all_end = f_end + r_start
        all_pm1 = f_pm1 + r_pm1

        all_counts = input_counts[0] + input_counts[16]

        sf.write("\t".join(row("FORWARD", input_counts[0], f_both, f_start, f_end, f_pm1)) + "\n")
        sf.write("\t".join(row("REVERSE", input_counts[16], r_both, rev_disp_start, rev_disp_end, r_pm1)) + "\n")
        sf.write("\t".join(row("ALL", all_counts, all_both, all_start, all_end, all_pm1)) + "\n")

    for p in delkey_paths.values():
        try:
            os.unlink(p)
        except OSError:
            pass

    for fh in del_fhs.values():
        if fh:
            fh.close()

    print("Done.")


if __name__ == "__main__":
    sys.exit(main())