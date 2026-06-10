#!/usr/bin/env python3

import argparse
import os
from typing import Dict, TextIO


HEADER = (
    "sample\tconsensus_id\tpos\tref\tdepth\tA_count\tC_count\tG_count\tT_count\t"
    "A_frac\tC_frac\tG_frac\tT_frac\tnonref_frac\n"
)
BASES = ("A", "C", "G", "T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Convert mpileup output into per-reference SNP frequency tables.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--pileup", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def count_bases(read_bases: str, ref_base: str) -> Dict[str, int]:
    counts = {base: 0 for base in BASES}
    ref_base = ref_base.upper()
    index = 0
    while index < len(read_bases):
        char = read_bases[index]

        if char == "^":
            index += 2
            continue
        if char == "$":
            index += 1
            continue
        if char in "+-":
            index += 1
            length_start = index
            while index < len(read_bases) and read_bases[index].isdigit():
                index += 1
            if length_start == index:
                continue
            indel_length = int(read_bases[length_start:index])
            index += indel_length
            continue
        if char in ".,":
            if ref_base in counts:
                counts[ref_base] += 1
            index += 1
            continue

        base = char.upper()
        if base in counts:
            counts[base] += 1
        index += 1

    return counts


def open_writer(output_dir: str, consensus_id: str) -> TextIO:
    output_path = os.path.join(output_dir, f"{consensus_id}.snp.tsv")
    handle = open(output_path, "w", encoding="utf-8", newline="")
    handle.write(HEADER)
    return handle


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    writers: Dict[str, TextIO] = {}
    try:
        with open(args.pileup, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) < 5:
                    continue

                consensus_id = fields[0]
                pos = fields[1]
                ref_base = fields[2].upper()
                depth = int(fields[3])
                read_bases = fields[4] if len(fields) >= 5 else ""
                counts = count_bases(read_bases, ref_base)

                if depth > 0:
                    freqs = {base: counts[base] / depth for base in BASES}
                    nonref = sum(counts[base] for base in BASES if base != ref_base) / depth
                else:
                    freqs = {base: 0.0 for base in BASES}
                    nonref = 0.0

                writer = writers.get(consensus_id)
                if writer is None:
                    writer = open_writer(args.output_dir, consensus_id)
                    writers[consensus_id] = writer

                writer.write(
                    f"{args.sample}\t{consensus_id}\t{pos}\t{ref_base}\t{depth}\t"
                    f"{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}\t"
                    f"{freqs['A']:.6f}\t{freqs['C']:.6f}\t{freqs['G']:.6f}\t{freqs['T']:.6f}\t{nonref:.6f}\n"
                )
    finally:
        for writer in writers.values():
            writer.close()


if __name__ == "__main__":
    main()
