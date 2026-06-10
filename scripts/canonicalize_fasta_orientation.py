#!/usr/bin/env python3

import argparse
from typing import Iterable, Iterator, Tuple


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
DEFAULT_FORWARD_PRIMER = "AGAGTTTGATCCTGGCTCAG"
DEFAULT_REVERSE_PRIMER = "ACGGCTACCTTGTTACGACTT"
DEFAULT_WINDOW = 80


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize each FASTA sequence to the 27F->1492R amplicon orientation before clustering."
    )
    parser.add_argument("--input", required=True, help="Input FASTA path")
    parser.add_argument("--output", required=True, help="Output FASTA path")
    parser.add_argument("--forward-primer", default=DEFAULT_FORWARD_PRIMER)
    parser.add_argument("--reverse-primer", default=DEFAULT_REVERSE_PRIMER)
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    return parser.parse_args()


def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    header = None
    seq_parts = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if header is not None:
        yield header, "".join(seq_parts)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def best_primer_match(primer: str, region: str) -> int:
    primer = primer.upper()
    region = region.upper()
    if not primer or not region:
        return 0
    if len(region) < len(primer):
        compare_len = len(region)
        return sum(1 for idx in range(compare_len) if region[idx] == primer[idx])

    best = 0
    limit = len(region) - len(primer) + 1
    for start in range(limit):
        score = sum(1 for idx, base in enumerate(primer) if region[start + idx] == base)
        if score > best:
            best = score
    return best


def orientation_score(sequence: str, forward_primer: str, reverse_primer_rc: str, window: int) -> int:
    head = sequence[:window]
    tail = sequence[-window:] if len(sequence) > window else sequence
    return best_primer_match(forward_primer, head) + best_primer_match(reverse_primer_rc, tail)


def canonicalize(sequence: str, forward_primer: str, reverse_primer: str, window: int) -> str:
    forward = sequence.upper()
    reverse = reverse_complement(forward)
    reverse_primer_rc = reverse_complement(reverse_primer.upper())

    forward_score = orientation_score(forward, forward_primer.upper(), reverse_primer_rc, window)
    reverse_score = orientation_score(reverse, forward_primer.upper(), reverse_primer_rc, window)

    if forward_score > reverse_score:
        return forward
    if reverse_score > forward_score:
        return reverse
    return forward if forward <= reverse else reverse


def write_fasta(path: str, records: Iterable[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


def main() -> None:
    args = parse_args()
    write_fasta(
        args.output,
        (
            (
                header,
                canonicalize(
                    sequence,
                    forward_primer=args.forward_primer,
                    reverse_primer=args.reverse_primer,
                    window=args.window,
                ),
            )
            for header, sequence in read_fasta(args.input)
        ),
    )


if __name__ == "__main__":
    main()
