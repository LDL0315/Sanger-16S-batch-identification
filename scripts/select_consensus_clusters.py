#!/usr/bin/env python3

import argparse
import os
import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Tuple


SIZE_PATTERN = re.compile(r";size=(\d+)(?:;|$)")


@dataclass
class FastaRecord:
    header: str
    sequence: str

    @property
    def label(self) -> str:
        return self.header.split()[0]

    @property
    def size(self) -> int:
        match = SIZE_PATTERN.search(self.header)
        return int(match.group(1)) if match else 0


def strip_size_suffix(label: str) -> str:
    return re.sub(r";size=\d+(?:;|$)", "", label)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Select consensus clusters by abundance thresholds.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--centroids", required=True)
    parser.add_argument("--uc", required=True)
    parser.add_argument("--total-reads", required=True, type=int)
    parser.add_argument("--min-reads", required=True, type=int)
    parser.add_argument("--min-frac", required=True, type=float)
    parser.add_argument("--passed-fasta", required=True)
    parser.add_argument("--split-dir", required=True)
    parser.add_argument("--abundance-tsv", required=True)
    return parser.parse_args()


def read_fasta(path: str) -> List[FastaRecord]:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return []

    records: List[FastaRecord] = []
    header = None
    seq_chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header=header, sequence="".join(seq_chunks)))
                header = line[1:]
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        records.append(FastaRecord(header=header, sequence="".join(seq_chunks)))
    return records


def parse_uc_counts(path: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return counts

    cluster_to_seed: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            record_type = fields[0]
            cluster_id = fields[1]
            query_label = fields[8]
            query_key = strip_size_suffix(query_label)
            query_size_match = SIZE_PATTERN.search(query_label)
            query_size = int(query_size_match.group(1)) if query_size_match else 1
            if record_type == "S":
                cluster_to_seed[cluster_id] = query_key
                counts.setdefault(query_key, 0)
                counts[query_key] += query_size
            elif record_type == "H":
                seed_label = cluster_to_seed.get(cluster_id)
                if seed_label is not None:
                    counts.setdefault(seed_label, 0)
                    counts[seed_label] += query_size
            elif record_type == "C":
                cluster_size = int(fields[2])
                seed_label = cluster_to_seed.get(cluster_id)
                if seed_label is not None:
                    counts[seed_label] = cluster_size
    return counts


def rename_passed_records(
    records: Iterable[FastaRecord],
    sample: str,
    total_reads: int,
    min_reads: int,
    min_frac: float,
    uc_counts: Dict[str, int],
) -> List[Tuple[str, FastaRecord, int, float]]:
    selected: List[Tuple[str, FastaRecord, int, float]] = []
    index = 1
    for record in records:
        cluster_reads = record.size or uc_counts.get(strip_size_suffix(record.label), 0)
        if total_reads > 0:
            cluster_frac = cluster_reads / total_reads
        else:
            cluster_frac = 0.0
        if cluster_reads < min_reads or cluster_frac < min_frac:
            continue
        consensus_id = f"consensus_{index:03d}"
        new_header = f"{consensus_id} sample={sample} reads={cluster_reads} fraction={cluster_frac:.6f} length={len(record.sequence)}"
        selected.append(
            (
                consensus_id,
                FastaRecord(header=new_header, sequence=record.sequence),
                cluster_reads,
                cluster_frac,
            )
        )
        index += 1
    return selected


def write_abundance(path: str, sample: str, selected: Iterable[Tuple[str, FastaRecord, int, float]]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("sample\tconsensus_id\tread_count\tread_fraction\tlength\n")
        for consensus_id, record, cluster_reads, cluster_frac in selected:
            handle.write(
                f"{sample}\t{consensus_id}\t{cluster_reads}\t{cluster_frac:.6f}\t{len(record.sequence)}\n"
            )


def write_passed_fasta(path: str, records: Iterable[Tuple[str, FastaRecord, int, float]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for _, record, _, _ in records:
            handle.write(f">{record.header}\n{record.sequence}\n")


def write_split_fastas(split_dir: str, records: Iterable[Tuple[str, FastaRecord, int, float]]) -> None:
    os.makedirs(split_dir, exist_ok=True)
    for existing in os.listdir(split_dir):
        if existing.startswith("consensus_") and existing.endswith(".fasta") and existing != "consensus_passed.fasta":
            os.remove(os.path.join(split_dir, existing))

    for consensus_id, record, _, _ in records:
        output_path = os.path.join(split_dir, f"{consensus_id}.fasta")
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write(f">{record.header}\n{record.sequence}\n")


def main() -> None:
    args = parse_args()
    os.makedirs(os.path.dirname(args.passed_fasta), exist_ok=True)
    os.makedirs(os.path.dirname(args.abundance_tsv), exist_ok=True)
    os.makedirs(args.split_dir, exist_ok=True)

    all_records = read_fasta(args.centroids)
    uc_counts = parse_uc_counts(args.uc)
    selected = rename_passed_records(
        records=all_records,
        sample=args.sample,
        total_reads=args.total_reads,
        min_reads=args.min_reads,
        min_frac=args.min_frac,
        uc_counts=uc_counts,
    )

    write_abundance(args.abundance_tsv, args.sample, selected)
    write_passed_fasta(args.passed_fasta, selected)
    write_split_fastas(args.split_dir, selected)


if __name__ == "__main__":
    main()
