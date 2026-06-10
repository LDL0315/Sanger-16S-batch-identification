#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


DEFAULT_SNP_FRAC = 0.2
DEFAULT_MIN_DEPTH = 20
DEFAULT_REPORT_MIN_FRAC = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a user-facing consensus report across all processed samples."
    )
    parser.add_argument("--output-root", required=True, help="Pipeline output root directory")
    parser.add_argument("--report-csv", required=True, help="Destination report CSV file")
    parser.add_argument(
        "--snp-frac-threshold",
        type=float,
        default=DEFAULT_SNP_FRAC,
        help=f"Count a position as SNP when nonref_frac >= this value (default: {DEFAULT_SNP_FRAC})",
    )
    parser.add_argument(
        "--min-depth",
        type=int,
        default=DEFAULT_MIN_DEPTH,
        help=f"Only count SNP positions with depth >= this value (default: {DEFAULT_MIN_DEPTH})",
    )
    parser.add_argument(
        "--report-min-frac",
        type=float,
        default=DEFAULT_REPORT_MIN_FRAC,
        help=f"Only include consensus rows with read_fraction >= this value in the final report (default: {DEFAULT_REPORT_MIN_FRAC})",
    )
    return parser.parse_args()


def read_abundance(path: Path) -> List[Dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_snp_positions(path: Path, snp_frac_threshold: float, min_depth: int) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0

    count = 0
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                depth = int(row["depth"])
                nonref_frac = float(row["nonref_frac"])
            except (KeyError, TypeError, ValueError):
                continue
            if depth >= min_depth and nonref_frac >= snp_frac_threshold:
                count += 1
    return count


def parse_consensus_number(consensus_id: str) -> str:
    if "_" not in consensus_id:
        return consensus_id
    suffix = consensus_id.rsplit("_", 1)[-1]
    try:
        return str(int(suffix))
    except ValueError:
        return consensus_id


def sample_display_name(sample_dir_name: str) -> str:
    return sample_dir_name.split("_", 1)[0]


def build_rows(
    output_root: Path, snp_frac_threshold: float, min_depth: int, report_min_frac: float
) -> List[Tuple[str, str, int, str, str]]:
    rows: List[Tuple[str, str, int, str, str]] = []
    for sample_dir in sorted(output_root.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name == "all_sample_fastas":
            continue

        abundance_path = sample_dir / "03_consensus" / "consensus_abundance.tsv"
        abundance_rows = read_abundance(abundance_path)
        if not abundance_rows:
            continue

        sample_name = sample_display_name(sample_dir.name)
        kept_rows = [row for row in abundance_rows if float(row["read_fraction"]) >= report_min_frac]
        if not kept_rows:
            kept_rows = abundance_rows
        kept_total_fraction = sum(float(row["read_fraction"]) for row in kept_rows)
        single_consensus = len(kept_rows) == 1

        for index, row in enumerate(kept_rows, start=1):
            consensus_id = row["consensus_id"]
            snp_path = sample_dir / "05_snp" / f"{consensus_id}.snp.tsv"
            snp_count = count_snp_positions(snp_path, snp_frac_threshold, min_depth)

            if snp_count == 0 and single_consensus:
                result_type = "\u5355\u4e00\u5e8f\u5217"
            else:
                result_type = f"SNP [{snp_count}]"

            rows.append(
                (
                    sample_name,
                    str(index),
                    int(row["length"]),
                    result_type,
                    f"{(float(row['read_fraction']) / kept_total_fraction) if kept_total_fraction > 0 else 0.0:.2f}",
                )
            )
    return rows


def write_report(path: Path, rows: Iterable[Tuple[str, str, int, str, str]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "\u6837\u54c1\u540d\u79f0",
                "\u7ec4\u88c5\u7ed3\u679c",
                "\u7ec4\u88c5\u957f\u5ea6",
                "\u7ed3\u679c\u7c7b\u578b",
                "\u5e8f\u5217\u5360\u6bd4",
            ]
        )
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    output_root = Path(args.output_root)
    report_path = Path(args.report_csv)
    report_path.parent.mkdir(parents=True, exist_ok=True)

    rows = build_rows(output_root, args.snp_frac_threshold, args.min_depth, args.report_min_frac)
    write_report(report_path, rows)


if __name__ == "__main__":
    main()
