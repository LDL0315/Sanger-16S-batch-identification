#!/usr/bin/env python3
"""Batch identify Sanger 16S FASTA sequences with a local NCBI BLAST database."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


DB_NAME = "16S_ribosomal_RNA"
FASTA_SUFFIXES = {".fasta", ".fa", ".fas", ".fna"}
VALID_DNA = set("ACGTRYSWKMBDHVN")

CSV_FIELDS = [
    "sample_file",
    "query_id",
    "query_length",
    "best_accession",
    "best_title",
    "organism_name",
    "percent_identity",
    "query_coverage",
    "alignment_length",
    "mismatches",
    "gaps",
    "evalue",
    "bitscore",
    "confidence",
    "note",
]


@dataclass(frozen=True)
class QueryRecord:
    internal_id: str
    sample_file: str
    query_id: str
    sequence: str
    note: str = ""

    @property
    def query_length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class BlastHit:
    internal_id: str
    accession: str
    title: str
    scientific_name: str
    percent_identity: float
    query_coverage: float
    alignment_length: int
    mismatches: int
    gaps: int
    evalue: float
    bitscore: float
    query_length: int


def eprint(message: str) -> None:
    print(message, file=sys.stderr)


def resolve_db_prefix(db_dir: Path) -> Path:
    return db_dir / DB_NAME


def database_exists(db_dir: Path) -> bool:
    prefix = resolve_db_prefix(db_dir)
    required_groups = [
        (".nhr", ".nin", ".nsq"),
        (".ndb", ".ntf", ".nto"),
    ]
    return any(all((prefix.with_suffix(ext)).exists() for ext in group) for group in required_groups)


def require_executable(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise RuntimeError(f"Required executable not found in PATH: {name}")
    return path


def discover_fasta_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        if input_path.suffix.lower() not in FASTA_SUFFIXES:
            raise ValueError(f"Input file does not look like FASTA: {input_path}")
        return [input_path]
    if not input_path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")
    return sorted(
        path
        for path in input_path.rglob("*")
        if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES
    )


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:].strip() or path.stem
                chunks = []
            elif header is None:
                # Keep a useful record for non-standard files that omit the first header.
                header = path.stem
                chunks = [line]
            else:
                chunks.append(line)

    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def sanitize_sequence(sequence: str) -> tuple[str, str]:
    cleaned = re.sub(r"\s+", "", sequence).upper()
    notes: list[str] = []
    if "U" in cleaned:
        cleaned = cleaned.replace("U", "T")
        notes.append("converted U to T")

    invalid_count = sum(1 for base in cleaned if base not in VALID_DNA)
    if invalid_count:
        cleaned = "".join(base if base in VALID_DNA else "N" for base in cleaned)
        notes.append(f"replaced {invalid_count} invalid character(s) with N")
    return cleaned, "; ".join(notes)


def load_queries(input_path: Path) -> tuple[list[QueryRecord], list[dict[str, str]]]:
    fasta_files = discover_fasta_files(input_path)
    if not fasta_files:
        raise ValueError(f"No FASTA files found in: {input_path}")

    queries: list[QueryRecord] = []
    skipped_rows: list[dict[str, str]] = []
    counter = 1

    for fasta_path in fasta_files:
        records = parse_fasta(fasta_path)
        sample_file = fasta_path.name
        if not records:
            skipped_rows.append(empty_result_row(sample_file, "", 0, "empty file or no FASTA records"))
            continue

        for header, raw_sequence in records:
            sequence, note = sanitize_sequence(raw_sequence)
            if not sequence:
                skipped_rows.append(empty_result_row(sample_file, header, 0, "empty sequence"))
                continue
            internal_id = f"q{counter:06d}"
            counter += 1
            queries.append(QueryRecord(internal_id, sample_file, header, sequence, note))

    return queries, skipped_rows


def write_query_fasta(records: Sequence[QueryRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.internal_id}\n")
            for index in range(0, len(record.sequence), 80):
                handle.write(record.sequence[index : index + 80] + "\n")


def run_command(args: Sequence[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        list(args),
        cwd=str(cwd) if cwd else None,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def run_blast(query_fasta: Path, db_dir: Path, max_target_seqs: int, threads: int) -> list[BlastHit]:
    db_prefix = resolve_db_prefix(db_dir)
    outfmt = "6 qseqid sacc stitle sscinames pident qcovs length mismatch gaps evalue bitscore qlen"
    cmd = [
        "blastn",
        "-task",
        "megablast",
        "-query",
        str(query_fasta),
        "-db",
        str(db_prefix),
        "-outfmt",
        outfmt,
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(threads),
    ]
    result = run_command(cmd)
    if result.returncode != 0:
        raise RuntimeError(
            "blastn failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"STDERR:\n{result.stderr.strip()}"
        )
    return parse_blast_table(result.stdout)


def parse_blast_table(table: str) -> list[BlastHit]:
    hits: list[BlastHit] = []
    for line_number, line in enumerate(table.splitlines(), start=1):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 12:
            raise ValueError(f"Unexpected BLAST outfmt line {line_number}: expected 12 fields, got {len(parts)}")
        try:
            hits.append(
                BlastHit(
                    internal_id=parts[0],
                    accession=parts[1],
                    title=parts[2],
                    scientific_name=parts[3],
                    percent_identity=float(parts[4]),
                    query_coverage=float(parts[5]),
                    alignment_length=int(parts[6]),
                    mismatches=int(parts[7]),
                    gaps=int(parts[8]),
                    evalue=float(parts[9]),
                    bitscore=float(parts[10]),
                    query_length=int(parts[11]),
                )
            )
        except ValueError as exc:
            raise ValueError(f"Could not parse BLAST outfmt line {line_number}: {line}") from exc
    return hits


def best_hit_key(hit: BlastHit) -> tuple[float, float, float, float]:
    return (hit.percent_identity, hit.query_coverage, -hit.evalue, hit.bitscore)


def group_best_hits(hits: Iterable[BlastHit]) -> dict[str, list[BlastHit]]:
    grouped: dict[str, list[BlastHit]] = {}
    for hit in hits:
        grouped.setdefault(hit.internal_id, []).append(hit)
    for internal_id, hit_list in grouped.items():
        grouped[internal_id] = sorted(hit_list, key=best_hit_key, reverse=True)
    return grouped


def confidence_for(hit: BlastHit) -> str:
    if hit.percent_identity >= 99.0 and hit.query_coverage >= 95.0:
        return "species_high"
    if hit.percent_identity >= 97.0 and hit.query_coverage >= 90.0:
        return "genus_or_close_species"
    return "low_confidence"


def organism_from_hit(hit: BlastHit) -> str:
    if hit.scientific_name and hit.scientific_name != "N/A":
        return hit.scientific_name
    title = hit.title
    # NCBI 16S titles commonly start with the organism name followed by strain/gene text.
    match = re.match(r"^([A-Z][a-zA-Z0-9_.-]+(?:\s+[a-z][a-zA-Z0-9_.-]+)?)\b", title)
    if match:
        return match.group(1)
    return title


def close_hit_note(best_hit: BlastHit, hits: Sequence[BlastHit]) -> str:
    close_hits = [
        hit
        for hit in hits[1:]
        if abs(hit.percent_identity - best_hit.percent_identity) <= 0.1
        and abs(hit.query_coverage - best_hit.query_coverage) <= 1.0
        and organism_from_hit(hit) != organism_from_hit(best_hit)
    ]
    if close_hits:
        return "multiple close hits; review manually"
    return ""


def result_row_for_query(record: QueryRecord, hits: Sequence[BlastHit]) -> dict[str, str]:
    if not hits:
        return empty_result_row(record.sample_file, record.query_id, record.query_length, combine_notes(record.note, "no BLAST hit"))

    best = hits[0]
    note = combine_notes(record.note, close_hit_note(best, hits))
    return {
        "sample_file": record.sample_file,
        "query_id": record.query_id,
        "query_length": str(record.query_length),
        "best_accession": best.accession,
        "best_title": best.title,
        "organism_name": organism_from_hit(best),
        "percent_identity": f"{best.percent_identity:.3f}",
        "query_coverage": f"{best.query_coverage:.1f}",
        "alignment_length": str(best.alignment_length),
        "mismatches": str(best.mismatches),
        "gaps": str(best.gaps),
        "evalue": f"{best.evalue:.3g}",
        "bitscore": f"{best.bitscore:.1f}",
        "confidence": confidence_for(best),
        "note": note,
    }


def combine_notes(*notes: str) -> str:
    return "; ".join(note for note in notes if note)


def empty_result_row(sample_file: str, query_id: str, query_length: int, note: str) -> dict[str, str]:
    return {
        "sample_file": sample_file,
        "query_id": query_id,
        "query_length": str(query_length),
        "best_accession": "",
        "best_title": "",
        "organism_name": "",
        "percent_identity": "",
        "query_coverage": "",
        "alignment_length": "",
        "mismatches": "",
        "gaps": "",
        "evalue": "",
        "bitscore": "",
        "confidence": "no_result",
        "note": note,
    }


def write_csv(rows: Sequence[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def command_check(args: argparse.Namespace) -> int:
    db_dir = Path(args.db_dir)
    ok = True

    for executable in ("blastn", "update_blastdb.pl"):
        path = shutil.which(executable)
        if path:
            print(f"OK: {executable} -> {path}")
        else:
            print(f"MISSING: {executable} not found in PATH")
            ok = False

    if database_exists(db_dir):
        print(f"OK: BLAST database found -> {resolve_db_prefix(db_dir)}")
    else:
        print(f"MISSING: BLAST database not found -> {resolve_db_prefix(db_dir)}")
        print("Run: python sanger16s.py setup-db --db-dir db")
        ok = False

    return 0 if ok else 1


def command_setup_db(args: argparse.Namespace) -> int:
    update_script = require_executable("update_blastdb.pl")
    db_dir = Path(args.db_dir)
    db_dir.mkdir(parents=True, exist_ok=True)
    cmd = [update_script, "--decompress", DB_NAME]
    eprint(f"Downloading {DB_NAME} into {db_dir} ...")
    result = run_command(cmd, cwd=db_dir)
    if result.stdout.strip():
        print(result.stdout, end="")
    if result.stderr.strip():
        eprint(result.stderr.strip())
    if result.returncode != 0:
        eprint("Database download failed.")
        return result.returncode
    if not database_exists(db_dir):
        eprint(f"Download finished, but database files were not detected at {resolve_db_prefix(db_dir)}")
        return 1
    print(f"OK: database ready at {resolve_db_prefix(db_dir)}")
    return 0


def command_run(args: argparse.Namespace) -> int:
    require_executable("blastn")
    input_path = Path(args.input)
    db_dir = Path(args.db_dir)
    out_path = Path(args.out)

    if not database_exists(db_dir):
        raise RuntimeError(f"BLAST database not found at {resolve_db_prefix(db_dir)}. Run setup-db first.")

    queries, skipped_rows = load_queries(input_path)
    rows: list[dict[str, str]] = list(skipped_rows)

    if queries:
        with tempfile.TemporaryDirectory(prefix="sanger16s_") as temp_dir:
            query_fasta = Path(temp_dir) / "queries.fasta"
            write_query_fasta(queries, query_fasta)
            hits = run_blast(query_fasta, db_dir, args.max_target_seqs, args.threads)
        grouped_hits = group_best_hits(hits)
        rows.extend(result_row_for_query(record, grouped_hits.get(record.internal_id, [])) for record in queries)

    write_csv(rows, out_path)
    print(f"Wrote {len(rows)} result row(s) to {out_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch identify Sanger 16S FASTA sequences using local NCBI BLAST+.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser("check", help="Check BLAST+ executables and local database files.")
    check.add_argument("--db-dir", default="db", help="Directory containing the 16S BLAST database.")
    check.set_defaults(func=command_check)

    setup_db = subparsers.add_parser("setup-db", help="Download and decompress the NCBI 16S BLAST database.")
    setup_db.add_argument("--db-dir", default="db", help="Directory where the database will be stored.")
    setup_db.set_defaults(func=command_setup_db)

    run = subparsers.add_parser("run", help="Run batch 16S identification and write one CSV file.")
    run.add_argument("--input", required=True, help="Input FASTA file or directory containing FASTA files.")
    run.add_argument("--db-dir", default="db", help="Directory containing the 16S BLAST database.")
    run.add_argument("--out", required=True, help="Output CSV path.")
    run.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1), help="Number of BLAST threads.")
    run.add_argument("--max-target-seqs", type=int, default=10, help="Top BLAST hits retained per query.")
    run.set_defaults(func=command_run)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (OSError, RuntimeError, ValueError) as exc:
        eprint(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
