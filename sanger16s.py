#!/usr/bin/env python3
"""Batch identify Sanger 16S FASTA sequences with a local NCBI BLAST database."""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


DB_NAME = "16S_ribosomal_RNA"
APP_VERSION = "0.5.0"
FASTA_SUFFIXES = (".fasta", ".fa", ".fas", ".fna")
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
CONSENSUS_HELPER_FILES = (
    "scripts/select_consensus_clusters.py",
    "scripts/parse_mpileup_snp.py",
    "scripts/canonicalize_fasta_orientation.py",
    "scripts/build_consensus_report.py",
)
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

WEB_PREVIEW_ROWS = 20


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


def path_has_suffix(path: Path | str, suffixes: Sequence[str]) -> bool:
    name = path.name if isinstance(path, Path) else str(path)
    lower_name = name.lower()
    return any(lower_name.endswith(suffix) for suffix in suffixes)


def strip_known_suffix(name: str, suffixes: Sequence[str]) -> str:
    lower_name = name.lower()
    for suffix in sorted(suffixes, key=len, reverse=True):
        if lower_name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def resolve_db_version(db_dir: Path, configured_version: str) -> str:
    version = configured_version.strip()
    if version:
        return version

    prefix = resolve_db_prefix(db_dir)
    candidates = [prefix.with_suffix(ext) for ext in (".nhr", ".nin", ".nsq", ".ndb", ".ntf", ".nto")]
    existing = [path for path in candidates if path.exists()]
    if not existing:
        return "unknown"

    latest_mtime = max(path.stat().st_mtime for path in existing)
    return datetime.fromtimestamp(latest_mtime).strftime("%Y-%m-%d")


def discover_fasta_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        if not path_has_suffix(input_path, FASTA_SUFFIXES):
            raise ValueError(f"Input file does not look like FASTA: {input_path}")
        return [input_path]
    if not input_path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {input_path}")
    return sorted(
        path
        for path in input_path.rglob("*")
        if path.is_file() and path_has_suffix(path, FASTA_SUFFIXES)
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


def analyze_to_csv(
    input_path: Path,
    db_dir: Path,
    out_path: Path,
    *,
    max_target_seqs: int,
    threads: int,
) -> int:
    require_executable("blastn")

    if not database_exists(db_dir):
        raise RuntimeError(f"BLAST database not found at {resolve_db_prefix(db_dir)}. Run setup-db first.")

    queries, skipped_rows = load_queries(input_path)
    rows: list[dict[str, str]] = list(skipped_rows)

    if queries:
        with tempfile.TemporaryDirectory(prefix="sanger16s_") as temp_dir:
            query_fasta = Path(temp_dir) / "queries.fasta"
            write_query_fasta(queries, query_fasta)
            hits = run_blast(query_fasta, db_dir, max_target_seqs, threads)
        grouped_hits = group_best_hits(hits)
        rows.extend(result_row_for_query(record, grouped_hits.get(record.internal_id, [])) for record in queries)

    write_csv(rows, out_path)
    return len(rows)


def detect_input_kind(file_names: Sequence[str]) -> str:
    detected_kinds: list[str] = []
    for name in file_names:
        if path_has_suffix(name, FASTA_SUFFIXES):
            detected_kinds.append("fasta")
        elif path_has_suffix(name, FASTQ_SUFFIXES):
            detected_kinds.append("fastq")
        else:
            raise ValueError(f"Unsupported file type: {name}")

    if not detected_kinds:
        raise ValueError("No supported sequence files were provided.")
    kinds = set(detected_kinds)
    if len(kinds) > 1:
        raise ValueError("Please upload either FASTA files or FASTQ files in one job, not both.")
    return detected_kinds[0]


def run_consensus_pipeline(
    *,
    pipeline_script: Path,
    input_dir: Path,
    output_dir: Path,
) -> None:
    bash_path = require_executable("bash")
    if not pipeline_script.exists():
        raise RuntimeError(f"Consensus pipeline script not found: {pipeline_script}")

    cmd = [
        bash_path,
        str(pipeline_script),
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
    ]
    result = run_command(cmd)
    if result.returncode != 0:
        stderr = result.stderr.strip()
        stdout = result.stdout.strip()
        details = stderr or stdout or "unknown error"
        raise RuntimeError(f"Consensus pipeline failed: {details}")


def missing_consensus_pipeline_files(pipeline_script: Path) -> list[str]:
    missing: list[str] = []
    if not pipeline_script.exists():
        missing.append(str(pipeline_script))
        return missing

    script_dir = pipeline_script.parent
    for relative_path in CONSENSUS_HELPER_FILES:
        candidate = script_dir / relative_path
        if not candidate.exists():
            missing.append(str(candidate))
    return missing


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def clamp_positive(value: Any, default: int, *, maximum: int | None = None) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError):
        return default
    if number < 1:
        return default
    if maximum is not None:
        return min(number, maximum)
    return number


def summarize_csv(out_path: Path, preview_limit: int = WEB_PREVIEW_ROWS) -> tuple[int, dict[str, int], list[dict[str, str]]]:
    total_rows = 0
    confidence_counts: dict[str, int] = {}
    preview_rows: list[dict[str, str]] = []

    with out_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            total_rows += 1
            confidence = row.get("confidence", "") or "unknown"
            confidence_counts[confidence] = confidence_counts.get(confidence, 0) + 1
            if len(preview_rows) < preview_limit:
                preview_rows.append(row)

    return total_rows, confidence_counts, preview_rows


def read_status_file(status_path: Path) -> dict[str, Any] | None:
    if not status_path.exists():
        return None
    with status_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_status_file(status_path: Path, payload: dict[str, Any]) -> None:
    status_path.parent.mkdir(parents=True, exist_ok=True)
    with status_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def update_job_status(status_path: Path, lock: threading.Lock, **changes: Any) -> dict[str, Any]:
    with lock:
        payload = read_status_file(status_path) or {}
        payload.update(changes)
        payload["updated_at"] = utc_now_iso()
        write_status_file(status_path, payload)
        return payload


def job_payload_for_api(job: dict[str, Any], *, has_result: bool) -> dict[str, Any]:
    payload = dict(job)
    payload["has_result"] = has_result
    return payload


def collect_job_payloads(jobs_root: Path) -> list[dict[str, Any]]:
    jobs: list[dict[str, Any]] = []
    for status_path in jobs_root.glob("*/status.json"):
        payload = read_status_file(status_path)
        if not payload:
            continue
        payload["has_result"] = (status_path.parent / "result.csv").exists()
        jobs.append(payload)
    jobs.sort(key=lambda item: item.get("created_at", ""), reverse=True)
    return jobs


def create_job_id() -> str:
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    return f"job-{timestamp}-{uuid.uuid4().hex[:8]}"


def create_web_app(
    *,
    db_dir: Path,
    work_dir: Path,
    db_version: str,
    pipeline_script: Path,
    default_threads: int,
    default_max_target_seqs: int,
    max_upload_mb: int,
    queue_workers: int,
):
    try:
        from flask import Flask, jsonify, render_template, request, send_file, url_for
        from werkzeug.utils import secure_filename
    except ImportError as exc:
        raise RuntimeError(
            "Web UI requires Flask. Install it first, for example with: conda install flask"
        ) from exc

    base_dir = Path(__file__).resolve().parent
    template_dir = base_dir / "web" / "templates"
    static_dir = base_dir / "web" / "static"
    jobs_root = work_dir / "jobs"
    jobs_root.mkdir(parents=True, exist_ok=True)

    status_lock = threading.Lock()
    executor = ThreadPoolExecutor(max_workers=max(1, queue_workers))

    def recover_interrupted_jobs() -> None:
        for status_path in jobs_root.glob("*/status.json"):
            payload = read_status_file(status_path)
            if payload and payload.get("status") in {"queued", "running"}:
                update_job_status(
                    status_path,
                    status_lock,
                    status="failed",
                    error="Server restarted while this job was still running.",
                    completed_at=utc_now_iso(),
                )

    def save_upload(upload_dir: Path, storage, seen_names: set[str], default_suffix: str) -> str:
        original_name = storage.filename or f"sequence{default_suffix}"
        filename = secure_filename(original_name) or f"sequence{default_suffix}"
        stem = strip_known_suffix(filename, (*FASTA_SUFFIXES, *FASTQ_SUFFIXES)) or "sequence"
        suffixes = Path(filename).suffixes
        suffix = "".join(suffixes) or default_suffix
        candidate = filename
        counter = 1
        while candidate in seen_names:
            candidate = f"{stem}_{counter}{suffix}"
            counter += 1
        seen_names.add(candidate)
        storage.save(upload_dir / candidate)
        return candidate

    def run_web_job(
        *,
        status_path: Path,
        input_dir: Path,
        result_path: Path,
        input_kind: str,
        threads: int,
        max_target_seqs: int,
    ) -> None:
        update_job_status(
            status_path,
            status_lock,
            status="running",
            stage="consensus" if input_kind == "fastq" else "blast",
            started_at=utc_now_iso(),
            error="",
        )
        try:
            analysis_input_dir = input_dir
            if input_kind == "fastq":
                pipeline_output_dir = status_path.parent / "pipeline-output"
                run_consensus_pipeline(
                    pipeline_script=pipeline_script,
                    input_dir=input_dir,
                    output_dir=pipeline_output_dir,
                )
                analysis_input_dir = pipeline_output_dir / "all_sample_fastas"
                if not analysis_input_dir.exists():
                    raise RuntimeError(
                        f"Consensus pipeline finished but no FASTA output directory was created: {analysis_input_dir}"
                    )
                update_job_status(
                    status_path,
                    status_lock,
                    stage="blast",
                    pipeline_output_dir=str(pipeline_output_dir),
                    generated_fasta_dir=str(analysis_input_dir),
                )

            row_count = analyze_to_csv(
                analysis_input_dir,
                db_dir,
                result_path,
                max_target_seqs=max_target_seqs,
                threads=threads,
            )
            total_rows, confidence_counts, preview_rows = summarize_csv(result_path)
            update_job_status(
                status_path,
                status_lock,
                status="completed",
                row_count=row_count,
                total_rows=total_rows,
                confidence_counts=confidence_counts,
                preview_rows=preview_rows,
                stage="completed",
                completed_at=utc_now_iso(),
                error="",
            )
        except Exception as exc:
            update_job_status(
                status_path,
                status_lock,
                status="failed",
                stage="failed",
                error=str(exc),
                completed_at=utc_now_iso(),
            )

    recover_interrupted_jobs()

    app = Flask(
        __name__,
        template_folder=str(template_dir),
        static_folder=str(static_dir),
    )
    app.config["MAX_CONTENT_LENGTH"] = max_upload_mb * 1024 * 1024

    @app.get("/")
    def index():
        return render_template(
            "index.html",
            config={
                "defaultThreads": default_threads,
                "defaultMaxTargetSeqs": default_max_target_seqs,
                "maxUploadMb": max_upload_mb,
                "dbPrefix": str(resolve_db_prefix(db_dir)),
                "csvFields": CSV_FIELDS,
                "appVersion": APP_VERSION,
                "dbVersion": db_version,
                "supportedInputText": "FASTA / FASTQ(.gz)",
            },
        )

    @app.get("/api/health")
    def api_health():
        blastn_path = shutil.which("blastn")
        bash_path = shutil.which("bash")
        missing_pipeline_files = missing_consensus_pipeline_files(pipeline_script)
        return jsonify(
            {
                "ok": bool(blastn_path) and database_exists(db_dir),
                "blastn": blastn_path or "",
                "bash": bash_path or "",
                "database_ready": database_exists(db_dir),
                "database_prefix": str(resolve_db_prefix(db_dir)),
                "consensus_pipeline_ready": not missing_pipeline_files and bool(bash_path),
                "consensus_pipeline_script": str(pipeline_script),
                "consensus_pipeline_missing_files": missing_pipeline_files,
                "jobs_root": str(jobs_root),
            }
        )

    @app.get("/api/jobs")
    def api_jobs():
        jobs = collect_job_payloads(jobs_root)
        return jsonify(
            [
                {
                    **job_payload_for_api(job, has_result=job.get("has_result", False)),
                    "download_url": url_for("download_result", job_id=job["job_id"]) if job.get("has_result") else "",
                }
                for job in jobs
            ]
        )

    @app.get("/api/jobs/<job_id>")
    def api_job(job_id: str):
        job_dir = jobs_root / job_id
        payload = read_status_file(job_dir / "status.json")
        if not payload:
            return jsonify({"error": "Job not found"}), 404
        has_result = (job_dir / "result.csv").exists()
        return jsonify(
            {
                **job_payload_for_api(payload, has_result=has_result),
                "download_url": url_for("download_result", job_id=job_id) if has_result else "",
            }
        )

    @app.post("/api/jobs")
    def api_create_job():
        if not database_exists(db_dir):
            return jsonify({"error": f"BLAST database not found at {resolve_db_prefix(db_dir)}"}), 400
        if not shutil.which("blastn"):
            return jsonify({"error": "blastn not found in PATH"}), 400

        uploads = request.files.getlist("files")
        if not uploads:
            return jsonify({"error": "Please upload at least one FASTA or FASTQ file."}), 400

        cpu_max = max(1, os.cpu_count() or 1)
        threads = clamp_positive(request.form.get("threads"), default_threads, maximum=cpu_max)
        max_target_seqs = clamp_positive(request.form.get("max_target_seqs"), default_max_target_seqs)

        job_id = create_job_id()
        job_dir = jobs_root / job_id
        input_dir = job_dir / "input"
        result_path = job_dir / "result.csv"
        status_path = job_dir / "status.json"
        input_dir.mkdir(parents=True, exist_ok=True)

        accepted_files: list[str] = []
        seen_names: set[str] = set()
        original_names = [storage.filename for storage in uploads if storage.filename]
        try:
            input_kind = detect_input_kind(original_names)
        except ValueError as exc:
            return jsonify({"error": str(exc)}), 400

        missing_pipeline_files = missing_consensus_pipeline_files(pipeline_script)
        if input_kind == "fastq" and (missing_pipeline_files or not shutil.which("bash")):
            details = ", ".join(missing_pipeline_files) if missing_pipeline_files else "bash not found in PATH"
            return jsonify({"error": f"FASTQ mode is unavailable because the consensus pipeline is incomplete: {details}"}), 400

        for storage in uploads:
            if not storage.filename:
                continue
            if not path_has_suffix(storage.filename, FASTA_SUFFIXES if input_kind == "fasta" else FASTQ_SUFFIXES):
                return jsonify({"error": f"Unsupported file type: {storage.filename}"}), 400
            default_suffix = ".fasta" if input_kind == "fasta" else ".fastq"
            accepted_files.append(save_upload(input_dir, storage, seen_names, default_suffix))

        if not accepted_files:
            return jsonify({"error": "No valid sequence files were uploaded."}), 400

        write_status_file(
            status_path,
            {
                "job_id": job_id,
                "status": "queued",
                "stage": "queued",
                "input_kind": input_kind,
                "created_at": utc_now_iso(),
                "updated_at": utc_now_iso(),
                "files": accepted_files,
                "file_count": len(accepted_files),
                "threads": threads,
                "max_target_seqs": max_target_seqs,
                "row_count": 0,
                "total_rows": 0,
                "confidence_counts": {},
                "preview_rows": [],
                "error": "",
            },
        )

        executor.submit(
            run_web_job,
            status_path=status_path,
            input_dir=input_dir,
            result_path=result_path,
            input_kind=input_kind,
            threads=threads,
            max_target_seqs=max_target_seqs,
        )

        payload = read_status_file(status_path) or {}
        return (
            jsonify(
                {
                    **job_payload_for_api(payload, has_result=False),
                    "download_url": "",
                }
            ),
            202,
        )

    @app.get("/jobs/<job_id>/result.csv")
    def download_result(job_id: str):
        result_path = jobs_root / job_id / "result.csv"
        if not result_path.exists():
            return jsonify({"error": "Result not found"}), 404
        return send_file(result_path, as_attachment=True, download_name=f"{job_id}.csv")

    return app


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
        print("Run: python sanger16s.py setup-db --db-dir /path/to/local-16s-db-dir")
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
    input_path = Path(args.input)
    db_dir = Path(args.db_dir)
    out_path = Path(args.out)
    row_count = analyze_to_csv(
        input_path,
        db_dir,
        out_path,
        max_target_seqs=args.max_target_seqs,
        threads=args.threads,
    )
    print(f"Wrote {row_count} result row(s) to {out_path}")
    return 0


def command_serve(args: argparse.Namespace) -> int:
    work_dir = Path(args.work_dir)
    db_dir = Path(args.db_dir)
    pipeline_script = Path(args.consensus_pipeline_script)
    app = create_web_app(
        db_dir=db_dir,
        work_dir=work_dir,
        db_version=resolve_db_version(db_dir, args.db_version),
        pipeline_script=pipeline_script,
        default_threads=args.threads,
        default_max_target_seqs=args.max_target_seqs,
        max_upload_mb=args.max_upload_mb,
        queue_workers=args.queue_workers,
    )
    app.run(host=args.host, port=args.port, debug=False)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch identify Sanger 16S FASTA sequences using local NCBI BLAST+.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    check = subparsers.add_parser("check", help="Check BLAST+ executables and local database files.")
    check.add_argument(
        "--db-dir",
        default="db",
        help="Directory containing the local 16S BLAST database files, for example /data/blast/16s-db.",
    )
    check.set_defaults(func=command_check)

    setup_db = subparsers.add_parser("setup-db", help="Download and decompress the NCBI 16S BLAST database.")
    setup_db.add_argument(
        "--db-dir",
        default="db",
        help="Directory where the local 16S BLAST database files will be downloaded and stored.",
    )
    setup_db.set_defaults(func=command_setup_db)

    run = subparsers.add_parser("run", help="Run batch 16S identification and write one CSV file.")
    run.add_argument("--input", required=True, help="Input FASTA file or directory containing FASTA files.")
    run.add_argument(
        "--db-dir",
        default="db",
        help="Directory containing the local 16S BLAST database files, for example /data/blast/16s-db.",
    )
    run.add_argument("--out", required=True, help="Output CSV path.")
    run.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1), help="Number of BLAST threads.")
    run.add_argument("--max-target-seqs", type=int, default=1, help="Top BLAST hits retained per query.")
    run.set_defaults(func=command_run)

    serve = subparsers.add_parser("serve", help="Run the LAN web UI for batch FASTA or FASTQ upload and CSV download.")
    serve.add_argument(
        "--db-dir",
        default="db",
        help="Directory containing the local 16S BLAST database files, for example /data/blast/16s-db.",
    )
    serve.add_argument("--db-version", default="", help="Displayed 16S database version label. Defaults to latest database file date.")
    serve.add_argument(
        "--consensus-pipeline-script",
        default=str(Path(__file__).resolve().with_name("run_16s_consensus_pipeline.sh")),
        help="Shell script used to turn uploaded FASTQ files into consensus FASTA files before BLAST analysis.",
    )
    serve.add_argument("--work-dir", default="web-work", help="Directory used for uploaded files, job state, and results.")
    serve.add_argument("--host", default="0.0.0.0", help="Host interface to bind.")
    serve.add_argument("--port", type=int, default=8080, help="Port for the web server.")
    serve.add_argument("--threads", type=int, default=max(1, os.cpu_count() or 1), help="BLAST threads per job.")
    serve.add_argument("--queue-workers", type=int, default=1, help="How many jobs may run at the same time.")
    serve.add_argument("--max-target-seqs", type=int, default=1, help="Top BLAST hits retained per query.")
    serve.add_argument("--max-upload-mb", type=int, default=64, help="Maximum total request size in MB.")
    serve.set_defaults(func=command_serve)

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
