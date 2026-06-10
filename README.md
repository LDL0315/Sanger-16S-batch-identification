# Sanger16S

Batch 16S bacterial identification with a local NCBI BLAST database, plus an optional LAN web UI for uploading FASTA or FASTQ files.

## Quick Links

- Chinese CLI guide: [README.zh-CN.md](./README.zh-CN.md)
- Chinese web deployment guide: [WEBAPP.zh-CN.md](./WEBAPP.zh-CN.md)
- Change history: [CHANGELOG.md](./CHANGELOG.md)
- Release draft: [RELEASE_NOTES_v0.5.0.md](./RELEASE_NOTES_v0.5.0.md)
- Contributing guide: [CONTRIBUTING.md](./CONTRIBUTING.md)
- Deployment examples: [deploy/README.md](./deploy/README.md)

## Features

- Batch identification from FASTA files using local `blastn`
- FASTQ-to-consensus preprocessing pipeline for ONT 16S data
- CSV summary output with confidence labels
- LAN web UI for upload, queueing, preview, and CSV download
- Single-file CLI entry point: `sanger16s.py`

## Repository Layout

```text
.
- sanger16s.py
- run_16s_consensus_pipeline.sh
- scripts/
- web/
- deploy/
- environment.yml
- LICENSE
- CHANGELOG.md
- README.md
- README.zh-CN.md
- WEBAPP.zh-CN.md
```

## Requirements

Recommended environment: Conda or Mamba on Linux.

The provided `environment.yml` includes:

- Python 3.10+
- NCBI BLAST+
- Flask
- seqkit
- vsearch
- minimap2
- samtools

System-level tools expected by the FASTQ consensus pipeline:

- `bash`
- `gzip`

## Quick Start

Create the environment:

```bash
conda env create -f environment.yml
conda activate sanger16s
```

Check BLAST and the local database:

```bash
python sanger16s.py check --db-dir /data/blast/16s-db
```

Download the NCBI 16S BLAST database:

```bash
python sanger16s.py setup-db --db-dir /data/blast/16s-db
```

Here `--db-dir` means the local directory that contains the BLAST database files such as
`16S_ribosomal_RNA.nsq` and `16S_ribosomal_RNA.nin`. It is not a fixed literal name like `db`.

Run direct FASTA identification:

```bash
python sanger16s.py run \
  --input ./fasta_files \
  --db-dir /data/blast/16s-db \
  --out ./results.csv
```

## Web UI

Start the LAN web server:

```bash
python sanger16s.py serve \
  --host 0.0.0.0 \
  --port 8080 \
  --db-dir /data/blast/16s-db \
  --db-version 2026-06-05 \
  --consensus-pipeline-script ./run_16s_consensus_pipeline.sh \
  --work-dir ./web-work \
  --threads 8 \
  --queue-workers 1 \
  --max-target-seqs 1 \
  --max-upload-mb 64
```

Supported upload modes:

- FASTA: direct BLAST identification
- FASTQ / FASTQ.GZ: consensus generation first, then BLAST identification

## CLI Commands

```text
check      Check executables and database files
setup-db   Download and decompress the NCBI 16S BLAST database
run        Run batch 16S identification from FASTA input
serve      Start the LAN web UI
```

## Output

The main CSV fields are:

- `sample_file`
- `query_id`
- `query_length`
- `best_accession`
- `best_title`
- `organism_name`
- `percent_identity`
- `query_coverage`
- `alignment_length`
- `mismatches`
- `gaps`
- `evalue`
- `bitscore`
- `confidence`
- `note`

Confidence labels:

- `species_high`
- `genus_or_close_species`
- `low_confidence`
- `no_result`

## Deployment Examples

Example deployment files are included in:

- `deploy/systemd/sanger16s-web.service.example`
- `deploy/nginx/sanger16s-web.conf.example`
- `deploy/sanger16s-web.env.example`

## Notes

- The web server stores uploaded files, job state, and output CSV files under `--work-dir`.
- FASTQ mode requires `run_16s_consensus_pipeline.sh` and all helper scripts under `scripts/`.
- 16S alone does not always resolve very closely related species. Treat low-confidence or ambiguous hits cautiously.

## License

MIT. See [LICENSE](./LICENSE).
