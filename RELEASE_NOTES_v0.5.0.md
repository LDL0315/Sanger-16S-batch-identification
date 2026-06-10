# Release Notes: v0.5.0

Sanger16S `v0.5.0` is the first public repository-ready release with both CLI and LAN web UI workflows.

## Highlights

- Added a browser-based upload and job board for LAN deployment
- Added FASTQ support through the consensus pipeline
- Unified `FASTA -> BLAST` and `FASTQ -> consensus -> BLAST` into one task system
- Added deployment examples for `systemd`, `nginx`, and env-based service startup
- Added English and Chinese documentation for public use

## Included in this release

### CLI

- `check` for dependency and database validation
- `setup-db` for downloading the NCBI 16S BLAST database
- `run` for direct FASTA batch identification
- `serve` for starting the LAN web interface

### Web UI

- Multi-file upload
- Job queueing and status polling
- FASTA and FASTQ workflow support
- Result preview and CSV download
- Software version and database version display

### FASTQ consensus workflow

- Consensus generation from uploaded FASTQ files
- Cluster abundance filtering
- Orientation normalization
- SNP parsing and consensus report generation

## Deployment notes

- Recommended platform: Linux
- Recommended environment manager: Conda or Mamba
- Web mode requires `Flask`
- FASTQ mode requires `bash`, `gzip`, and the helper scripts under `scripts/`
- See `deploy/` for example service and reverse-proxy configuration

## Upgrade notes

- Default `max-target-seqs` is now `1`
- Public repository layout now includes `README.md`, `LICENSE`, `CHANGELOG.md`, and deployment examples
- Ignore rules now exclude generated work directories and cache files

## Suggested GitHub release title

```text
v0.5.0 - First public web-enabled release
```

## Suggested GitHub release summary

```text
Sanger16S v0.5.0 is the first public release with both CLI and LAN web UI support. It adds FASTQ consensus-to-BLAST workflows, deployment examples, repository documentation, and a production-oriented upload dashboard for batch 16S identification.
```
