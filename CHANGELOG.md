# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2026-06-05

### Added

- Added a LAN web UI for batch sequence submission, job tracking, result preview, and CSV download.
- Added FASTQ support through `run_16s_consensus_pipeline.sh`, enabling `FASTQ -> consensus -> BLAST` workflows.
- Added helper scripts under `scripts/` for consensus orientation, cluster selection, SNP parsing, and consensus reporting.
- Added health checks for BLAST, database readiness, consensus pipeline availability, and helper file completeness.
- Added deployment examples for `systemd`, `nginx`, and environment-based startup configuration.
- Added English and Chinese documentation for CLI usage and web deployment.
- Added MIT license and repository publishing metadata files.

### Changed

- Refactored the CLI so the web workflow reuses the same 16S analysis logic as the direct FASTA workflow.
- Updated default `max-target-seqs` to `1` for both CLI and web UI.
- Expanded `environment.yml` to cover FASTQ consensus pipeline dependencies.
- Improved the web frontend with a dashboard layout, version display, footer attribution, and stable polling behavior.

### Fixed

- Reduced task-board flicker by only re-rendering job cards when job data actually changes.
- Ensured result preview tables always render all CSV columns in a fixed field order.
- Improved upload validation by rejecting mixed FASTA and FASTQ inputs in a single job.

## [0.3.0] - 2026-06-05

### Added

- Added the first web-accessible interface on top of the local BLAST workflow.
- Added upload, background job execution, CSV download, and job status polling.
- Added initial deployment notes for LAN use.

## [0.1.0] - 2026-05-28

### Added

- Initial CLI release for batch 16S identification from FASTA files using a local NCBI BLAST database.
- Added database setup, environment checks, FASTA parsing, sequence sanitization, BLAST execution, and CSV reporting.
