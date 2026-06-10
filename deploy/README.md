# Deployment Examples

This directory contains example files for deploying the Sanger16S web service on a Linux server.

## Files

- `sanger16s-web.env.example`
  - Environment variables for service startup
- `systemd/sanger16s-web.service.example`
  - Example `systemd` unit file
- `nginx/sanger16s-web.conf.example`
  - Example `nginx` reverse-proxy configuration

## Typical Usage

1. Copy `sanger16s-web.env.example` to a real environment file.
2. Adjust paths for:
   - the project directory
   - the local BLAST database directory
   - the consensus pipeline script
   - the web working directory
3. Copy the `systemd` unit example and update it for your server.
4. Optionally add `nginx` as a reverse proxy.

## Notes

- The example files are templates, not drop-in universal defaults.
- Always verify the `--db-dir` path points to the directory that contains files such as:
  - `16S_ribosomal_RNA.nsq`
  - `16S_ribosomal_RNA.nin`
  - `16S_ribosomal_RNA.nhr`
- FASTQ mode also requires `run_16s_consensus_pipeline.sh` and the helper scripts under `scripts/`.
