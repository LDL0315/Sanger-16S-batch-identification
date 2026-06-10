# Contributing

Thanks for contributing to Sanger16S.

## Scope

This repository contains:

- CLI code for 16S identification
- FASTQ consensus preprocessing helpers
- LAN web UI code
- Deployment examples and user documentation

## Development Setup

Recommended platform: Linux.

Create the environment:

```bash
conda env create -f environment.yml
conda activate sanger16s
```

## Suggested Workflow

1. Create a feature branch.
2. Keep changes focused on one topic.
3. Update documentation when behavior or deployment changes.
4. Run basic checks before opening a pull request.

## Before Opening a Pull Request

Please verify:

- Python files compile successfully
- Command-line help still makes sense
- README or deployment docs are updated when needed
- Generated files such as `web-work/`, `db/`, and `__pycache__/` are not committed

Example local check:

```bash
python -m compileall sanger16s.py scripts
```

## Pull Request Notes

When applicable, include:

- What changed
- Why it changed
- Any deployment impact
- Any database or environment assumptions

## Documentation Style

- Keep command examples copy-pasteable
- Use explicit paths for BLAST database directories
- Avoid ambiguous examples like plain `db` when a real directory path is clearer

## Issues

Bug reports are more useful when they include:

- Input type: FASTA or FASTQ
- Example command or web action
- Error message
- Environment details
- Whether the issue is in CLI mode, web mode, or consensus preprocessing
