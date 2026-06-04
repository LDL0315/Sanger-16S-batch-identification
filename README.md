# Sanger 16S batch identification

This tool batch-processes Sanger 16S FASTA files and reports the best local
NCBI 16S BLAST hit for each sequence in one CSV file.

## 1. Install on a Linux server

Use Conda or Mamba. No `sudo` is required.

```bash
conda env create -f environment.yml
conda activate sanger16s
```

If you use Mamba:

```bash
mamba env create -f environment.yml
conda activate sanger16s
```

## 2. Check the environment

```bash
python sanger16s.py check --db-dir db
```

Before downloading the database, this command should find `blastn` and
`update_blastdb.pl`, but it will report that the database is missing.

## 3. Download the local NCBI 16S database

The server must be able to connect to NCBI.

```bash
python sanger16s.py setup-db --db-dir db
```

This downloads and decompresses the NCBI `16S_ribosomal_RNA` BLAST database
into `db/`.

## 4. Run batch identification

Put your FASTA files into one folder, for example:

```text
fasta_files/
  sample_001.fasta
  sample_002.fasta
  sample_003.fa
```

Then run:

```bash
python sanger16s.py run --input fasta_files --db-dir db --out results.csv
```

The script reads `.fasta`, `.fa`, `.fas`, and `.fna` files recursively. The
file name is used as the sample name. If one FASTA file contains multiple
records, each record is reported as one CSV row.

## CSV columns

The output CSV contains:

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

## Confidence labels

The best BLAST hit is selected by:

1. highest `percent_identity`
2. highest `query_coverage`
3. lowest `evalue`
4. highest `bitscore`

Confidence labels are:

- `species_high`: identity >= 99.0 and coverage >= 95.
- `genus_or_close_species`: identity >= 97.0 and coverage >= 90.
- `low_confidence`: below those thresholds, but the best hit is still reported.
- `no_result`: empty sequence, empty file, or no BLAST hit.

If several species-level hits are extremely close, the `note` column reports
`multiple close hits; review manually`.

## Useful options

Use more or fewer CPU threads:

```bash
python sanger16s.py run --input fasta_files --db-dir db --out results.csv --threads 8
```

Keep more BLAST candidate hits internally when deciding whether close hits
exist:

```bash
python sanger16s.py run --input fasta_files --db-dir db --out results.csv --max-target-seqs 20
```

## Notes

16S identity is not always enough to distinguish very closely related bacterial
species. This tool reports the highest identity hit plus a confidence label, so
low-confidence results are not forced into a definitive species call.
