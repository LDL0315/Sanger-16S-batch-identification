#!/usr/bin/env bash

set -euo pipefail
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLUSTER_HELPER="${SCRIPT_DIR}/scripts/select_consensus_clusters.py"
PILEUP_HELPER="${SCRIPT_DIR}/scripts/parse_mpileup_snp.py"
ORIENT_HELPER="${SCRIPT_DIR}/scripts/canonicalize_fasta_orientation.py"
REPORT_HELPER="${SCRIPT_DIR}/scripts/build_consensus_report.py"

MIN_LEN=1300
MAX_LEN=1700
IDENTITY=0.99
MIN_READS=20
MIN_FRAC=0.01
MAPQ=20
INPUT_DIR=""
OUTPUT_DIR=""
FORWARD_PRIMER="AGAGTTTGATCCTGGCTCAG"
REVERSE_PRIMER="ACGGCTACCTTGTTACGACTT"

usage() {
    cat <<'EOF'
Usage:
  run_16s_consensus_pipeline.sh -i <input_dir> -o <output_dir> [options]

Required:
  -i, --input-dir     Directory containing per-sample FASTQ/FASTQ.GZ files
  -o, --output-dir    Output directory

Options:
      --min-len       Minimum read length to keep (default: 1300)
      --max-len       Maximum read length to keep (default: 1700)
      --id            Clustering identity threshold (default: 0.99)
      --min-reads     Minimum reads per cluster to keep (default: 20)
      --min-frac      Minimum cluster fraction to keep (default: 0.01)
      --mapq          Minimum MAPQ for SNP pileup alignments (default: 20)
      --forward-primer Forward primer sequence used for 5'-to-3' orientation
      --reverse-primer Reverse primer sequence used for 5'-to-3' orientation
  -h, --help          Show this help message
EOF
}

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

fail() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

check_dependency() {
    local tool="$1"
    command -v "${tool}" >/dev/null 2>&1 || fail "Missing required dependency: ${tool}"
}

strip_fastq_extension() {
    local name="$1"
    name="${name%.fastq.gz}"
    name="${name%.fq.gz}"
    name="${name%.fastq}"
    name="${name%.fq}"
    printf '%s\n' "${name}"
}

append_summary() {
    local summary_file="$1"
    local sample="$2"
    local status="$3"
    local filtered_reads="$4"
    local passed_consensus="$5"
    local reason="$6"
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${sample}" "${status}" "${filtered_reads}" "${passed_consensus}" "${reason}" >> "${summary_file}"
}

append_collection_csv() {
    local csv_file="$1"
    local sample="$2"
    local consensus_id="$3"
    local fasta_name="$4"
    local read_count="$5"
    local read_fraction="$6"
    local length="$7"
    printf '%s,%s,%s,%s,%s,%s,%s\n' \
        "${sample}" "${consensus_id}" "${fasta_name}" "${read_count}" "${read_fraction}" "${length}" "${OUTPUT_DIR}/all_sample_fastas/${fasta_name}" >> "${csv_file}"
}

write_empty_reference_counts() {
    local output_path="$1"
    cat > "${output_path}" <<'EOF'
sample	consensus_id	mapped_reads	mapped_fraction
EOF
}

reset_sample_outputs() {
    local sample_dir="$1"
    if [[ -d "${sample_dir}" ]]; then
        rm -rf "${sample_dir}"
    fi
}

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -i|--input-dir)
                INPUT_DIR="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --min-len)
                MIN_LEN="$2"
                shift 2
                ;;
            --max-len)
                MAX_LEN="$2"
                shift 2
                ;;
            --id)
                IDENTITY="$2"
                shift 2
                ;;
            --min-reads)
                MIN_READS="$2"
                shift 2
                ;;
            --min-frac)
                MIN_FRAC="$2"
                shift 2
                ;;
            --mapq)
                MAPQ="$2"
                shift 2
                ;;
            --forward-primer)
                FORWARD_PRIMER="$2"
                shift 2
                ;;
            --reverse-primer)
                REVERSE_PRIMER="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                fail "Unknown argument: $1"
                ;;
        esac
    done

    [[ -n "${INPUT_DIR}" ]] || fail "Input directory is required"
    [[ -n "${OUTPUT_DIR}" ]] || fail "Output directory is required"
    [[ -d "${INPUT_DIR}" ]] || fail "Input directory not found: ${INPUT_DIR}"
    [[ -f "${CLUSTER_HELPER}" ]] || fail "Cluster helper not found: ${CLUSTER_HELPER}"
    [[ -f "${PILEUP_HELPER}" ]] || fail "Pileup helper not found: ${PILEUP_HELPER}"
    [[ -f "${ORIENT_HELPER}" ]] || fail "Orientation helper not found: ${ORIENT_HELPER}"
    [[ -f "${REPORT_HELPER}" ]] || fail "Report helper not found: ${REPORT_HELPER}"
}

check_dependencies() {
    check_dependency gzip
    check_dependency seqkit
    check_dependency vsearch
    check_dependency minimap2
    check_dependency samtools
    check_dependency python3
}

main() {
    parse_args "$@"
    check_dependencies

    mkdir -p "${OUTPUT_DIR}"

    local summary_file="${OUTPUT_DIR}/pipeline_summary.tsv"
    local collection_dir="${OUTPUT_DIR}/all_sample_fastas"
    local collection_csv="${collection_dir}/consensus_summary.csv"
    local report_csv="${collection_dir}/consensus_report.csv"
    cat > "${summary_file}" <<'EOF'
sample	status	filtered_reads	passed_consensus	reason
EOF
    mkdir -p "${collection_dir}"
    rm -f "${collection_dir}"/*.fasta "${collection_dir}"/*.csv
    cat > "${collection_csv}" <<'EOF'
sample,consensus_id,fasta_name,read_count,read_fraction,length,fasta_path
EOF

    local input_files=(
        "${INPUT_DIR}"/*.fastq
        "${INPUT_DIR}"/*.fastq.gz
        "${INPUT_DIR}"/*.fq
        "${INPUT_DIR}"/*.fq.gz
    )

    if [[ ${#input_files[@]} -eq 0 ]]; then
        fail "No FASTQ or FASTQ.GZ files found in ${INPUT_DIR}"
    fi

    local total_samples=0
    local success_samples=0

    for input_fastq in "${input_files[@]}"; do
        [[ -f "${input_fastq}" ]] || continue
        total_samples=$((total_samples + 1))

        local sample_name
        sample_name="$(strip_fastq_extension "$(basename "${input_fastq}")")"
        log "Processing sample: ${sample_name}"

        local sample_dir="${OUTPUT_DIR}/${sample_name}"
        reset_sample_outputs "${sample_dir}"
        local filter_dir="${sample_dir}/01_filter"
        local cluster_dir="${sample_dir}/02_cluster"
        local consensus_dir="${sample_dir}/03_consensus"
        local align_dir="${sample_dir}/04_align"
        local snp_dir="${sample_dir}/05_snp"
        mkdir -p "${filter_dir}" "${cluster_dir}" "${consensus_dir}" "${align_dir}" "${snp_dir}"

        local filtered_fastq="${filter_dir}/${sample_name}.filtered.fastq.gz"
        local filtered_fasta="${filter_dir}/${sample_name}.filtered.fasta"
        local oriented_fasta="${filter_dir}/${sample_name}.filtered.oriented.fasta"
        local filtered_stats="${filter_dir}/${sample_name}.filtered.stats.tsv"
        local derep_fasta="${cluster_dir}/${sample_name}.derep.fasta"
        local all_consensus_fasta="${cluster_dir}/all_consensus.fasta"
        local uc_file="${cluster_dir}/clusters.uc"
        local passed_consensus_fasta="${consensus_dir}/consensus_passed.fasta"
        local abundance_tsv="${consensus_dir}/consensus_abundance.tsv"
        local split_dir="${consensus_dir}"
        local reference_counts_tsv="${align_dir}/reference_read_counts.tsv"
        local pileup_file="${snp_dir}/${sample_name}.mpileup.tsv"
        local filtered_reads=0
        local passed_consensus=0

        seqkit seq -m "${MIN_LEN}" -M "${MAX_LEN}" "${input_fastq}" | gzip -c > "${filtered_fastq}"
        seqkit stats -T "${filtered_fastq}" > "${filtered_stats}"
        filtered_reads="$(awk 'NR==2 {print $4+0}' "${filtered_stats}")"

        if [[ "${filtered_reads}" -eq 0 ]]; then
            log "Sample ${sample_name}: no reads remained after length filtering"
            write_empty_reference_counts "${reference_counts_tsv}"
            append_summary "${summary_file}" "${sample_name}" "failed" "0" "0" "no_reads_after_length_filter"
            continue
        fi

        seqkit fq2fa "${filtered_fastq}" -o "${filtered_fasta}"
        python3 "${ORIENT_HELPER}" \
            --input "${filtered_fasta}" \
            --output "${oriented_fasta}" \
            --forward-primer "${FORWARD_PRIMER}" \
            --reverse-primer "${REVERSE_PRIMER}"
        vsearch \
            --derep_fulllength "${oriented_fasta}" \
            --sizeout \
            --fasta_width 0 \
            --output "${derep_fasta}" >/dev/null 2>&1

        vsearch \
            --cluster_size "${derep_fasta}" \
            --id "${IDENTITY}" \
            --sizein \
            --sizeout \
            --fasta_width 0 \
            --centroids "${all_consensus_fasta}" \
            --uc "${uc_file}" >/dev/null 2>&1

        python3 "${CLUSTER_HELPER}" \
            --sample "${sample_name}" \
            --centroids "${all_consensus_fasta}" \
            --uc "${uc_file}" \
            --total-reads "${filtered_reads}" \
            --min-reads "${MIN_READS}" \
            --min-frac "${MIN_FRAC}" \
            --passed-fasta "${passed_consensus_fasta}" \
            --split-dir "${split_dir}" \
            --abundance-tsv "${abundance_tsv}"

        passed_consensus="$(awk 'NR>1 {count++} END {print count+0}' "${abundance_tsv}")"

        if [[ "${passed_consensus}" -eq 0 ]]; then
            log "Sample ${sample_name}: no clusters passed abundance thresholds"
            write_empty_reference_counts "${reference_counts_tsv}"
            append_summary "${summary_file}" "${sample_name}" "failed" "${filtered_reads}" "0" "no_clusters_passed_threshold"
            continue
        fi

        while IFS=$'\t' read -r sample_value consensus_id read_count read_fraction length_value; do
            [[ "${sample_value}" == "sample" ]] && continue
            [[ -n "${consensus_id}" ]] || continue

            local source_fasta="${consensus_dir}/${consensus_id}.fasta"
            local target_name
            if [[ "${passed_consensus}" -eq 1 ]]; then
                target_name="${sample_name}.fasta"
            else
                target_name="${sample_name}__${consensus_id}.fasta"
            fi

            cp "${source_fasta}" "${collection_dir}/${target_name}"
            append_collection_csv \
                "${collection_csv}" \
                "${sample_value}" \
                "${consensus_id}" \
                "${target_name}" \
                "${read_count}" \
                "${read_fraction}" \
                "${length_value}"
        done < "${abundance_tsv}"

        rm -f "${passed_consensus_fasta}.fai"
        samtools faidx "${passed_consensus_fasta}"

        local sorted_bam="${align_dir}/${sample_name}.primary.q${MAPQ}.bam"
        minimap2 -ax map-ont "${passed_consensus_fasta}" "${filtered_fastq}" \
            | samtools view -b -F 2308 -q "${MAPQ}" - \
            | samtools sort -o "${sorted_bam}" -
        samtools index "${sorted_bam}"

        {
            printf 'sample\tconsensus_id\tmapped_reads\tmapped_fraction\n'
            samtools idxstats "${sorted_bam}" \
                | awk -v sample="${sample_name}" -v total="${filtered_reads}" 'BEGIN {OFS="\t"} $1 != "*" {frac=(total>0 ? $3/total : 0); printf "%s\t%s\t%d\t%.6f\n", sample, $1, $3, frac}'
        } > "${reference_counts_tsv}"

        samtools mpileup -aa -A -B -Q 0 -f "${passed_consensus_fasta}" "${sorted_bam}" > "${pileup_file}"

        python3 "${PILEUP_HELPER}" \
            --sample "${sample_name}" \
            --pileup "${pileup_file}" \
            --output-dir "${snp_dir}"

        append_summary "${summary_file}" "${sample_name}" "success" "${filtered_reads}" "${passed_consensus}" "ok"
        success_samples=$((success_samples + 1))
    done

    python3 "${REPORT_HELPER}" \
        --output-root "${OUTPUT_DIR}" \
        --report-csv "${report_csv}"

    log "Completed ${success_samples}/${total_samples} samples successfully"
    log "Summary: ${summary_file}"
    log "Consensus report: ${report_csv}"
}

main "$@"
