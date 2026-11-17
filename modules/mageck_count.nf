/*
 * MAGECK2 COUNT module — UMI-aware sgRNA counting from FASTQ files
 *
 * Counts sgRNA reads (or unique molecules, when --umi is set) from FASTQ
 * files using the provided library. Supports single-end and paired-end
 * reads: the input list may contain both read pairs, and they are split
 * here from the Illumina file-name convention (..._R1_... / ..._R2_...)
 * into --fastq (R1) and --fastq-2 (R2). In the common "R1 = guide,
 * R2 = UMI" layout the UMI window is located in R2.
 *
 * Emits count table, UMI count table (only when --umi is used),
 * normalized counts, summary, and log.
 */

process MAGECK_COUNT {
    publishDir "${params.output_dir}/mageck_count", mode: 'copy'
    tag "${prefix}"
    label 'mageck_count'

    input:
        path library_file
        path fastqs             // all FASTQ files (R1 + R2 mixed); split by name below
        val sample_labels
        val params_obj

    output:
        path "${prefix}.count.txt",            emit: count_table
        path "${prefix}.umi_count.txt",        optional: true, emit: umi_count
        path "${prefix}.count_normalized.txt",  emit: count_normalized
        path "${prefix}.countsummary.txt",      emit: count_summary
        path "mageck2.count.log",               emit: log

    script:
    prefix = params_obj.output_prefix ?: 'mageck2'

    // Pre-compute UMI arguments in Groovy (params are plain values here)
    def umi_args = ''
    if (params_obj.umi && params_obj.umi != 'none') {
        umi_args = "--umi ${params_obj.umi}"
        if (params_obj.umi_start   != null && params_obj.umi_start   != -1) umi_args += " --umi-start ${params_obj.umi_start}"
        if (params_obj.umi_end     != null && params_obj.umi_end     != -1) umi_args += " --umi-end ${params_obj.umi_end}"
        if (params_obj.umi_start_2 != null && params_obj.umi_start_2 != -1) umi_args += " --umi-start-2 ${params_obj.umi_start_2}"
        if (params_obj.umi_end_2   != null && params_obj.umi_end_2   != -1) umi_args += " --umi-end-2 ${params_obj.umi_end_2}"
    }
    def umi_echo = umi_args ? "echo \"[MAGECK_COUNT] UMI: ${umi_args}\"" : ''

    """
    #! /usr/bin/env bash
    echo "[MAGECK_COUNT] Starting sgRNA counting (mageck2)..."
    echo "[MAGECK_COUNT] Library: ${library_file}"
    echo "[MAGECK_COUNT] Sample labels: ${sample_labels}"
    ${umi_echo}

    # Split the FASTQ list into first/second read pairs from the Illumina
    # file-name convention (..._R1_... / ..._R2_...). Files without a pair
    # marker are treated as single-end (first pair).
    R1_UNSORTED=()
    R2_UNSORTED=()
    for fq in ${fastqs}; do
        case "\${fq}" in
            *_R2_*) R2_UNSORTED+=("\${fq}") ;;
            *)      R1_UNSORTED+=("\${fq}") ;;
        esac
    done

    # Deterministic (sorted) order — required so --fastq / --fastq-2 stay
    # consistent with each other and with --sample-label
    R1_FILES=()
    while IFS= read -r line; do [ -n "\${line}" ] && R1_FILES+=("\${line}"); done < <(printf '%s\\n' "\${R1_UNSORTED[@]}" | sort)
    R2_FILES=()
    while IFS= read -r line; do [ -n "\${line}" ] && R2_FILES+=("\${line}"); done < <(printf '%s\\n' "\${R2_UNSORTED[@]}" | sort)

    if [ \${#R2_FILES[@]} -gt 0 ]; then
        if [ \${#R2_FILES[@]} -ne \${#R1_FILES[@]} ]; then
            echo "[MAGECK_COUNT] ERROR: R1/R2 file count mismatch (\${#R1_FILES[@]} R1 vs \${#R2_FILES[@]} R2)"
            exit 1
        fi
        echo "[MAGECK_COUNT] Paired-end: \${#R1_FILES[@]} R1 + \${#R2_FILES[@]} R2 files"
    else
        echo "[MAGECK_COUNT] Single-end: \${#R1_FILES[@]} files"
    fi

    FASTQ_ARG=""
    for fq in "\${R1_FILES[@]}"; do
        FASTQ_ARG="\${FASTQ_ARG} \${fq}"
    done
    FASTQ2_OPT=""
    if [ \${#R2_FILES[@]} -gt 0 ]; then
        FASTQ2_ARG=""
        for fq in "\${R2_FILES[@]}"; do
            FASTQ2_ARG="\${FASTQ2_ARG} \${fq}"
        done
        FASTQ2_OPT="--fastq-2 \${FASTQ2_ARG}"
    fi

    # Build optional arguments
    OPT_ARGS=""

    if [ -n "${params_obj.sample_labels}" ]; then
        OPT_ARGS="\${OPT_ARGS} --sample-label ${params_obj.sample_labels}"
    fi

    if [ "${params_obj.norm_method}" != "median" ]; then
        OPT_ARGS="\${OPT_ARGS} --norm-method ${params_obj.norm_method}"
    fi

    if [ -n "${params_obj.control_sgrna}" ] && [ "${params_obj.control_sgrna}" != "false" ]; then
        OPT_ARGS="\${OPT_ARGS} --control-sgrna ${params_obj.control_sgrna}"
    fi

    if [ -n "${params_obj.control_gene}" ] && [ "${params_obj.control_gene}" != "false" ]; then
        OPT_ARGS="\${OPT_ARGS} --control-gene ${params_obj.control_gene}"
    fi

    if [ "${params_obj.trim5}" != "AUTO" ]; then
        OPT_ARGS="\${OPT_ARGS} --trim-5 ${params_obj.trim5}"
    fi

    if [ "${params_obj.sgrna_len}" != "20" ]; then
        OPT_ARGS="\${OPT_ARGS} --sgrna-len ${params_obj.sgrna_len}"
    fi

    if [ "${params_obj.count_n}" = "true" ] || [ "${params_obj.count_n}" = true ]; then
        OPT_ARGS="\${OPT_ARGS} --count-n"
    fi

    if [ "${params_obj.reverse_complement}" = "true" ] || [ "${params_obj.reverse_complement}" = true ]; then
        OPT_ARGS="\${OPT_ARGS} --reverse-complement"
    fi

    if [ "${params_obj.html_report}" = "true" ] || [ "${params_obj.html_report}" = true ]; then
        OPT_ARGS="\${OPT_ARGS} --html-report"
    fi

    OPT_ARGS="\${OPT_ARGS} ${umi_args}"

    echo "[MAGECK_COUNT] Running: mageck2 count -l ${library_file} -n ${prefix} \${OPT_ARGS} --fastq \${FASTQ_ARG} \${FASTQ2_OPT}"

    mageck2 count \\
        -l "${library_file}" \\
        -n "${prefix}" \\
        \${OPT_ARGS} \\
        --fastq \${FASTQ_ARG} \\
        \${FASTQ2_OPT}

    # Fixed, unique name so the log never collides with the MLE log
    # when both share the same --output_prefix (KGE_REPORT stages them flat)
    mv "${prefix}.log" "mageck2.count.log"

    echo "[MAGECK_COUNT] Counting complete."
    """
}