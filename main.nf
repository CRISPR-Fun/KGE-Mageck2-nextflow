/*
 * ==========================================================================
 * kge-mageck2 — UMI-aware CRISPR KO Screen Analysis Pipeline (MAGeCK2)
 * ==========================================================================
 *
 * Driven by the files you provide:
 *   --count_table    → required, always the starting point
 *   --comparisons    → triggers RRA (one run per treatment/control pair)
 *   --design_matrix  → triggers MLE
 *   --comparisons + --design_matrix → runs both RRA and MLE
 *
 * UMI-aware counting (FASTQ/demux modes only): set --umi plus the UMI
 * window to count unique molecules instead of raw reads. For the common
 * "R1 = guide, R2 = UMI" layout:
 *   --umi secondpair --umi_start_2 0 --umi_end_2 10
 *
 * Usage examples:
 *
 *   # RRA only (7 comparisons), from an existing count table:
 *   nextflow run kge-mageck2 -profile conda \
 *     --count_table all_samples.count.txt \
 *     --comparisons comparisons.txt
 *
 *   # UMI count from paired-end FASTQ + RRA + MLE:
 *   nextflow run kge-mageck2 -profile conda \
 *     --fastq_dir ./fastq --library_file library.csv \
 *     --sample-label ctrl,treat1,treat2 \
 *     --umi secondpair --umi_start_2 0 --umi_end_2 10 \
 *     --comparisons comparisons.txt \
 *     --design_matrix design_matrix.txt
 */

nextflow.enable.dsl = 2

// ============================================================================
// INCLUDES
// ============================================================================

include { DEMULTIPLEX        } from './modules/demultiplex'
include { FASTQC_READS       } from './modules/fastqc'
include { MAGECK_COUNT        } from './modules/mageck_count'
include { MAGECK_TEST_RRA     } from './modules/mageck_test'
include { MAGECK_MLE_ANALYSIS } from './modules/mageck_mle'
include { KGE_REPORT           } from './modules/kge_report'
include { MULTIQC_REPORT      } from './modules/multiqc'

// ============================================================================
// WORKFLOW
// ============================================================================

workflow {

    // ------------------------------------------------------------------
    // Normalise booleans: the Nextflow CLI passes boolean overrides as
    // strings, and the string "false" is TRUTHY in Groovy — without this
    // `--fastqc false` would still run FastQC
    // ------------------------------------------------------------------
    def doDemux        = params.demultiplex.toString().toBoolean()
    def doHtmlReport   = params.html_report.toString().toBoolean()
    def skipEnrichment = params.skip_enrichment.toString().toBoolean()
    def runFastqc      = params.fastqc.toString().toBoolean()

    // ------------------------------------------------------------------
    // Validate required inputs
    // ------------------------------------------------------------------
    validate_inputs()

    // MAGeCK2 step logs (count / per-comparison test / MLE) — collected and
    // passed to KGE_REPORT so the HTML report can recover each comparison's
    // samples and embed the logs in the combined log download
    run_logs_ch = Channel.empty()

    // UMI count table from MAGECK_COUNT (only when counting runs with --umi)
    umi_count_collected = []

    // ==================================================================
    // INPUT STAGE: resolve count table from one of three sources
    // ==================================================================

    if (params.count_table) {
        // Mode A: start from existing count table — skip demux + count
        count_table_ch    = Channel.fromPath(params.count_table, checkIfExists: true)
        count_norm_ch     = Channel.empty()
        count_summary_ch  = Channel.empty()

    } else if (doDemux) {
        // Mode B: start from BCL — demux → fastqc → count
        // (the R1/R2 split happens inside MAGECK_COUNT, from the Illumina
        //  file-name convention — bcl2fastq: S1_L001_R1_001.fastq.gz,
        //  bcl-convert: S1_R1_001.fastq.gz)
        demux_result = DEMULTIPLEX(
            file(params.run_folder, checkIfExists: true),
            file(params.sample_sheet, checkIfExists: true)
        )

        fastq_stream = demux_result.out.fastqs
            .flatMap { dir -> dir.listFiles().collect { it } }
            .filter { it.getName() =~ /.*\.(fastq\.gz|fastq|fq\.gz|fq)$/ }

        if (runFastqc) {
            FASTQC_READS(fastq_stream.map { it })
        }

        MAGECK_COUNT(
            file(params.library_file, checkIfExists: true),
            fastq_stream.collect(),
            params.sample_labels,
            params
        )
        count_table_ch    = MAGECK_COUNT.out.count_table
        count_norm_ch     = MAGECK_COUNT.out.count_normalized
        count_summary_ch  = MAGECK_COUNT.out.count_summary
        run_logs_ch       = run_logs_ch.concat(MAGECK_COUNT.out.log)
        umi_count_collected = MAGECK_COUNT.out.umi_count.collect()

    } else if (params.fastq_dir) {
        // Mode C: start from FASTQ files (plain-Groovy listing so the R1/R2
        // patterns can be checked up front; MAGECK_COUNT re-splits by name)
        fastq_dir_f = new File(params.fastq_dir)
        if (!fastq_dir_f.isDirectory()) {
            log.error "ERROR: --fastq_dir '${params.fastq_dir}' is not a directory"
            System.exit(1)
        }
        r1_list = listGlob(fastq_dir_f, params.fastq_pattern)
        r2_list = params.fastq_pattern_2 ? listGlob(fastq_dir_f, params.fastq_pattern_2) : []

        if (r1_list.size() == 0) {
            log.error "ERROR: no FASTQ files matched '${params.fastq_pattern}' in ${params.fastq_dir}"
            System.exit(1)
        }
        if (r2_list.size() > 0 && r2_list.size() != r1_list.size()) {
            log.error "ERROR: R1/R2 file count mismatch (${r1_list.size()} R1 vs ${r2_list.size()} R2) — check --fastq_pattern / --fastq_pattern_2"
            System.exit(1)
        }
        if (r2_list.size() > 0) {
            log.info "  [kge-mageck2] Paired-end: ${r1_list.size()} R1 + ${r2_list.size()} R2 files"
        } else {
            log.info "  [kge-mageck2] Single-end: ${r1_list.size()} files (no R2 files matched)"
        }

        // Channels built from the same glob patterns (dataflow .collect()
        // passed as a collection input — same pattern as kge-nextflow,
        // stable across Nextflow 24.x and 26.x)
        def r1_ch = Channel.fromPath("${params.fastq_dir}/${params.fastq_pattern}")
        def r2_ch = r2_list.size() > 0 ? Channel.fromPath("${params.fastq_dir}/${params.fastq_pattern_2}") : null

        if (runFastqc) {
            // one call (DSL2 processes are single-use per workflow) over
            // the combined R1(+R2) channel; each element is a single file
            FASTQC_READS((r2_ch ? r1_ch.concat(r2_ch) : r1_ch).map { it })
        }

        MAGECK_COUNT(
            file(params.library_file, checkIfExists: true),
            (r2_ch ? r1_ch.concat(r2_ch) : r1_ch).collect(),
            params.sample_labels,
            params
        )
        count_table_ch    = MAGECK_COUNT.out.count_table
        count_norm_ch     = MAGECK_COUNT.out.count_normalized
        count_summary_ch  = MAGECK_COUNT.out.count_summary
        run_logs_ch       = run_logs_ch.concat(MAGECK_COUNT.out.log)
        umi_count_collected = MAGECK_COUNT.out.umi_count.collect()

    } else {
        log.error "ERROR: provide one of --count_table, --fastq_dir, or --demultiplex"
        System.exit(1)
    }

    // ==================================================================
    // RRA: iterate over comparison pairs from --comparisons file
    // ==================================================================

    gene_summaries_ch  = Channel.empty()
    sgrna_summaries_ch = Channel.empty()

    if (params.comparisons) {
        // Parse the comparisons TSV: one (treatment, control) pair per row
        comparisons_ch = Channel.fromPath(params.comparisons, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row['treatment'].trim(), row['control'].trim()) }

        // Combine: (count_table) × (treatment, control) → one MAGECK_TEST_RRA per pair
        MAGECK_TEST_RRA(
            count_table_ch.combine(comparisons_ch),
            params
        )

        gene_summaries_ch  = gene_summaries_ch.concat(MAGECK_TEST_RRA.out.gene_summary)
        sgrna_summaries_ch = sgrna_summaries_ch.concat(MAGECK_TEST_RRA.out.sgrna_summary)
        run_logs_ch        = run_logs_ch.concat(MAGECK_TEST_RRA.out.test_log)
    }

    // ==================================================================
    // MLE: requires --design_matrix
    // ==================================================================

    design_matrix_ch = Channel.empty()

    if (params.design_matrix) {
        MAGECK_MLE_ANALYSIS(
            count_table_ch,
            file(params.design_matrix, checkIfExists: true),
            params
        )

        gene_summaries_ch  = gene_summaries_ch.concat(MAGECK_MLE_ANALYSIS.out.gene_summary)
        sgrna_summaries_ch = sgrna_summaries_ch.concat(MAGECK_MLE_ANALYSIS.out.sgrna_summary)
        run_logs_ch        = run_logs_ch.concat(MAGECK_MLE_ANALYSIS.out.mle_log)

        design_matrix_ch = Channel.fromPath(params.design_matrix)
    }

    // Collect all summaries + step logs into lists for KGE_REPORT
    // (step logs are staged next to the gene summaries so the report can
    //  recover each comparison's samples for the sgRNA line plots and
    //  embed them in the combined log download)
    gene_summaries_collected   = gene_summaries_ch.collect()
    sgrna_summaries_collected  = sgrna_summaries_ch.collect()
    run_logs_collected         = run_logs_ch.collect()

    // Design matrix: real file or NO_DESIGN_MATRIX marker for RRA-only mode
    if (params.design_matrix) {
        // Already assigned above — keep it as-is (Channel.fromPath)
    } else {
        design_matrix_ch = Channel.fromPath('assets/NO_DESIGN_MATRIX')
    }

    // ==================================================================
    // KGE Interactive HTML Report
    // ==================================================================

    if (doHtmlReport) {
        KGE_REPORT(
            count_table_ch,
            gene_summaries_collected,
            sgrna_summaries_collected,
            run_logs_collected,
            design_matrix_ch,
            umi_count_collected,
            skipEnrichment,
            params.enrichment_top_n,
            params.enrichment_fdr,
            params.organism,
            params.fdr_threshold,
            params.top_n,
            params.output_prefix
        )
    }
}


// ============================================================================
// HELPERS
// ============================================================================

/**
 * List the files in a directory whose names match a glob pattern
 * (same syntax as --fastq_pattern: *, ? and {a,b,c} alternation).
 * Uses java.nio.file.PathMatcher — no external dependencies.
 */
def listGlob(File dir, String pattern) {
    def matcher = java.nio.file.FileSystems.getDefault().getPathMatcher("glob:${pattern}")
    return dir.listFiles().findAll { matcher.matches(java.nio.file.Path.of(it.name)) }.sort { it.name }
}

// ============================================================================
// INPUT VALIDATION
// ============================================================================

def validate_inputs() {
    def inputModes = [params.count_table ? 1 : 0, params.fastq_dir ? 1 : 0, params.demultiplex.toString().toBoolean() ? 1 : 0].sum()
    if (inputModes == 0) {
        log.error "ERROR: provide one of --count_table, --fastq_dir, or --demultiplex"
        System.exit(1)
    }
    if (inputModes > 1) {
        log.error "ERROR: provide only one input mode: --count_table, --fastq_dir, or --demultiplex"
        System.exit(1)
    }

    if (params.demultiplex.toString().toBoolean()) {
        if (!params.run_folder)    { log.error "ERROR: --run_folder is required with --demultiplex"; System.exit(1) }
        if (!params.sample_sheet)  { log.error "ERROR: --sample_sheet is required with --demultiplex"; System.exit(1) }
        if (!params.library_file)  { log.error "ERROR: --library_file is required when starting from FASTQ"; System.exit(1) }
    }

    if (params.fastq_dir && !params.library_file) {
        log.error "ERROR: --library_file is required when starting from FASTQ"
        System.exit(1)
    }

    if (!params.comparisons && !params.design_matrix) {
        log.error "ERROR: provide at least one of --comparisons (for RRA) or --design_matrix (for MLE)"
        System.exit(1)
    }

    // UMI window validation (only meaningful when counting runs in this launch)
    if (params.umi && params.umi != 'none') {
        if (params.count_table) {
            log.info "  [kge-mageck2] --umi set but starting from a count table — UMI options are ignored (counting already happened)"
        } else {
            if (params.umi == 'firstpair' && (params.umi_start == -1 || params.umi_end == -1)) {
                log.error "ERROR: --umi firstpair requires --umi_start and --umi_end (offsets relative to the guide in R1)"
                System.exit(1)
            }
            if (params.umi == 'secondpair' && (params.umi_start_2 == -1 || params.umi_end_2 == -1)) {
                log.error "ERROR: --umi secondpair requires --umi_start_2 and --umi_end_2 (offsets from the start of R2)"
                System.exit(1)
            }
        }
    }

    if (params.design_matrix && !params.comparisons) {
        log.info "  [kge-mageck2] MLE-only mode (--design_matrix provided, no --comparisons)"
    }
    if (params.comparisons && !params.design_matrix) {
        log.info "  [kge-mageck2] RRA-only mode (--comparisons provided, no --design_matrix)"
    }

    def mode = params.count_table ? 'count_table' : (params.demultiplex.toString().toBoolean() ? 'BCL' : 'FASTQ')
    def analysis = params.comparisons && params.design_matrix ? 'RRA+MLE' : (params.comparisons ? 'RRA' : 'MLE')
    def umi = (params.umi && params.umi != 'none') ? "${params.umi}" : 'none'

    log.info """
    ============================================
      kge-mageck2: UMI-aware CRISPR KO Screen Pipeline
    ============================================
      Input mode     : ${mode}
      Analysis       : ${analysis}
      Comparisons    : ${params.comparisons ?: 'N/A (MLE only)'}
      Design matrix  : ${params.design_matrix ?: 'N/A (RRA only)'}
      UMI            : ${umi}
      HTML report    : ${params.html_report.toString().toBoolean()}
      Output dir     : ${params.output_dir}
    ============================================
    """.stripIndent()
}