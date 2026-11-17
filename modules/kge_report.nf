/*
 * KGE REPORT module — Interactive Plotly HTML Report Generation (mageck2)
 *
 * mageck2 has no `report` CLI subcommand (unlike mageck 0.5.9.5), so this
 * process calls the report entry point directly:
 *
 *     python -c "from mageck2.htmlReport import mageck_report_main; ..."
 *
 * producing the same combined report (`<prefix>.report.html`) as `mageck
 * report` in the original KGE pipeline: interactive Plotly charts, dynamic
 * gene tables, pathway enrichment, and — when a UMI count table is staged
 * — the per-sample UMI (unique molecule) figure.
 *
 * The MAGeCK2 step logs (count + one per RRA comparison + MLE) are staged
 * next to the gene summaries so the report can recover the samples
 * compared in each section (sgRNA line plots) and embed them in the
 * combined log download.
 *
 * All argument values are pre-computed in Groovy and injected as literals
 * into a quoted Python heredoc (no bash variable expansion in the Python
 * code).
 */

process KGE_REPORT {
    publishDir "${params.output_dir}/kge_report", mode: 'copy'
    tag "${report_prefix}"
    label 'kge_report'

    input:
        path count_table_file
        path gene_summaries
        path sgrna_summaries
        path run_logs              // MAGeCK2 step logs (count + per-comparison test [+ MLE])
        path design_matrix_file
        path umi_count             // <prefix>.umi_count.txt from MAGECK_COUNT (empty in count_table mode)
        val skip_enrichment
        val enrichment_top_n
        val enrichment_fdr
        val organism
        val fdr_threshold
        val top_n_value
        val output_prefix

    output:
        path "${report_prefix}.report.html", emit: html_report

    script:
    report_prefix = output_prefix ?: 'mageck2'

    // Design matrix: real file path or Python None for RRA-only mode
    dm_value = (design_matrix_file.name == 'NO_DESIGN_MATRIX') ? 'None' : "'${design_matrix_file}'"

    // Python booleans
    skip_enrichment_py = skip_enrichment ? 'True' : 'False'

    // Stage the UMI count table so the report can draw the UMI figure
    // (htmlReport looks for <output_prefix>.umi_count.txt in its base dir)
    umi_staging = (umi_count && umi_count.size() > 0) ?
        "cp \"${umi_count[0]}\" \"${report_prefix}.umi_count.txt\" && echo \"[KGE_REPORT] Staged UMI count table: ${umi_count[0].name}\"" :
        'echo "[KGE_REPORT] No UMI count table (count-table mode or --umi none)"'

    """
    #! /usr/bin/env bash
    set -euo pipefail

    echo "[KGE_REPORT] Generating interactive HTML report (mageck2 + KGE)..."

    ${umi_staging}

    python - <<'PYEOF'
import glob
from argparse import Namespace
from mageck2.htmlReport import mageck_report_main

gene_summary = sorted(glob.glob('*.gene_summary.txt'))
sgrna_summary = sorted(glob.glob('*.sgrna_summary.txt'))

args = Namespace(
    output_prefix='${report_prefix}',
    count_table='${count_table_file}',
    title=None,
    gene_summary=gene_summary,
    sgrna_summary=sgrna_summary,
    design_matrix=${dm_value},
    fdr_threshold=${fdr_threshold},
    top_n=${top_n_value},
    norm_method='median',
    skip_qc=False,
    skip_results=False,
    enrichment_top_n=${enrichment_top_n},
    organism='${organism}',
    skip_enrichment=${skip_enrichment_py},
    enrichment_fdr=${enrichment_fdr},
    library=None,
    treatment=None,
    control=None,
    output_file=None,
)
mageck_report_main(args)
PYEOF

    echo "[KGE_REPORT] Done."
    ls -la ${report_prefix}.report.html
    """
}