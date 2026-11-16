// Stage 6: combine all samples into one HTML report + summary table
process REPORT {
    label 'process_low'
    publishDir "${params.outdir}/summary", mode: 'copy'

    input:
    path files, stageAs: 'inputs/*'

    output:
    path "report.html",         emit: html
    path "summary_metrics.csv", emit: metrics
    path "summary_ari.png",     emit: plot

    script:
    """
    report.py --indir inputs
    """
}
