// Stage 1: QC metrics and spot/gene filtering
process QC_FILTER {
    tag "$sample_id"
    label 'process_low'
    publishDir path: { "${params.outdir}/${sample_id}/01_qc" }, mode: 'copy', pattern: "*.{png,csv}"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}_qc.h5ad"), emit: h5ad
    path "${sample_id}_qc_summary.csv",                emit: summary
    path "*.png",                                      emit: plots

    script:
    """
    qc_filter.py --input ${h5ad} \\
        --min_counts ${params.min_counts} --min_genes ${params.min_genes} \\
        --max_pct_mt ${params.max_pct_mt} --min_cells ${params.min_cells}
    """
}
