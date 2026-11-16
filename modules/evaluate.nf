// Stage 4: ARI / NMI of both methods against the manual layer annotation
process EVALUATE {
    tag "$sample_id"
    label 'process_low'
    publishDir path: { "${params.outdir}/${sample_id}/04_evaluate" }, mode: 'copy', pattern: "*.{png,csv}"

    input:
    tuple val(sample_id), path(h5ad), path(expr_csv), path(spatial_csv)

    output:
    path "${sample_id}_metrics.csv", emit: metrics
    path "*.png",                    emit: plots

    script:
    """
    evaluate.py --input ${h5ad} --expr ${expr_csv} --spatial ${spatial_csv}
    """
}
