// Stage 3a: expression-only clustering (kNN on PCs + Leiden)
process CLUSTER_EXPR {
    tag "$sample_id"
    label 'process_low'
    publishDir path: { "${params.outdir}/${sample_id}/03_cluster_expr" }, mode: 'copy', pattern: "*.{png,csv}"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}_expr_clusters.csv"), emit: clusters
    path "*.png",                                                emit: plots

    script:
    """
    cluster_expr.py --input ${h5ad} --n_pcs ${params.n_pcs} \\
        --n_clusters ${params.n_clusters} --n_neighbors ${params.n_neighbors}
    """
}
