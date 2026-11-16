// Stage 3b: spatially-aware clustering (squidpy spatial graph, BANKSY-style
// neighbour-averaged PCs + Leiden)
process CLUSTER_SPATIAL {
    tag "$sample_id"
    label 'process_low'
    publishDir path: { "${params.outdir}/${sample_id}/03_cluster_spatial" }, mode: 'copy', pattern: "*.{png,csv}"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}_spatial_clusters.csv"), emit: clusters
    path "*.png",                                                   emit: plots

    script:
    """
    cluster_spatial.py --input ${h5ad} --n_pcs ${params.n_pcs} \\
        --n_clusters ${params.n_clusters} --n_neighbors ${params.n_neighbors} \\
        --n_rings ${params.n_rings} --lambda_ ${params.lambda}
    """
}
