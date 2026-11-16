// Stage 2: normalisation, log1p, highly variable genes, PCA
process NORMALIZE {
    tag "$sample_id"
    label 'process_low'
    publishDir path: { "${params.outdir}/${sample_id}/02_normalize" }, mode: 'copy', pattern: "*.{png,h5ad}"

    input:
    tuple val(sample_id), path(h5ad)

    output:
    tuple val(sample_id), path("${sample_id}_norm.h5ad"), emit: h5ad
    path "*.png",                                        emit: plots

    script:
    """
    normalize.py --input ${h5ad} --n_hvg ${params.n_hvg} --n_pcs 50
    """
}
