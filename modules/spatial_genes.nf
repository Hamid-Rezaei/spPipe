// Stage 5: spatially variable genes (Moran's I), layer markers, domain markers
process SPATIAL_GENES {
    tag "$sample_id"
    label 'process_low'
    publishDir path: { "${params.outdir}/${sample_id}/05_spatial_genes" }, mode: 'copy', pattern: "*.{png,csv}"

    input:
    tuple val(sample_id), path(h5ad), path(spatial_csv)

    output:
    path "${sample_id}_moran_top.csv",      emit: moran
    path "${sample_id}_domain_markers.csv", emit: markers
    path "*.png",                           emit: plots

    script:
    """
    spatial_genes.py --input ${h5ad} --spatial ${spatial_csv} --n_top ${params.n_top_svg}
    """
}
