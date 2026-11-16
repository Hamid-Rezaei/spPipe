#!/usr/bin/env nextflow
/*
 * spPipe: does spatial information help unsupervised clustering recover the
 * layers of the human dorsolateral prefrontal cortex (10x Visium, spatialLIBD)?
 *
 *   FETCH_DATA -> QC_FILTER -> NORMALIZE -+-> CLUSTER_EXPR ----+-> EVALUATE ------+
 *                                         +-> CLUSTER_SPATIAL -+-> SPATIAL_GENES -+-> REPORT
 *
 * Usage:  nextflow run main.nf -profile conda      (or -profile docker)
 */

include { FETCH_DATA      } from './modules/fetch_data'
include { QC_FILTER       } from './modules/qc_filter'
include { NORMALIZE       } from './modules/normalize'
include { CLUSTER_EXPR    } from './modules/cluster_expr'
include { CLUSTER_SPATIAL } from './modules/cluster_spatial'
include { EVALUATE        } from './modules/evaluate'
include { SPATIAL_GENES   } from './modules/spatial_genes'
include { REPORT          } from './modules/report'

workflow {
    log.info """
    spPipe: spatial vs expression-only clustering of DLPFC layers
    =============================================================
    samplesheet : ${params.samplesheet}
    outdir      : ${params.outdir}
    clustering  : ${params.n_clusters} clusters, ${params.n_pcs} PCs, spatial rings=${params.n_rings}, lambda=${params.lambda}
    """.stripIndent()

    samples = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> row.sample_id.toString().trim() }
        .filter { sid -> sid }

    FETCH_DATA(samples)
    QC_FILTER(FETCH_DATA.out.h5ad)
    NORMALIZE(QC_FILTER.out.h5ad)

    CLUSTER_EXPR(NORMALIZE.out.h5ad)
    CLUSTER_SPATIAL(NORMALIZE.out.h5ad)

    EVALUATE(
        NORMALIZE.out.h5ad
            .join(CLUSTER_EXPR.out.clusters)
            .join(CLUSTER_SPATIAL.out.clusters)
    )
    SPATIAL_GENES(NORMALIZE.out.h5ad.join(CLUSTER_SPATIAL.out.clusters))

    report_inputs = QC_FILTER.out.summary
        .mix(QC_FILTER.out.plots.flatten())
        .mix(NORMALIZE.out.plots.flatten())
        .mix(CLUSTER_EXPR.out.plots.flatten())
        .mix(CLUSTER_SPATIAL.out.plots.flatten())
        .mix(EVALUATE.out.metrics)
        .mix(EVALUATE.out.plots.flatten())
        .mix(SPATIAL_GENES.out.moran)
        .mix(SPATIAL_GENES.out.plots.flatten())
        .collect()
    REPORT(report_inputs)
}

