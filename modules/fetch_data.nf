// Stage 0 (setup): download one sample from public repositories and build AnnData
process FETCH_DATA {
    tag "$sample_id"
    label 'process_low'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'terminate' }

    input:
    val sample_id

    output:
    tuple val(sample_id), path("${sample_id}_raw.h5ad"), emit: h5ad

    script:
    """
    fetch_data.py --sample_id ${sample_id}
    """
}
