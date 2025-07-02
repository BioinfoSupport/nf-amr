
process AMRFINDERPLUS_UPDATE {
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.22--hf69ffd2_0'
    memory '4 GB'
    cpus 1
    time '1h'
    output:
		    path('amrfinder_db', type: 'dir')
    script:
    """
    amrfinder_update -d db && mv \$(realpath db/latest) amrfinder_db
    """
}

