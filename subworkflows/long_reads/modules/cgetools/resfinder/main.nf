process RESFINDER {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    time '1h'
    input:
        tuple val(meta), path('long_reads.fastq.gz')
    output:
				tuple val(meta), path("resfinder/", type: 'dir')
    script:
		    """
				python -m resfinder ${task.ext.args?:''} --nanopore -ifq long_reads.fastq.gz -acq -d -j resfinder/data.json -o 'resfinder/'
		    """    
}

