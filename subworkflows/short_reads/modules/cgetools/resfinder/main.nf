process RESFINDER {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    time '1h'
    input:
        tuple val(meta), path(seq)
    output:
				tuple val(meta), path("resfinder/", type: 'dir')
    script:
		    """
				python -m resfinder ${task.ext.args?:''} -ifq ${seq.join(' ')} -acq -d -j resfinder/data.json -o 'resfinder/'
		    """    
}

