process RESFINDER {
		label 'cgetools'
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:v2.0"
    memory '4 GB'
    cpus 3
    time '1h'
    input:
        tuple val(meta), path(seq)
        val(mode)
    output:
				tuple val(meta), path("resfinder/", type: 'dir')
    script:
    		def input = ''
    		if (mode=='nanopore') {
    			input = "--nanopore -ifq ${seq.join(' ')}"
    		} else if (mode=='illumina') {
    			input = "-ifq ${seq.join(' ')}"
    		} else {
    			input = "-ifa ${seq.join(' ')}"
    		}
		    """
				python -m resfinder ${task.ext.args?:''} ${input} --kma_threads ${task.cpus} -acq -d -j resfinder/data.json -o 'resfinder/'
		    """    
}

