#!/usr/bin/env nextflow


process RESFINDER {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
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
				python -m resfinder ${task.ext.args?:''} ${input} -acq -d -j resfinder/data.json -o 'resfinder/'
		    """    
}

