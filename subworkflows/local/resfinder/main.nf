#!/usr/bin/env nextflow

params.skip_resfinder = false

process RUN_RESFINDER {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna)
    output:
				tuple val(meta), path("resfinder/", type: 'dir')
    script:
		    """
				python -m resfinder ${task.ext.args?:''} -ifa '${assembly_fna}' -acq -d -j resfinder/data.json -o 'resfinder/'
		    """    
}

workflow RESFINDER {
		take:
	    	fa_ch
		main:
			if (!params.skip_resfinder) {
				out_ch = RUN_RESFINDER(fa_ch)
			} else {
				out_ch = Channel.empty()
			}
		emit:
				out_ch
}

