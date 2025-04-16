#!/usr/bin/env nextflow

params.skip_resfinder = false

process RUN_RESFINDER_FA {
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

process RUN_RESFINDER_FQ {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(fqs)
    output:
				tuple val(meta), path("resfinder/", type: 'dir')
    script:
    		def ifq = fqs.reduce({v0,x -> v0 + " '${x}'"})
		    """
				python -m resfinder ${task.ext.args?:''} -ifq ${ifq} -acq -d -j resfinder/data.json -o 'resfinder/'
		    """    
}

workflow RESFINDER {
		take:
	    	fa_ch // [meta,fasta] or [meta,[fq]] or [meta,[fq1,fq2]]
		main:
			if (params.skip_resfinder) {
				out_ch = Channel.empty()
			} else {
				out_ch = RUN_RESFINDER_FA(fa_ch)
			}
		emit:
				out_ch
}

