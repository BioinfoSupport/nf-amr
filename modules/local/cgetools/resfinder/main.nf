#!/usr/bin/env nextflow


process RESFINDER_FA_RUN {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path('assembly.fna'), val(args)
    output:
				tuple val(meta), path("resfinder/", type: 'dir')
    script:
		    """
				python -m resfinder ${task.ext.args?:''} ${args} -ifa 'assembly.fna' -acq -d -j resfinder/data.json -o 'resfinder/'
		    """    
}

process RESFINDER_FQ_RUN {
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

