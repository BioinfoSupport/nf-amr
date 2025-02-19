#!/usr/bin/env nextflow

process CGE_PLASMIDFINDER_RUN {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna), val(org_name)
    output:
				tuple val(meta), path("*.plasmidfinder.json")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
    		def plf_flags = params.organisms[org_name]["plasmidfinder_flags"]
		    """
	  		mkdir -p '${prefix}.plasmidfinder'
		  	plasmidfinder.py -q -x -p /db/plasmidfinder_db ${plf_flags} ${args} -i '${assembly_fna}' -o '${prefix}.plasmidfinder/'
			  cp '${prefix}.plasmidfinder/data.json' '${prefix}.plasmidfinder.json'
		    """    
}

process CGE_PLASMIDFINDER_FORMAT {
	  label 'Rscript'
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path('plasmidfinder.json')
    output:
				path('*.plasmidfinder.rds')
    script:
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
				#!/usr/bin/env Rscript
				library(tidyverse)
				jsonlite::fromJSON("plasmidfinder.json")\$plasmidfinder\$results\$Enterobacteriales\$enterobacteriales |>
		    tibble() |>
		    unnest_wider(1) |>
		    mutate(contig_id = str_replace(contig_name," .*",""))	|>
				mutate(assembly_id = "${meta.id}") |>
				saveRDS(file="${prefix}.plasmidfinder.rds")
		    """
}

workflow PLASMIDFINDER {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				rds_ch = fa_ch 
				  | CGE_PLASMIDFINDER_RUN 
					| CGE_PLASMIDFINDER_FORMAT
		emit:
				rds_ch    // channel: [ path(rds_file) ]
}



