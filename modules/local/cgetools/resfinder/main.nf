
process CGE_RESFINDER_RUN {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna)
    output:
				tuple val(meta), path("*.resfinder.json")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
					python -m resfinder -ifa '${assembly_fna}' -acq -d -o '${prefix}.resfinder/'
					cp '${prefix}.resfinder/'*.json '${prefix}.resfinder.json'
		    """    
}

process CGE_RESFINDER_FORMAT {
	  label 'Rscript'
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path('resfinder.json')
    output:
				path('*.resfinder.rds')
    script:
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
				#!/usr/bin/env Rscript
				library(tidyverse)
				jsonlite::fromJSON("resfinder.json")\$seq_regions |>
					tibble() |>
					unnest_wider(1) |>
					mutate(contig_id = str_replace(query_id," .*","")) |>
					mutate(assembly_id = "${meta.id}") |>
					saveRDS(file="${prefix}.resfinder.rds")
		    """
}

