#!/usr/bin/env nextflow

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def setAssemblyId(json_file,assembly_id) {
    def slurp = new JsonSlurper().parseText(json_file.text)
    def builder = new JsonBuilder(slurp)
		builder.content.assembly_id = assembly_id
		json_file.text = builder.toPrettyString()
}

process CGE_RESFINDER {
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

workflow RESFINDER {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				out_ch = fa_ch 
				  | CGE_RESFINDER
				  | map({meta,json -> setAssemblyId(json,meta.id);[meta,json]})
		emit:
				out_ch    // channel: [ val(meta), path(json_file) ]
}



