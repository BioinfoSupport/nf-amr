#!/usr/bin/env nextflow

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def setAssemblyId(json_file,assembly_id) {
    def slurp = new JsonSlurper().parseText(json_file.text)
    def builder = new JsonBuilder(slurp)
		builder.content.assembly_id = assembly_id
		json_file.text = builder.toPrettyString()
}

process CGE_PLASMIDFINDER {
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

workflow PLASMIDFINDER {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				out_ch = fa_ch 
				  | CGE_PLASMIDFINDER 
				  | map({meta,json -> setAssemblyId(json,meta.id);[meta,json]})
		emit:
				out_ch    // channel: [ val(meta), path(json_file) ]
}



