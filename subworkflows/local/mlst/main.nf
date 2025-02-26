#!/usr/bin/env nextflow

import groovy.json.JsonSlurper
import groovy.json.JsonBuilder

def setAssemblyId(json_file,assembly_id) {
    def slurp = new JsonSlurper().parseText(json_file.text)
    def builder = new JsonBuilder(slurp)
		builder.content.assembly_id = assembly_id
		json_file.text = builder.toPrettyString()
}

process CGE_MLST {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna), val(org_name)
    output:
    		tuple val(meta), path("*.mlst.json")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
    		def mlst_flags = params.organisms[org_name]["mlst_flags"]
		    """
					mkdir -p '${prefix}.mlst/tmp'
					mlst.py -p '/db/mlst_db' -i '${assembly_fna}' -o '${prefix}.mlst' ${mlst_flags} ${args} --tmp_dir '${prefix}.mlst/tmp' -x &> ${prefix}.mlst/tmp/mslt.log
					cp '${prefix}.mlst/data.json' '${prefix}.mlst.json'
		    """
}

workflow MLST {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				out_ch = fa_ch 
				  | CGE_MLST
				  | map({meta,json -> setAssemblyId(json,meta.id);[meta,json]})
		emit:
				out_ch    // channel: [ val(meta), path(json_file) ]
}

