#!/usr/bin/env nextflow

process PLASMIDFINDER {
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
