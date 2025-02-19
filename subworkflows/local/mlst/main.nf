

process CGE_MLST_RUN {
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


process CGE_MLST_FORMAT {
	  label 'Rscript'
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path('mlst.json')
    output:
				path('*.mlst.rds')
    script:
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
				#!/usr/bin/env Rscript
				library(tidyverse)
				jsonlite::fromJSON("mlst.json")\$mlst\$results\$sequence_type |>
		    enframe(value="mlst_type") |>
				mutate(assembly_id = "${meta.id}") |>
				saveRDS(file="${prefix}.mlst.rds")
		    """
}


workflow MLST {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				rds_ch = fa_ch 
				  | CGE_MLST_RUN 
					| CGE_MLST_FORMAT
		emit:
				rds_ch    // channel: [ path(rds_file) ]
}

