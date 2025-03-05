#!/usr/bin/env nextflow

process MLST {
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

