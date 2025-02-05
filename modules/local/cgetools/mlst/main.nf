
process CGE_MLST {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna)
    output:
    		tuple val(meta), path("*.mlst.json")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
    		
    		//params.organisms[ params.genome ].containsKey("mlst_flags")
    		if (params.organisms[meta.species]
    		
    		def species = meta.species ?: "N/A"
    		if (species!="N/A") {
			    """
					mkdir -p '${prefix}.mlst/tmp'
					mlst.py -p '/db/mlst_db' -i '${assembly_fna}' -o '${prefix}.mlst' --species ${species} ${args} --tmp_dir '${prefix}.mlst/tmp' -x &> ${prefix}.mlst/tmp/mslt.log
					cp '${prefix}.mlst/data.json' '${prefix}.mlst.json'
			    """
    		} else {
    			error "invalid species"
    		}
}


