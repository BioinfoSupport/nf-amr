
process CGE_RESFINDER {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), file(assembly_fna)
    output:
				tuple val(meta), file("*.resfinder.json")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
					python -m resfinder -ifa '${assembly_fna}' -acq -d -o '${prefix}.resfinder/'
					cp '${prefix}.resfinder/'*.json '${prefix}.resfinder.json'
		    """    
}


