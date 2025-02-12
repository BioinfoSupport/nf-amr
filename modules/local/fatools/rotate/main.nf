
process FATOOLS_ROTATE {
	  container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), file(assembly_fna)
    output:
    		tuple val(meta), file("*.rot.fasta")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
		    fa_rotate --out '${prefix}.rot.fasta' '${assembly_fna}'
		    """
}


