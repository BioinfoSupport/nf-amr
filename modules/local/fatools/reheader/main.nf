
process FATOOLS_REHEADER {
	  container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna)
    output:
    		tuple val(meta), path("*.hdr.fasta")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
		    fa_reheader --out '${prefix}.hdr.fasta' '${assembly_fna}'
		    """
}


