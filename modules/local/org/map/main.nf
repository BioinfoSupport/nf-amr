
process ORG_MAP {
    container 'registry.gitlab.unige.ch/amr-genomics/species_profiler:main'
    memory '20 GB'
    cpus 4
    input:
        tuple val(meta), path(assembly_fna)
    output:
        tuple val(meta), env("ORG_NAME"), path("*.org.map"), path("*.org.env")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
		    species_profiler ${args} --thread=${task.cpus} --out '${prefix}.org.map' '${assembly_fna}'
		    r -e 'readr::read_tsv("${prefix}.org.map",show_col_types=FALSE) |> 
		    	dplyr::slice_max(ANI,n=1,with_ties = FALSE) |> 
		    	dplyr::mutate(stringr::str_glue("ORG_ACC=\\"{assembly_acc}\\"\\nORG_NAME=\\"{org_name}\\"\\nORG_ANI=\\"{ANI}\\"")) |> 
		    	dplyr::pull() |> 
		    	writeLines("${prefix}.org.env")
		    '
		    source "${prefix}.org.env"
		    """
}
