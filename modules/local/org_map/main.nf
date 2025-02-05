
process ORGANISM_MAP {
    container "registry.gitlab.unige.ch/amr-genomics/species_profiler:main"
    memory '20 GB'
    cpus 4
    input:
        tuple val(meta), path(assembly_fna)
    output:
        tuple val(meta), path(assembly_fna), path("*.org.map"), path("*.org.env"), env("ORG_NAME")
    script:
				def args = task.ext.args ?: ''
    		def prefix = task.ext.prefix ?: (meta.id?:assembly_fna.baseName)
		    """
		    species_profiler ${args} --thread=${task.cpus} --out '${prefix}.org.map' '${assembly_fna}'
		    r -e 'readr::read_tsv("${prefix}.org.map",show_col_types=FALSE) |> 
		    	dplyr::slice_max(ANI,n=1,with_ties = FALSE) |> 
		    	dplyr::mutate(stringr::str_glue("ORG_ACC=\\"{assembly_acc}\\"\\nORG_NAME=\\"{org_name}\\"")) |> 
		    	dplyr::pull() |> 
		    	writeLines("${prefix}.org.env")
		    '
		    source "${prefix}.org.env"
		    """
}


workflow ORGANISM_BEST_HIT {
	take:
	  assembly_fna
	main:
	fa_ch = Channel.fromPath(params.input)
			.map({x -> tuple(["id":x.baseName],x)})
	org_ch = ORGANISM_MAP(fa_ch)
	
	org_ch
		.map({meta,org_map,org_env -> [meta.id,tuple(meta,org_env)]})
		.view()
	
	
	mlst_ch = fa_ch
	 //params.organisms[ params.genome ].containsKey("mlst_flags")
	 
}