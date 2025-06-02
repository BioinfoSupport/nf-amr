
process ORG_DETECT {
    container 'registry.gitlab.unige.ch/amr-genomics/species_profiler:main'
    memory '20 GB'
    cpus 8
    input:
        tuple val(meta), path('assembly.fna')
        path('org_db')
    output:
        tuple val(meta), path("org.ani"), emit: all_ani
        tuple val(meta), env("ORG_NAME"), emit: org_name
        tuple val(meta), env("ORG_ACC"), emit: org_acc
        tuple val(meta), env("ORG_ANI"), emit: org_ani
    script:
		    """
		    species_profiler ${task.ext.args?:''} --thread=${task.cpus} --db 'org_db' --out 'org.ani' 'assembly.fna'
		    r -e 'readr::read_tsv("org.ani",show_col_types=FALSE) |> 
		    	dplyr::slice_max(ANI,n=1,with_ties = FALSE) |> 
		    	dplyr::mutate(stringr::str_glue("ORG_ACC=\\"{assembly_acc}\\"\\nORG_NAME=\\"{org_name}\\"\\nORG_ANI=\\"{ANI}\\"")) |> 
		    	dplyr::pull() |> 
		    	writeLines("org.env")
		    '
		    source "org.env"
		    """
}

