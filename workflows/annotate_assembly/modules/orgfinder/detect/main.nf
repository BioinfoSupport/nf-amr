
process ORGFINDER_DETECT {
    container 'registry.gitlab.unige.ch/amr-genomics/orgfinder:v0.2'
    memory '20 GB'
    cpus 8
    time '15 min'
    input:
        tuple val(meta), path('assembly.fna')
    output:
    		tuple val(meta), path("orgfinder",type: 'dir'), emit: orgfinder
        tuple val(meta), env("ORG_NAME"), emit: org_name
        tuple val(meta), env("ORG_ACC"), emit: org_acc
        tuple val(meta), env("ORG_ANI"), emit: ani
    script:
		    """
		    mkdir ./orgfinder && cp /app/db/tax.tsv /app/db/tax.rds /app/db/db.tsv ./orgfinder
		    orgfinder ${task.ext.args?:''} --thread=${task.cpus} --out 'orgfinder/ani.tsv' 'assembly.fna'
		    r -e 'readr::read_tsv("orgfinder/ani.tsv",show_col_types=FALSE) |> 
		    	dplyr::slice_max(ANI,n=1,with_ties = FALSE) |> 
		    	dplyr::mutate(stringr::str_glue("ORG_ACC=\\"{assembly_acc}\\"\\nORG_NAME=\\"{org_name}\\"\\nORG_ANI=\\"{ANI}\\"")) |> 
		    	dplyr::pull() |> 
		    	writeLines("org.env")
		    '
		    source "org.env"
		    """
}

