#!/usr/bin/env nextflow

process RUN_ORG_MAP {
    container 'registry.gitlab.unige.ch/amr-genomics/species_profiler:main'
    memory '20 GB'
    cpus 4
    input:
        tuple val(meta), path(assembly_fna)
    output:
        tuple val(meta), path("org.ani"), env("ORG_NAME"), env("ORG_ACC"), env("ORG_ANI")
    script:
		    """
		    species_profiler ${task.ext.args?:''} --thread=${task.cpus} --out 'org.ani' '${assembly_fna}'
		    r -e 'readr::read_tsv("org.ani",show_col_types=FALSE) |> 
		    	dplyr::slice_max(ANI,n=1,with_ties = FALSE) |> 
		    	dplyr::mutate(stringr::str_glue("ORG_ACC=\\"{assembly_acc}\\"\\nORG_NAME=\\"{org_name}\\"\\nORG_ANI=\\"{ANI}\\"")) |> 
		    	dplyr::pull() |> 
		    	writeLines("org.env")
		    '
		    source "org.env"
		    """
}

workflow ORG_MAP {
		take:
	    	fa_ch
		main:
				org_ch = fa_ch
					.filter({meta,fa -> !meta.containsKey('org_name')})
					| RUN_ORG_MAP
				org_ch = org_ch
					.map({meta,ani,org_name,org_acc,org_ani -> 
							[meta,['org_name':org_name, 'org_acc':org_acc, 'org_ani':org_ani],ani]
					})
					
				out_ch = fa_ch
					.join(org_ch,remainder:true)
					.map({meta,fa,meta_org,ani -> 
							if (meta.org_name!=null) {
								[meta,['org_name':meta.org_name],ani]
							} else {
								[meta, meta_org, fa]
							}
					})
		emit:
				out_ch // [meta, meta_org, fa]
}

