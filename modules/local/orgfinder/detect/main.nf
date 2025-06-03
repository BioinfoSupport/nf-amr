
process ORGFINDER_FASTANI_RUN {
		container 'docker.io/staphb/fastani:1.34'
    memory '20.GB'
    cpus 4
		input:
        tuple val(meta), path('assembly.fna')
        each path('org_db')
    output:
        tuple val(meta), path("fastANI.tsv")
    script:
		    """
		    ls org_db/fna/*.fna > refList.txt
				fastANI --threads ${task.cpus} --refList refList.txt --query assembly.fna --output fastANI.tsv
		    """
}


process ORGFINDER_FASTANI_REFORMAT {
    container 'registry.gitlab.unige.ch/amr-genomics/rscript:main'
    memory '2.GB'
    cpus 1
		input:
        tuple val(meta), path('fastANI.tsv')
        each path('org_db')
    output:
        tuple val(meta), path("org.ani"), emit: all_ani
        tuple val(meta), env("ORG_NAME"), emit: org_name
        tuple val(meta), env("ORG_ACC"), emit: org_acc
        tuple val(meta), env("ORG_ANI"), emit: org_ani
    script:
		    """
		    #!/usr/bin/env Rscript
		    db <- readr::read_tsv("org_db/db.tsv",col_types = "ccc")
				ani <- readr::read_tsv('fastANI.tsv',col_names = c('query','ref','ANI','bi_frag','query_frag')) |>
						mutate(assembly_acc = str_replace(basename(ref),".fna$",""),ref=NULL) |>
						left_join(db,by="assembly_acc",relationship="many-to-one") |>
						relocate(query,assembly_acc)
				readr::write_tsv(ani,file='org.ani')
		    ani |> 
		    	dplyr::slice_max(ANI,n=1,with_ties = FALSE) |>
          group_walk(~{Sys.setenv(ORG_ANI=.x$ANI,ORG_ACC=.x$assembly_acc,ORG_NAME=.x$org_name)})
		    """
}


workflow ORGFINDER_DETECT {
    take:
    		fa_ch    // channel: [ val(meta), path(assembly_fna) ]
    		org_db   // path to orgfinder_db
    main:
		    fastani = ORGFINDER_FASTANI_RUN(fa_ch,org_db)
		    orgfinder = ORGFINDER_FASTANI_REFORMAT(fastani.all_ani,org_db)
		emit:
		    all_ani   = orgfinder.all_ani  // channel: [ val(meta), path(resfinder) ]
		    org_name  = orgfinder.org_name // channel: [ val(meta), path(fai) ]
		    org_acc   = orgfinder.org_acc  // channel: [ val(meta), path(prokka) ]
		    org_ani   = orgfinder.org_ani  // channel: path(amrfinderplus_db) ]
}



