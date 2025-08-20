
include { RMD_RENDER        } from '../../modules/rmd/render'
include { ORGANIZE_FILES    } from '../../modules/organize_files'



process MULTITABLE {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
    memory '8 GB'
    cpus 2
    time '30 min'
    input:
    		tuple val(meta),path("db")
    		each path("lib_typing.R")
    output:
        tuple val(meta), path('multitable.xlsx'), emit:"xlsx"
    script:
        """
				#!/usr/bin/env Rscript
				source("lib_typing.R")
				db <- db_load("db")
				tbl <- list(
					assemblies = summarise_assembly(db),
					contigs = summarise_contigs(db),
					resistances = summarise_resistances(db),
					plasmidfinder_hits = summarise_plasmidfinder_hits(db)
				)
				tbl\$contigs <- left_join(tbl\$contigs,select(tbl\$assemblies,assembly_id,species_name,mlst_type),by='assembly_id')
				openxlsx::write.xlsx(tbl,file="multitable.xlsx")
        """
}


workflow MULTIREPORT {
		take:
	    	fa_ch
	    	fai_ch
	    	runinfo_ch
	    	orgfinder_ch
	    	amrfinderplus_ch
	    	resfinder_ch
	    	mobtyper_ch
	    	plasmidfinder_ch
	    	cgemlst_ch
	    	MLST_ch
	    	prokka_ch
	    	
		main:
			ORGANIZE_FILES(
				Channel.empty().mix(
					  fa_ch.map(           {meta,file -> [file,"${meta.sample_id}/input_assembly/assembly.fasta"]}),
						fai_ch.map(          {meta,file -> [file,"${meta.sample_id}/input_assembly/assembly.fasta.fai"]}),
						runinfo_ch.map(      {meta,file -> [file,"${meta.sample_id}/input_assembly/anninfo.json"]}),
						orgfinder_ch.map(    {meta,file -> [file,"${meta.sample_id}/input_assembly/${file.name}"]}),
						amrfinderplus_ch.map({meta,file -> [file,"${meta.sample_id}/input_assembly/${file.name}"]}),
						resfinder_ch.map(    {meta,file -> [file,"${meta.sample_id}/input_assembly/${file.name}"]}),
						plasmidfinder_ch.map({meta,file -> [file,"${meta.sample_id}/input_assembly/${file.name}"]}),
						cgemlst_ch.map(      {meta,file -> [file,"${meta.sample_id}/input_assembly/${file.name}"]}),
						mobtyper_ch.map(     {meta,file -> [file,"${meta.sample_id}/input_assembly/${file.name}"]}),
						MLST_ch.map(         {meta,file -> [file,"${meta.sample_id}/input_assembly/${file.name}"]}),
						prokka_ch.map(       {meta,file -> [file,"${meta.sample_id}/input_assembly/${file.name}"]})
				)
				.collect({[it]})
				.map({["unused",it]})
			)
			RMD_RENDER(
				ORGANIZE_FILES.out.map({m,x -> ["multireport.html",x,"indir='${x}'"]}),
				file("${moduleDir}/assets/multireport.Rmd"),
				file("${moduleDir}/assets/lib_typing.R")
			)
			MULTITABLE(
				ORGANIZE_FILES.out.map({m,x -> ["multitable.xlsx",x]}),
				file("${moduleDir}/assets/lib_typing.R")
			)
			
		emit:
	    html = RMD_RENDER.out  // channel: [ meta, path(html) ]
	    xlsx = MULTITABLE.out  // channel: [ meta, path(xlsx) ]
}




