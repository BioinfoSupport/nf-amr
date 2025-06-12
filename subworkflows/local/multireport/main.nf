
include { RMD_RENDER        } from '../../../modules/local/rmd/render'
include { ORGANIZE_FILES    } from '../../../modules/local/organize_files'



process MULTITABLE {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
    memory '8 GB'
    cpus 2
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
				list(
					assemblies = summarise_assembly(db),
					contigs = summarise_contigs(db),
					resistances = summarise_resistances(db)
				) |>
				openxlsx::write.xlsx(file="multitable.xlsx")
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
				Channel.empty()
					.mix(
					  fa_ch.map(           {meta,file -> [file,"${meta.id}/assembly.fasta"]}),
						fai_ch.map(          {meta,file -> [file,"${meta.id}/assembly.fasta.fai"]}),
						runinfo_ch.map(      {meta,file -> [file,"${meta.id}/runinfo.json"]}),
						orgfinder_ch.map(    {meta,file -> [file,"${meta.id}/${file.baseName}"]}),
						amrfinderplus_ch.map({meta,file -> [file,"${meta.id}/${file.baseName}"]}),
						resfinder_ch.map(    {meta,file -> [file,"${meta.id}/${file.baseName}"]}),
						mobtyper_ch.map(     {meta,file -> [file,"${meta.id}/${file.baseName}"]}),
						plasmidfinder_ch.map({meta,file -> [file,"${meta.id}/${file.baseName}"]}),
						cgemlst_ch.map(      {meta,file -> [file,"${meta.id}/${file.baseName}"]}),
						MLST_ch.map(         {meta,file -> [file,"${meta.id}/${file.baseName}"]}),
						prokka_ch.map(       {meta,file -> [file,"${meta.id}/${file.baseName}"]})
					)
					.collect({x -> [x]})
			)
			RMD_RENDER(
				ORGANIZE_FILES.out.map({x -> ["multireport.html",x,"indir='${x}'"]}),
				file("${moduleDir}/assets/multireport.Rmd"),
				file("${moduleDir}/assets/lib_typing.R")
			)
			MULTITABLE(
				ORGANIZE_FILES.out.map({x -> ["multitable.xlsx",x]}),
				file("${moduleDir}/assets/lib_typing.R")
			)
			
		emit:
	    html = RMD_RENDER.out  // channel: [ meta, path(html) ]
	    xlsx = MULTITABLE.out  // channel: [ meta, path(xlsx) ]
}




