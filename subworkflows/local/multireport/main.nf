
include { RMD_RENDER        } from '../../../modules/local/rmd/render'
include { COLLECT_FILES     } from '../../../modules/local/collect_files'

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
			COLLECT_FILES(
				"assemblies_annotations",
				fa_ch.map({meta,file -> [meta,file,"assembly.fasta"]})
					.concat(fai_ch.map({meta,file -> [meta,file,"assembly.fasta.fai"]}))
					.concat(runinfo_ch.map({meta,file -> [meta,file,"runinfo.json"]}))
					.concat(orgfinder_ch)
					.concat(amrfinderplus_ch)
					.concat(resfinder_ch)
					.concat(mobtyper_ch)
					.concat(plasmidfinder_ch)
					.concat(cgemlst_ch)
					.concat(MLST_ch)
					.concat(prokka_ch)
					.collect({x -> [x]})
			)
			RMD_RENDER(
				COLLECT_FILES.out.map({x -> ["multireport.html",x,"indir='./assemblies_annotations'"]}),
				file("${moduleDir}/assets/multireport.Rmd"),
				file("${moduleDir}/assets/lib_typing.R")
			)
			
		emit:
	    html    = RMD_RENDER.out  // channel: [ meta, path(html) ]
}




