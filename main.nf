#!/usr/bin/env nextflow

include { AMR_ANNOTATE     } from './subworkflows/local/amr_annotate'
//include { FATOOLS_ROTATE   } from './modules/local/fatools/rotate'
//include { FATOOLS_REHEADER } from './modules/local/fatools/reheader'

workflow {
	main:
			println("""
			 8888b.  88888b.d88b.  888d888 
			    "88b 888 "888 "88b 888P"   
			.d888888 888  888  888 888     
			888  888 888  888  888 888     
			"Y888888 888  888  888 888     
			""")
			
			fa_ch = Channel.fromPath(params.input)
					.map({x -> tuple(["id":x.baseName],x)})
			amr_ch = AMR_ANNOTATE(fa_ch)
			amr_ch.resfinder.view()
	publish:
			amr_ch.resfinder >> 'resfinder'
}





