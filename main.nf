#!/usr/bin/env nextflow

nextflow.preview.output = true


include { AMR_ANNOTATE } from './subworkflows/local/amr_annotate'

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
			
	publish:
			amr_ch.resfinder >> 'resfinder'
			amr_ch.org_map >> 'org_map'
			amr_ch.org_db >> 'org_db'
			amr_ch.plasmidfinder >> 'plasmidfinder'
			amr_ch.mlst >> 'mlst'
			//amr_ch.report >> '.'
}

output {} // required to publish the output !



