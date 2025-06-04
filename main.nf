#!/usr/bin/env nextflow

nextflow.preview.output = true

params.orgfinder_db = "data/db/org_db"

include { ASSEMBLE_READS    } from './workflows/assemble_reads'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'
include { RMD_RENDER        } from './modules/local/rmd/render'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow {
	main:
			println("""
			 8888b.  88888b.d88b.  888d888 
			    "88b 888 "888 "88b 888P"   
			.d888888 888  888  888 888     
			888  888 888  888  888 888     
			"Y888888 888  888  888 888     
			""")
			
			// Validate parameters and print summary of supplied ones
			validateParameters()
			log.info(paramsSummaryLog(workflow))
			
			//ch_input = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_input.json"))
			//ASSEMBLE_READS(Channel.empty())
			
			fa_ch = Channel.fromPath(params.input)
					.map({x -> tuple(["id":x.baseName],x)})
			amr_ch = ANNOTATE_ASSEMBLY(fa_ch)

	publish:
			amr_ch.results          >> 'results'
			amr_ch.report_html      >> 'report_html'
}


output {

	results {
		path({x -> {filename -> "samples}"}})
		mode 'copy'
	}
	
	report_html {
		path({x -> {filename -> "multireport.html"}})
		mode 'copy'
	}
	
} // required to publish the output !



