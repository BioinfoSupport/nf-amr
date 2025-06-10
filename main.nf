#!/usr/bin/env nextflow

nextflow.preview.output = true

params.orgfinder_db = "data/db/org_db"

include { ASSEMBLE_READS    } from './workflows/assemble_reads'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

workflow {
	main:
			// Validate parameters and print summary of supplied ones
			validateParameters()
			log.info(paramsSummaryLog(workflow))
			
			//ch_input = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_input.json"))
			//ASSEMBLE_READS(Channel.empty())
			
			fa_ch = Channel.fromPath(params.input)
					.map({x -> tuple(["id":x.baseName],x)})
			amr_ch = ANNOTATE_ASSEMBLY(fa_ch)

	publish:
      results = amr_ch.results
      report_html = amr_ch.report_html
}


output {
	
	results {
		path {x -> x >> "samples"}
		mode 'copy'
	}
	
	report_html {
		path {x -> x[1] >> "${x[0]}"}
		mode 'copy'
	}
	
} // required to publish the output !



