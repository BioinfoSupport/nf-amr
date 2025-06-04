#!/usr/bin/env nextflow

nextflow.preview.output = true

params.orgfinder_db = "data/db/org_db"

include { ASSEMBLE_READS } from './workflows/assemble_reads'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'
include { MULTI_REPORT } from './workflows/multi_report'
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
			MULTI_REPORT(amr_ch.runinfo)
			
	publish:
			amr_ch.runinfo >> 'runinfo'
			amr_ch.faidx >> 'faidx'
			amr_ch.resfinder >> 'resfinder'
			amr_ch.mobtyper >> 'mobtyper'
			amr_ch.org_ani >> 'org_ani'
			amr_ch.org_db >> 'org_db'
			amr_ch.plasmidfinder >> 'plasmidfinder'
			amr_ch.mlst >> 'mlst'
			amr_ch.report_html >> 'report_html'
			amr_ch.prokka >> 'prokka'
			amr_ch.amrfinderplus_db >> 'amrfinderplus_db'
			amr_ch.amrfinderplus >> 'amrfinderplus'
			fa_ch >> 'assembly_fa'
}


output {
	org_db {
		path({x -> {filename -> "db/${filename}"}})
		mode 'copy'
	}
	amrfinderplus_db {
		path({x -> {filename -> "db/${filename}"}})
		mode 'copy'
	}
	
	resfinder {
		path({x -> {filename -> "samples/${x[0].id}/${filename}"}})
		mode 'copy'
	}
	mobtyper {
		path({x -> {filename -> "samples/${x[0].id}/${filename}"}})
		mode 'copy'
	}
	prokka {
		path({x -> {filename -> "samples/${x[0].id}/${filename}"}})
		mode 'copy'
	}	
	plasmidfinder {
		path({x -> {filename -> "samples/${x[0].id}/${filename}"}})
		mode 'copy'
	}
	mlst {
		path({x -> {filename -> "samples/${x[0].id}/${filename}"}})
		mode 'copy'
	}
	org_ani {
		path({x -> {filename -> "samples/${x[0].id}/${filename}"}})
		mode 'copy'
	}
	report_html {
		path({x -> {filename -> "samples/${x[0].id}/${filename}"}})
		mode 'copy'
	}
	runinfo {
		path({x -> {filename -> "samples/${x[0].id}/runinfo.json"}})
		mode 'copy'
	}
	assembly_fa {
		path({x -> {filename -> "samples/${x[0].id}/assembly.fasta"}})
		mode 'copy'
	}
	faidx {
		path({x -> {filename -> "samples/${x[0].id}/assembly.fasta.fai"}})
		mode 'copy'
	}
	
	amrfinderplus {
		path({x -> {filename -> "samples/${x[0].id}/${filename}"}})
		mode 'copy'
	}
} // required to publish the output !



