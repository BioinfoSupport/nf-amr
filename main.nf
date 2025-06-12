#!/usr/bin/env nextflow

nextflow.preview.output = true

//include { ASSEMBLE_READS    } from './workflows/assemble_reads'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'
include { MULTIREPORT       } from './subworkflows/local/multireport'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'


workflow {
	main:
			// Validate parameters and print summary of supplied ones
			//validateParameters()
			//log.info(paramsSummaryLog(workflow))
			
			//ch_input = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
			//ASSEMBLE_READS(Channel.empty())
			
			fa_ch = Channel.fromPath(params.input)
					.map({x -> tuple(["id":x.baseName],x)})
			ann_ch = ANNOTATE_ASSEMBLY(fa_ch)
			
			MULTIREPORT(
				fa_ch,
				ann_ch.fai,
	    	ann_ch.runinfo,
	    	ann_ch.orgfinder,
	    	ann_ch.amrfinderplus,
	    	ann_ch.resfinder,
	    	ann_ch.mobtyper,
	    	ann_ch.plasmidfinder,
	    	ann_ch.cgemlst,
	    	ann_ch.MLST,
	    	ann_ch.prokka
			)

	publish:
      fai           = ann_ch.fai
	    runinfo       = ann_ch.runinfo
	    orgfinder     = ann_ch.orgfinder
      amrfinderplus = ann_ch.amrfinderplus
    	resfinder     = ann_ch.resfinder
    	mobtyper      = ann_ch.mobtyper
    	plasmidfinder = ann_ch.plasmidfinder
    	cgemlst       = ann_ch.cgemlst
    	MLST          = ann_ch.MLST
    	prokka        = ann_ch.prokka
    	html_report   = MULTIREPORT.out.html
    	xlsx_report   = MULTIREPORT.out.xlsx
}


output {

	fai {
		path { x -> x[1] >> "samples/${x[0].id}/assembly.fasta.fai" }
		mode 'copy'
	}
	
	runinfo {
		path { x -> x[1] >> "samples/${x[0].id}/runinfo.json" }
		mode 'copy'
	}

	orgfinder {
		path { x -> "samples/${x[0].id}/" }
		mode 'copy'
	}

	amrfinderplus {
		path { x -> "samples/${x[0].id}/" }
		mode 'copy'
	}
	
	resfinder {
		path { x -> "samples/${x[0].id}/" }
		mode 'copy'
	}
	
	mobtyper {
		path { x -> "samples/${x[0].id}/" }
		mode 'copy'
	}

	plasmidfinder {
		path { x -> "samples/${x[0].id}/" }
		mode 'copy'
	}

	cgemlst {
		path { x -> "samples/${x[0].id}/" }
		mode 'copy'
	}

	MLST {
		path { x -> "samples/${x[0].id}/" }
		mode 'copy'
	}

	prokka {
		path { x -> "samples/${x[0].id}/" }
		mode 'copy'
	}
	
	html_report {
		path { x -> x[1] >> "${x[0]}" }
		mode 'copy'
	}
	xlsx_report {
		path { x -> x[1] >> "${x[0]}" }
		mode 'copy'
	}
	
}



