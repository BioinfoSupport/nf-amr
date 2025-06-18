#!/usr/bin/env nextflow

nextflow.preview.output = true

include { KRAKEN2_DB        } from './modules/local/kraken2/db'
include { KRAKEN2_CLASSIFY  } from './modules/local/kraken2/classify'

//include { ASSEMBLE_READS    } from './workflows/assemble_reads'
include { IDENTITY          } from './modules/local/identity'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'
include { MULTIREPORT       } from './subworkflows/local/multireport'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.kraken2_db = null
params.fastq_long = null

workflow {
	main:
			// Validate parameters and print summary of supplied ones
			//validateParameters()
			//log.info(paramsSummaryLog(workflow))
			
			// -------------------
			// Prepare databases
			// -------------------
			k2_db = Channel.empty()
			if (params.kraken2_db=="download") {
				k2_db = KRAKEN2_DB()
			} else if (params.kraken2_db) {
				k2_db = Channel.fromPath(params.kraken2_db)
			}
			
			// -------------------
			// Prepare sequences
			// -------------------
			//ch_ss = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
			fql_ch = Channel.empty()
			fq_ch = Channel.empty()
			fa_ch = Channel.empty()
			if (params.fastq_long) {
				fql_ch = Channel.fromPath(params.fastq_long)
						.map({x -> tuple(["id":x.baseName],x)})
			}

			//ASSEMBLE_READS(Channel.empty())
			fa_ch = Channel.fromPath(params.input)
					.map({x -> tuple(["id":x.baseName],x)})

			k2_ch = KRAKEN2_CLASSIFY(k2_db,fa_ch)

			ann_ch = ANNOTATE_ASSEMBLY(fa_ch)
			IDENTITY(fa_ch)
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
			fasta         = IDENTITY.out
      fai           = ann_ch.fai
      kraken2_db    = k2_db
      kraken2       = k2_ch
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
	fasta {
		path { x -> x[1] >> "samples/${x[0].id}/assembly/assembly.fasta" }
		mode 'copy'
	}

	fai {
		path { x -> x[1] >> "samples/${x[0].id}/assembly/assembly.fasta.fai" }
		mode 'copy'
	}
	
	runinfo {
		path { x -> x[1] >> "samples/${x[0].id}/assembly/anninfo.json" }
		mode 'copy'
	}

	orgfinder {
		path { x -> "samples/${x[0].id}/assembly/" }
		mode 'copy'
	}

	kraken2_db {
		path { x -> x[1] >> "db/kraken2" }
		mode 'copy'
	}
	kraken2 {
		path { x -> "samples/${x[0].id}/assembly/" }
		mode 'copy'
	}

	amrfinderplus {
		path { x -> "samples/${x[0].id}/assembly/" }
		mode 'copy'
	}
	
	resfinder {
		path { x -> "samples/${x[0].id}/assembly/" }
		mode 'copy'
	}
	
	mobtyper {
		path { x -> "samples/${x[0].id}/assembly/" }
		mode 'copy'
	}

	plasmidfinder {
		path { x -> "samples/${x[0].id}/assembly/" }
		mode 'copy'
	}

	cgemlst {
		path { x -> "samples/${x[0].id}/assembly/" }
		mode 'copy'
	}

	MLST {
		path { x -> "samples/${x[0].id}/assembly/" }
		mode 'copy'
	}

	prokka {
		path { x -> "samples/${x[0].id}/assembly/" }
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



