#!/usr/bin/env nextflow

nextflow.preview.output = true

include { MINIMAP2_ALIGN_ONT } from './modules/local/minimap2/align_ont'
include { SAMTOOLS_STATS     } from './modules/local/samtools/stats'
include { NANOPLOT           } from './modules/local/nanoplot'
include { MULTIQC            } from './modules/local/multiqc'


//include { ASSEMBLE_READS    } from './workflows/assemble_reads'
include { IDENTITY          } from './modules/local/identity'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'
include { MULTIREPORT       } from './subworkflows/local/multireport'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.fastq_long = null
params.fastq_short = null

workflow {
	main:
			// Validate parameters and print summary of supplied ones
			//validateParameters()
			//log.info(paramsSummaryLog(workflow))
			
			// -------------------
			// Prepare sequences
			// -------------------
			//ch_ss = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
			fql_ch = Channel.empty()
			fqs_ch = Channel.empty()
			fa_ch = Channel.empty()
			if (params.fastq_long) {
				fql_ch = Channel.fromPath(params.fastq_long)
						.map({x -> tuple(["id":x.name.replaceAll(/(\.fastq|.fq)\.gz$/, "")],x)})
			}
			if (params.fastq_short) {
				fqs_ch = Channel
						.fromFilePairs(params.fastq_short,size:-1) { file -> file.name.replaceAll(/(.*)(_R?[12])?(_[0-9][0-9][0-9])?(\.fastq|.fq)\.gz$/, '$1') }
						.map({id,x -> [["id":id],x]})
			}
			if (params.input) {
				fa_ch = Channel.fromPath(params.input)
						.map({x -> tuple(["id":x.baseName],x)})
			}



			// -------------------
			// Run long read tools
			// -------------------
			NANOPLOT(fql_ch)
			MINIMAP2_ALIGN_ONT(fa_ch.join(fql_ch))
			SAMTOOLS_STATS(MINIMAP2_ALIGN_ONT.out.cram)
			MULTIQC(SAMTOOLS_STATS.out.map({x,y->y}),file("${moduleDir}/assets/multiqc/config.yml"))
			

			// -------------------
			// Run short read tools
			// -------------------
			//FASTQC(fqs_ch)
			//BWA_MEM(fa_ch.join(fqs_ch))

			

			// -------------------
			// Run assembly tools
			// -------------------
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
    	multiqc       = MULTIQC.out.html
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
	multiqc {
		path { x -> "./" }
		mode 'copy'
	}
}



