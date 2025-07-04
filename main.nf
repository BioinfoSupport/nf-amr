#!/usr/bin/env nextflow

nextflow.preview.output = true

include { MULTIQC            } from './modules/multiqc'
include { ORGANIZE_FILES     } from './modules/organize_files'

include { ASSEMBLE_READS    } from './subworkflows/assemble_reads'
include { IDENTITY          } from './modules/identity'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'

include { ASSEMBLY_QC       } from './subworkflows/assembly_qc'
include { LONG_READS        } from './subworkflows/long_reads'
include { SHORT_READS       } from './subworkflows/short_reads'

include { MULTIREPORT       } from './subworkflows/multireport'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.fastq_long = []
params.fastq_short = []


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
						.fromFilePairs(params.fastq_short,size:-1) { file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') }
						.map({id,x -> [["id":id],x]})
			}
			if (params.input) {
				fa_ch = Channel.fromPath(params.input)
						.map({x -> tuple(["id":x.baseName],x)})
			}
			IDENTITY(fa_ch)

			// -------------------
			// QC
			// -------------------
			LONG_READS(fql_ch)
			SHORT_READS(fqs_ch)
			ASSEMBLY_QC(fa_ch,fql_ch,fqs_ch)

			// Reads assembly
			ASSEMBLE_READS(fql_ch,fqs_ch)	
			

			// MultiQC
			ORGANIZE_FILES(
				Channel.empty().mix(
					ASSEMBLY_QC.out.long_cram_stats.map({meta,file -> [file,"${meta.id}_long.cram.stats"]}),
					ASSEMBLY_QC.out.short_cram_stats.map({meta,file -> [file,"${meta.id}_short.cram.stats"]}),
					LONG_READS.out.nanostat.map({meta,file -> [file,"${meta.id}_long.nanostat"]}),
					SHORT_READS.out.fastqc_zip.map({meta,files -> [files[0],"${meta.id}_short_fastqc.zip"]}),
					SHORT_READS.out.fastqc_zip.map({meta,files -> [files[1],"${meta.id}_short_R2_fastqc.zip"]})
				)
				.collect({x -> [x]})
			)
			MULTIQC(ORGANIZE_FILES.out,file("${moduleDir}/assets/multiqc/config.yml"))


			// -------------------
			// Run assembly tools
			// -------------------
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
			// Input assembly
			fasta            = IDENTITY.out
      fai              = ann_ch.fai
	    runinfo          = ann_ch.runinfo
	    orgfinder        = ann_ch.orgfinder
      amrfinderplus    = ann_ch.amrfinderplus
    	resfinder        = ann_ch.resfinder
    	mobtyper         = ann_ch.mobtyper
    	plasmidfinder    = ann_ch.plasmidfinder
    	cgemlst          = ann_ch.cgemlst
    	MLST             = ann_ch.MLST
    	prokka           = ann_ch.prokka
			
			// Input assembly QC
    	long_cram           = ASSEMBLY_QC.out.long_cram
    	long_crai           = ASSEMBLY_QC.out.long_crai
    	long_cram_stats     = ASSEMBLY_QC.out.long_cram_stats
			short_cram          = ASSEMBLY_QC.out.short_cram
			short_crai          = ASSEMBLY_QC.out.short_crai
			short_cram_stats    = ASSEMBLY_QC.out.short_cram_stats

			// Long-reads
			long_qc             = LONG_READS.out.nanoplot
			long_resfinder      = LONG_READS.out.resfinder
			long_plasmidfinder  = LONG_READS.out.plasmidfinder
			long_flye           = ASSEMBLE_READS.out.long_flye
			long_unicycler      = ASSEMBLE_READS.out.long_unicycler

			// Short-reads
			short_qc            = SHORT_READS.out.fastqc_html
			short_resfinder     = SHORT_READS.out.resfinder
			short_plasmidfinder = SHORT_READS.out.plasmidfinder
			short_spades        = ASSEMBLE_READS.out.short_spades
			short_unicycler     = ASSEMBLE_READS.out.short_unicycler

			// Hybrid assembly
			hybrid_unicycler          = ASSEMBLE_READS.out.hybrid_unicycler
			hybrid_hybracter          = ASSEMBLE_READS.out.hybrid_hybracter
			
			// Summary reports
    	html_report      = MULTIREPORT.out.html
    	xlsx_report      = MULTIREPORT.out.xlsx
    	multiqc          = Channel.empty() //MULTIQC.out.html
}


output {
	fasta {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/assembly.fasta" }
		mode 'copy'
	}

	fai {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/assembly.fasta.fai" }
		mode 'copy'
	}
	
	runinfo {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/anninfo.json" }
		mode 'copy'
	}

	orgfinder {
		path { x -> "samples/${x[0].id}/input_assembly/" }
		mode 'copy'
	}

	amrfinderplus {
		path { x -> "samples/${x[0].id}/input_assembly/" }
		mode 'copy'
	}
	
	resfinder {
		path { x -> "samples/${x[0].id}/input_assembly/" }
		mode 'copy'
	}
	
	mobtyper {
		path { x -> "samples/${x[0].id}/input_assembly/" }
		mode 'copy'
	}

	plasmidfinder {
		path { x -> "samples/${x[0].id}/input_assembly/" }
		mode 'copy'
	}

	cgemlst {
		path { x -> "samples/${x[0].id}/input_assembly/" }
		mode 'copy'
	}

	MLST {
		path { x -> "samples/${x[0].id}/input_assembly/" }
		mode 'copy'
	}

	prokka {
		path { x -> "samples/${x[0].id}/input_assembly/" }
		mode 'copy'
	}
	
	long_cram {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/long_reads.cram" }
		mode 'copy'
	}
	long_crai {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/long_reads.cram.crai" }
		mode 'copy'
	}
	long_cram_stats {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/long_reads.cram.stats" }
		mode 'copy'
	}
	short_cram {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/short_reads.cram" }
		mode 'copy'
	}
	short_crai {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/short_reads.cram.crai" }
		mode 'copy'
	}
	short_cram_stats {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/short_reads.cram.stats" }
		mode 'copy'
	}
	
	
	// -------------------
	// Long-reads
	// -------------------
	long_qc {
		path { x -> "samples/${x[0].id}/long_reads/qc" }
		mode 'copy'
	}
	long_resfinder {
		path { x -> "samples/${x[0].id}/long_reads/" }
		mode 'copy'
	}
	long_plasmidfinder {
		path { x -> "samples/${x[0].id}/long_reads/" }
		mode 'copy'
	}
	long_flye {
		path { x -> "samples/${x[0].id}/long_reads/" }
		mode 'copy'
	}
	long_unicycler {
		path { x -> "samples/${x[0].id}/long_reads/" }
		mode 'copy'
	}

	// -------------------
	// Short-reads
	// -------------------
	short_qc {
		path { x -> "samples/${x[0].id}/short_reads/qc/" }
		mode 'copy'
	}
	short_resfinder {
		path { x -> "samples/${x[0].id}/short_reads/" }
		mode 'copy'
	}
	short_plasmidfinder {
		path { x -> "samples/${x[0].id}/short_reads/" }
		mode 'copy'
	}
	short_spades {
		path { x -> "samples/${x[0].id}/short_reads/" }
		mode 'copy'
	}
	short_unicycler {
		path { x -> "samples/${x[0].id}/short_reads/" }
		mode 'copy'
	}


	hybrid_unicycler {
		path { x -> "samples/${x[0].id}/hybrid/" }
		mode 'copy'
	}
	hybrid_hybracter {
		path { x -> "samples/${x[0].id}/hybrid/" }
		mode 'copy'
	}
	
	// -------------------
	// Summary reports
	// -------------------	
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



