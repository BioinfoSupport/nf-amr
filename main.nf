#!/usr/bin/env nextflow

nextflow.preview.output = true

include { MULTIQC            } from './modules/multiqc'
include { ORGANIZE_FILES     } from './modules/organize_files'

//include { ASSEMBLE_READS    } from './workflows/assemble_reads'
include { IDENTITY          } from './modules/identity'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'

include { ASSEMBLY_QC       } from './subworkflows/assembly_qc'
include { LONG_READS        } from './subworkflows/long_reads'
include { SHORT_READS       } from './subworkflows/short_reads'

include { MULTIREPORT       } from './subworkflows/multireport'
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
			

			// MultiQC
			ORGANIZE_FILES(
				Channel.empty().mix(
					ASSEMBLY_QC.out.long_reads_cram_stats.map({meta,file -> [file,"${meta.id}_long.cram.stats"]}),
					ASSEMBLY_QC.out.short_reads_cram_stats.map({meta,file -> [file,"${meta.id}_short.cram.stats"]}),
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
    	multiqc          = Channel.empty() //MULTIQC.out.html
			
			// Assembly QC
    	long_reads_cram        = ASSEMBLY_QC.out.long_reads_cram
    	long_reads_crai        = ASSEMBLY_QC.out.long_reads_crai
    	long_reads_cram_stats  = ASSEMBLY_QC.out.long_reads_cram_stats

			short_reads_cram       = ASSEMBLY_QC.out.short_reads_cram
			short_reads_crai       = ASSEMBLY_QC.out.short_reads_crai
			short_reads_cram_stats = ASSEMBLY_QC.out.short_reads_cram_stats

			// long-reads
			long_reads_qc             = LONG_READS.out.nanoplot
			long_reads_resfinder      = LONG_READS.out.resfinder
			long_reads_plasmidfinder  = LONG_READS.out.plasmidfinder
			
			// short-reads
			short_reads_qc            = SHORT_READS.out.fastqc_html
			short_reads_resfinder     = SHORT_READS.out.resfinder
			short_reads_plasmidfinder = SHORT_READS.out.plasmidfinder

			
    	html_report      = MULTIREPORT.out.html
    	xlsx_report      = MULTIREPORT.out.xlsx
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
	
	long_reads_cram {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/long_reads.cram" }
		mode 'copy'
	}
	long_reads_crai {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/long_reads.cram.crai" }
		mode 'copy'
	}
	long_reads_cram_stats {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/long_reads.cram.stats" }
		mode 'copy'
	}
	short_reads_cram {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/short_reads.cram" }
		mode 'copy'
	}
	short_reads_crai {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/short_reads.cram.crai" }
		mode 'copy'
	}
	short_reads_cram_stats {
		path { x -> x[1] >> "samples/${x[0].id}/input_assembly/short_reads.cram.stats" }
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
	
	
	
	long_reads_qc {
		path { x -> "samples/${x[0].id}/long_reads/qc" }
		mode 'copy'
	}
	long_reads_resfinder {
		path { x -> "samples/${x[0].id}/long_reads/" }
		mode 'copy'
	}
	long_reads_plasmidfinder {
		path { x -> "samples/${x[0].id}/long_reads/" }
		mode 'copy'
	}
	
	short_reads_qc {
		path { x -> "samples/${x[0].id}/short_reads/qc/" }
		mode 'copy'
	}
	short_reads_resfinder {
		path { x -> "samples/${x[0].id}/short_reads/" }
		mode 'copy'
	}
	short_reads_plasmidfinder {
		path { x -> "samples/${x[0].id}/short_reads/" }
		mode 'copy'
	}
	
	
	
}



