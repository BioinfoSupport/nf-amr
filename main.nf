#!/usr/bin/env nextflow

nextflow.preview.output = true

include { SAMTOOLS_FASTQ     } from './modules/samtools/fastq'
include { MULTIQC            } from './modules/multiqc'
include { ORGANIZE_FILES     } from './modules/organize_files'

include { IDENTITY          } from './modules/identity'
include { ASSEMBLE_READS    } from './workflows/assemble_reads'
include { ANNOTATE_ASSEMBLY } from './workflows/annotate_assembly'

include { ASSEMBLY_QC       } from './subworkflows/assembly_qc'
include { LONG_READS        } from './subworkflows/long_reads'
include { SHORT_READS       } from './subworkflows/short_reads'

include { MULTIREPORT       } from './subworkflows/multireport'
include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

params.long_reads = []
params.short_reads = []
params.input_assembly = []

workflow {
	main:
			// Validate parameters and print summary of supplied ones
			validateParameters()
			log.info(paramsSummaryLog(workflow))
			
			// -------------------
			// Prepare SampleSheet
			// -------------------
			ss = [
				asm_ch: Channel.empty(),
				lr_ch: Channel.empty(),
				sr_ch: Channel.empty()
			]
			if (params.samplesheet) {
				SS = Channel.fromList(samplesheetToList(params.samplesheet, "assets/schema_samplesheet.json"))
					.multiMap({x ->
						asm_ch: [[sample_id:x[0].sample_id],x[0].input_assembly]
						lr_ch: [[sample_id:x[0].sample_id],x[0].long_reads]
						sr_ch: [[sample_id:x[0].sample_id],[x[0].short_reads_1,x[0].short_reads_2]]
					})
				ss.lr_ch = SS.lr_ch
				ss.sr_ch = SS.sr_ch
				ss.asm_ch = SS.asm_ch
			} else {
				if (params.long_reads) {
					ss.lr_ch = Channel.fromPath(params.long_reads)
							.map({x -> tuple(["sample_id":x.name.replaceAll(/\.(fastq\.gz|fq\.gz|bam|cram)$/, "")],x)})
				}
				if (params.short_reads) {
					ss.sr_ch = Channel
							.fromFilePairs(params.short_reads,size:-1) { file -> file.name.replaceAll(/_(R?[12])(_001)?\.(fq|fastq)\.gz$/, '') }
							.map({id,x -> [["sample_id":id],x]})
				}
				if (params.input_assembly) {
					ss.asm_ch = Channel.fromPath(params.input_assembly)
							.map({x -> tuple(["sample_id":x.baseName],x)})
				}				
			}
			ss.asm_ch = ss.asm_ch.filter({x,y -> y})
			ss.sr_ch = ss.sr_ch.map({x,y -> [x,y.findAll({v->v})]}).filter({x,y -> y})
			ss.lr_ch = ss.lr_ch.filter({x,y -> y})
			
			
			// Ideally this should not be called but is needed to publish the assembly in the output
			IDENTITY(ss.asm_ch)

			// ------------------------------------------------------------------
			// CONVERT long_reads given in BAM/CRAM format into FASTQ format
			// ------------------------------------------------------------------
			ss.lr_ch = ss.lr_ch.branch({meta,f -> 
				bam: f.name =~ /\.(bam|cram)$/
				fq: true
			})
			ss.lr_ch = ss.lr_ch.fq.mix(SAMTOOLS_FASTQ(ss.lr_ch.bam))
			
			// -------------------
			// QC
			// -------------------
			LONG_READS(ss.lr_ch)
			SHORT_READS(ss.sr_ch)
			ASSEMBLY_QC(ss.asm_ch,ss.lr_ch,ss.sr_ch)

			// Reads assembly
			ASSEMBLE_READS(ss.lr_ch,ss.sr_ch)	
			
			// MultiQC
			ORGANIZE_FILES(
				Channel.empty().mix(
					ASSEMBLY_QC.out.long_bam_stats.map({meta,file -> [file,"${meta.sample_id}_long.bam.stats"]}),
					ASSEMBLY_QC.out.short_bam_stats.map({meta,file -> [file,"${meta.sample_id}_short.bam.stats"]}),
					LONG_READS.out.nanostat.map({meta,file -> [file,"${meta.sample_id}_long.nanostat"]}),
					SHORT_READS.out.fastqc_zip.map({meta,files -> [files[0],"${meta.sample_id}_short_fastqc.zip"]}),
					SHORT_READS.out.fastqc_zip.map({meta,files -> [files[1],"${meta.sample_id}_short_R2_fastqc.zip"]})
				)
				.collect({x -> [x]})
			)
			MULTIQC(ORGANIZE_FILES.out,file("${moduleDir}/assets/multiqc/config.yml"))


			// -------------------
			// Run assembly tools
			// -------------------
			ann_ch = ANNOTATE_ASSEMBLY(ss.asm_ch)
			MULTIREPORT(
				ss.asm_ch,
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
    	long_bam           = ASSEMBLY_QC.out.long_bam
    	long_bai           = ASSEMBLY_QC.out.long_bai
    	long_bam_stats     = ASSEMBLY_QC.out.long_bam_stats
			short_bam          = ASSEMBLY_QC.out.short_bam
			short_bai          = ASSEMBLY_QC.out.short_bai
			short_bam_stats    = ASSEMBLY_QC.out.short_bam_stats

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
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/assembly.fasta" }
		mode 'copy'
	}

	fai {
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/assembly.fasta.fai" }
		mode 'copy'
	}
	
	runinfo {
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/anninfo.json" }
		mode 'copy'
	}

	orgfinder {
		path { x -> "samples/${x[0].sample_id}/input_assembly/" }
		mode 'copy'
	}

	amrfinderplus {
		path { x -> "samples/${x[0].sample_id}/input_assembly/" }
		mode 'copy'
	}
	
	resfinder {
		path { x -> "samples/${x[0].sample_id}/input_assembly/" }
		mode 'copy'
	}
	
	mobtyper {
		path { x -> "samples/${x[0].sample_id}/input_assembly/" }
		mode 'copy'
	}

	plasmidfinder {
		path { x -> "samples/${x[0].sample_id}/input_assembly/" }
		mode 'copy'
	}

	cgemlst {
		path { x -> "samples/${x[0].sample_id}/input_assembly/" }
		mode 'copy'
	}

	MLST {
		path { x -> "samples/${x[0].sample_id}/input_assembly/" }
		mode 'copy'
	}

	prokka {
		path { x -> "samples/${x[0].sample_id}/input_assembly/" }
		mode 'copy'
	}
	
	long_bam {
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/long_reads.bam" }
		mode 'copy'
	}
	long_bai {
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/long_reads.bam.bai" }
		mode 'copy'
	}
	long_bam_stats {
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/long_reads.bam.stats" }
		mode 'copy'
	}
	short_bam {
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/short_reads.bam" }
		mode 'copy'
	}
	short_bai {
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/short_reads.bam.bai" }
		mode 'copy'
	}
	short_bam_stats {
		path { x -> x[1] >> "samples/${x[0].sample_id}/input_assembly/short_reads.bam.stats" }
		mode 'copy'
	}
	
	
	// -------------------
	// Long-reads
	// -------------------
	long_qc {
		path { x -> "samples/${x[0].sample_id}/long_reads/qc" }
		mode 'copy'
	}
	long_resfinder {
		path { x -> "samples/${x[0].sample_id}/long_reads/" }
		mode 'copy'
	}
	long_plasmidfinder {
		path { x -> "samples/${x[0].sample_id}/long_reads/" }
		mode 'copy'
	}
	long_flye {
		path { x -> "samples/${x[0].sample_id}/long_reads/" }
		mode 'copy'
	}
	long_unicycler {
		path { x -> "samples/${x[0].sample_id}/long_reads/" }
		mode 'copy'
	}

	// -------------------
	// Short-reads
	// -------------------
	short_qc {
		path { x -> "samples/${x[0].sample_id}/short_reads/qc/" }
		mode 'copy'
	}
	short_resfinder {
		path { x -> "samples/${x[0].sample_id}/short_reads/" }
		mode 'copy'
	}
	short_plasmidfinder {
		path { x -> "samples/${x[0].sample_id}/short_reads/" }
		mode 'copy'
	}
	short_spades {
		path { x -> "samples/${x[0].sample_id}/short_reads/" }
		mode 'copy'
	}
	short_unicycler {
		path { x -> "samples/${x[0].sample_id}/short_reads/" }
		mode 'copy'
	}


	hybrid_unicycler {
		path { x -> "samples/${x[0].sample_id}/hybrid/" }
		mode 'copy'
	}
	hybrid_hybracter {
		path { x -> "samples/${x[0].sample_id}/hybrid/" }
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



