#!/usr/bin/env nextflow

include { HYBRACTER_HYBRID  } from '../modules/local/hybracter/hybrid'
include { HYBRACTER_LONG  } from '../modules/local/hybracter/long'

workflow ASSEMBLE_READS {
		take:
	    	reads_ch    // channel: [ val(meta), path(long_reads), path(read1), path(read2) ]
		main:
				println("Start read assembly")
				
				// TODO: decide the assembler to use: HYBRACTER_LONG, HYBRACTER_HYBRID or SPADES ?
				HYBRACTER_HYBRID(['RH1',file('data/hybracter/RH1/barcode01.fastq.gz'),file('data/hybracter/RH1/RH1_S15_L002_R1_001.fastq.gz'),file('data/hybracter/RH1/RH1_S15_L002_R2_001.fastq.gz')])
		emit:
				fasta = Channel.empty()
				contig_stat = Channel.empty()
}




