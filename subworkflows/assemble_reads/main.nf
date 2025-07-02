
include { HYBRACTER_HYBRID  } from '../../modules/hybracter/hybrid'
include { HYBRACTER_LONG  } from '../../modules/hybracter/long'

workflow ASSEMBLE_READS {
		take:
	    	fql_ch    // channel: [ val(meta), path(long_reads) ]
	    	fqs_ch    // channel: [ val(meta), path(short_reads) ]
		main:
		
				// Short reads only assemblies
				//SPADES(fqs_ch)
				
				// Long reads only assemblies
				//HYBRACTER_LONG(fql_ch)
				//FLYE(fql_ch)
				

				// Hybrid assemblies
				//HYBRACTER_HYBRID(['RH1',file('data/hybracter/RH1/barcode01.fastq.gz'),file('data/hybracter/RH1/RH1_S15_L002_R1_001.fastq.gz'),file('data/hybracter/RH1/RH1_S15_L002_R2_001.fastq.gz')])

		emit:
		    short_spades     = Channel.empty()
				long_hybracter   = Channel.empty()
				long_flye        = Channel.empty()
				long_flye_medaka = Channel.empty()
				hybrid_hybracter = Channel.empty()
}




