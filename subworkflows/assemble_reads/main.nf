
include { HYBRACTER as HYBRACTER_LONG   } from '../../modules/hybracter'
include { HYBRACTER as HYBRACTER_HYBRID } from '../../modules/hybracter'
include { SPADES    as SPADES_SHORT     } from '../../modules/spades'
include { FLYE      as FLYE_LONG        } from '../../modules/flye'

workflow ASSEMBLE_READS {
		take:
	    	fql_ch    // channel: [ val(meta), path(long_reads) ]
	    	fqs_ch    // channel: [ val(meta), path(short_reads) ]
		main:
		
				// Short reads only assemblies
				SPADES_SHORT(fqs_ch.map({meta,fqs -> [meta,fqs,[]]}))
				
				// Long reads only assemblies
				//HYBRACTER_LONG(fql_ch.map({meta,fql -> [meta,fql,[]]})
				FLYE_LONG(fql_ch)

				// Hybrid assemblies
				//HYBRACTER_HYBRID(fql_ch.join(fqs_ch).map({meta,fql,fqs -> [meta,fql,fqs]})

		emit:
		    short_spades     = SPADES_SHORT.out
		    long_flye        = FLYE_LONG.out
		    hybrid_hybracter = Channel.empty()
				long_hybracter   = Channel.empty()
}




