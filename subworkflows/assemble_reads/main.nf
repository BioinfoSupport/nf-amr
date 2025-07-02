
include { UNICYCLER as UNICYCLER_LONG   } from '../../modules/assembly/unicycler'
include { UNICYCLER as UNICYCLER_SHORT  } from '../../modules/assembly/unicycler'
include { UNICYCLER as UNICYCLER_HYBRID } from '../../modules/assembly/unicycler'
include { HYBRACTER as HYBRACTER_LONG   } from '../../modules/assembly/hybracter'
include { HYBRACTER as HYBRACTER_HYBRID } from '../../modules/assembly/hybracter'
include { SPADES    as SPADES_SHORT     } from '../../modules/assembly/spades'
include { FLYE      as FLYE_LONG        } from '../../modules/assembly/flye'


params.unicycler_long = false
params.unicycler_short = false
params.unicycler_hybrid = true
params.hybracter_long = false
params.hybracter_hybrid = false
params.spades_short = false
params.flye_long = false


workflow ASSEMBLE_READS {
		take:
	    	fql_ch    // channel: [ val(meta), path(long_reads) ]
	    	fqs_ch    // channel: [ val(meta), path(short_reads) ]
		main:
				// Short reads only assemblies
				SPADES_SHORT(fqs_ch.map({meta,fqs -> [meta,fqs,[]]}))
				UNICYCLER_SHORT(fqs_ch.map({meta,fqs -> [meta,fqs,[]]}))
				
				// Long reads only assemblies
				FLYE_LONG(fql_ch)
				//HYBRACTER_LONG(fql_ch.map({meta,fql -> [meta,[],fql]})
				UNICYCLER_LONG(fql_ch.map({meta,fql -> [meta,[],fql]}))

				// Hybrid assemblies
				//HYBRACTER_HYBRID(fql_ch.join(fqs_ch).map({meta,fql,fqs -> [meta,fqs,fql]})
				UNICYCLER_HYBRID(fql_ch.join(fqs_ch).map({meta,fql,fqs -> [meta,fqs,fql]}))

				// TODO: Run assemblies individual QC 
				// TODO: Run QC summary report

		emit:
		    short_spades     = SPADES_SHORT.out
		    short_unicycler  = UNICYCLER_SHORT.out
		    
		    long_flye        = FLYE_LONG.out
		    long_unicycler   = UNICYCLER_LONG.out
		    long_hybracter   = Channel.empty()
		    
		    hybrid_unicycler = UNICYCLER_HYBRID.out
		    hybrid_hybracter = Channel.empty()
}

