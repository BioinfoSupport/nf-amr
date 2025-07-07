
include { UNICYCLER as UNICYCLER_LONG   } from '../../modules/assembly/unicycler'
include { UNICYCLER as UNICYCLER_SHORT  } from '../../modules/assembly/unicycler'
include { UNICYCLER as UNICYCLER_HYBRID } from '../../modules/assembly/unicycler'
include { HYBRACTER as HYBRACTER_LONG   } from '../../modules/assembly/hybracter'
include { HYBRACTER as HYBRACTER_HYBRID } from '../../modules/assembly/hybracter'
include { SPADES    as SPADES_SHORT     } from '../../modules/assembly/spades'
include { FLYE      as FLYE_LONG        } from '../../modules/assembly/flye'


params.unicycler_long = false
params.unicycler_short = false
params.unicycler_hybrid = false
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
				SPADES_SHORT(fqs_ch.filter({params.spades_short}).map({meta,fqs -> [meta,fqs,[]]}))
				UNICYCLER_SHORT(fqs_ch.filter({params.unicycler_short}).map({meta,fqs -> [meta,fqs,[]]}))
				
				// Long reads only assemblies
				FLYE_LONG(fql_ch.filter({params.flye_long}))
				HYBRACTER_LONG(fql_ch.filter({params.hybracter_long}).map({meta,fql -> [meta,[],fql]}))
				UNICYCLER_LONG(fql_ch.filter({params.unicycler_long}).map({meta,fql -> [meta,[],fql]}))

				// Hybrid assemblies
				HYBRACTER_HYBRID(fql_ch.join(fqs_ch).filter({params.hybracter_hybrid}).map({meta,fql,fqs -> [meta,fqs,fql]}))
				UNICYCLER_HYBRID(fql_ch.join(fqs_ch).filter({params.unicycler_hybrid}).map({meta,fql,fqs -> [meta,fqs,fql]}))

				// TODO: Run assemblies individual QC 
				// TODO: Run QC summary report
		emit:
		    short_spades     = SPADES_SHORT.out
		    short_unicycler  = UNICYCLER_SHORT.out
		    
		    long_flye        = FLYE_LONG.out
		    long_unicycler   = UNICYCLER_LONG.out
		    long_hybracter   = HYBRACTER_LONG.out
		    
		    hybrid_unicycler = UNICYCLER_HYBRID.out
		    hybrid_hybracter = HYBRACTER_HYBRID.out
}

