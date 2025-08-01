
include { UNICYCLER    as UNICYCLER_LONG   } from './modules/unicycler'
include { UNICYCLER    as UNICYCLER_SHORT  } from './modules/unicycler'
include { UNICYCLER    as UNICYCLER_HYBRID } from './modules/unicycler'
include { HYBRACTER    as HYBRACTER_LONG   } from './modules/hybracter'
include { HYBRACTER    as HYBRACTER_HYBRID } from './modules/hybracter'
include { SPADES       as SPADES_SHORT     } from './modules/spades'
include { FLYE_MEDAKA  as FLYE_LONG        } from './subworkflows/flye_medaka'
include { PILON_POLISH                     } from './subworkflows/pilon_polish'


params.unicycler_long = false
params.unicycler_short = false
params.unicycler_hybrid = false
params.hybracter_long = false
params.hybracter_hybrid = false
params.spades_short = false
params.flye_long = false
params.flye_hybrid = false


workflow ASSEMBLE_READS {
		take:
	    	fql_ch    // channel: [ val(meta), path(long_reads) ]
	    	fqs_ch    // channel: [ val(meta), path(short_reads) ]
		main:
				// Short reads only assemblies
				SPADES_SHORT(fqs_ch.filter({params.spades_short}).map({meta,fqs -> [meta,fqs,[]]}))
				UNICYCLER_SHORT(fqs_ch.filter({params.unicycler_short}).map({meta,fqs -> [meta,fqs,[]]}))
				
				// Long reads only assemblies
				FLYE_LONG(fql_ch.filter({params.flye_long|params.flye_hybrid}))
				HYBRACTER_LONG(fql_ch.filter({params.hybracter_long}).map({meta,fql -> [meta,[],fql]}))
				UNICYCLER_LONG(fql_ch.filter({params.unicycler_long}).map({meta,fql -> [meta,[],fql]}))

				// Hybrid assemblies
				HYBRACTER_HYBRID(fql_ch.join(fqs_ch).filter({params.hybracter_hybrid}).map({meta,fql,fqs -> [meta,fqs,fql]}))
				UNICYCLER_HYBRID(fql_ch.join(fqs_ch).filter({params.unicycler_hybrid}).map({meta,fql,fqs -> [meta,fqs,fql]}))
				PILON_POLISH(
					FLYE_LONG.out.map({meta,flye -> [meta,flye/'02_medaka/consensus.fasta']}).filter({params.flye_hybrid}),
					fqs_ch.filter({params.flye_hybrid})
				)

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

