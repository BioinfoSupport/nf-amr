
include { UNICYCLER    as LONG_UNICYCLER   } from './modules/unicycler'
include { UNICYCLER    as SHORT_UNICYCLER  } from './modules/unicycler'
include { UNICYCLER    as HYBRID_UNICYCLER } from './modules/unicycler'
include { HYBRACTER    as LONG_HYBRACTER   } from './modules/hybracter'
include { HYBRACTER    as HYBRID_HYBRACTER } from './modules/hybracter'
include { SPADES       as SHORT_SPADES     } from './modules/spades'
include { FLYE_MEDAKA  as LONG_FLYE_MEDAKA } from './subworkflows/flye_medaka'
include { PILON_POLISH                     } from './subworkflows/pilon_polish'


params.long_unicycler = false
params.short_unicycler = false
params.hybrid_unicycler = false
params.long_hybracter = false
params.hybrid_hybracter = false
params.short_spades = false
params.long_flye_medaka = false
params.hybrid_flye_medaka_pilon = false


workflow ASSEMBLE_READS {
		take:
	    	fql_ch    // channel: [ val(meta), path(long_reads) ]
	    	fqs_ch    // channel: [ val(meta), path(short_reads) ]
		main:
				// Short reads only assemblies
				SHORT_SPADES(fqs_ch.filter({params.short_spades}).map({meta,fqs -> [meta,fqs,[]]}))
				SHORT_UNICYCLER(fqs_ch.filter({params.short_unicycler}).map({meta,fqs -> [meta,fqs,[]]}))
				
				// Long reads only assemblies
				LONG_FLYE_MEDAKA(fql_ch.filter({params.long_flye_medaka|params.hybrid_flye_medaka_pilon}))
				LONG_HYBRACTER(fql_ch.filter({params.long_hybracter}).map({meta,fql -> [meta,[],fql]}))
				LONG_UNICYCLER(fql_ch.filter({params.long_unicycler}).map({meta,fql -> [meta,[],fql]}))

				// Hybrid assemblies
				HYBRID_HYBRACTER(fql_ch.join(fqs_ch).filter({params.hybrid_hybracter}).map({meta,fql,fqs -> [meta,fqs,fql]}))
				HYBRID_UNICYCLER(fql_ch.join(fqs_ch).filter({params.hybrid_unicycler}).map({meta,fql,fqs -> [meta,fqs,fql]}))
				PILON_POLISH(
					LONG_FLYE_MEDAKA.out.map({meta,flye -> [meta,flye/'02_medaka/consensus.fasta']}).filter({params.hybrid_flye_medaka_pilon}),
					fqs_ch.filter({params.hybrid_flye_medaka_pilon})
				)

				// TODO: Run assemblies individual QC
				// TODO: Run QC summary report
		emit:
		    short_spades     = SHORT_SPADES.out
		    short_unicycler  = SHORT_UNICYCLER.out
		    
		    long_flye_medaka = LONG_FLYE_MEDAKA.out
		    long_unicycler   = LONG_UNICYCLER.out
		    long_hybracter   = LONG_HYBRACTER.out
		    
		    hybrid_unicycler = HYBRID_UNICYCLER.out
		    hybrid_hybracter = HYBRID_HYBRACTER.out
		    hybrid_flye_medaka_pilon = PILON_POLISH.out
}

