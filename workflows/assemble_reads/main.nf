
include { UNICYCLER    as LONG_UNICYCLER   } from './modules/unicycler'
include { UNICYCLER    as SHORT_UNICYCLER  } from './modules/unicycler'
include { UNICYCLER    as HYBRID_UNICYCLER } from './modules/unicycler'
include { HYBRACTER    as LONG_HYBRACTER   } from './modules/hybracter'
include { HYBRACTER    as HYBRID_HYBRACTER } from './modules/hybracter'
include { SPADES       as SHORT_SPADES     } from './modules/spades'
include { FLYE_MEDAKA  as LONG_FLYE_MEDAKA } from './subworkflows/flye_medaka'
include { PILON_POLISH as PILON_POLISH_ROUND1 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND2 } from './subworkflows/pilon_polish'
include { PILON_POLISH as PILON_POLISH_ROUND3 } from './subworkflows/pilon_polish'


params.long_unicycler = false
params.short_unicycler = false
params.hybrid_unicycler = false
params.long_hybracter = false
params.hybrid_hybracter = false
params.short_spades = false
params.long_flye_medaka = false
params.hybrid_flye_medaka_pilon = false



process HYBRID_FLYE_MEDAKA_PILON_FOLDER {
    input:
    	tuple val(meta),path('flye_medaka'),path('flye_medaka_pilon/03_pilon_round1'),path('flye_medaka_pilon/04_pilon_round2'),path('flye_medaka_pilon/05_pilon_round3')
    output:
    	tuple val(meta),path("flye_medaka_pilon",type: 'dir')
    script:
    """
    	cp --no-dereference flye_medaka/* flye_medaka_pilon/
    """
}

workflow HYBRID_FLYE_MEDAKA_PILON {
	take:
		flye_medaka_ch
		fqs_ch
	main:
		pilon1_ch = PILON_POLISH_ROUND1(
			flye_medaka_ch.map({meta,asm -> [meta,asm/'02_medaka/consensus.fasta']}),
			fqs_ch
		)
		pilon2_ch = PILON_POLISH_ROUND2(
			pilon1_ch.map({meta,asm -> [meta,asm/'pilon.fasta']}),
			fqs_ch
		)
		pilon3_ch = PILON_POLISH_ROUND3(
			pilon2_ch.map({meta,asm -> [meta,asm/'pilon.fasta']}),
			fqs_ch
		)
		HYBRID_FLYE_MEDAKA_PILON_FOLDER(flye_medaka_ch.join(pilon1_ch).join(pilon2_ch).join(pilon3_ch))
	emit:
		HYBRID_FLYE_MEDAKA_PILON_FOLDER.out
}



workflow ASSEMBLE_READS {
		take:
	    	fql_ch    // channel: [ val(meta), path(long_reads) ]
	    	fqs_ch    // channel: [ val(meta), path(short_reads) ]
		main:
				// Short reads only assemblies
				SHORT_SPADES(
					fqs_ch
						.filter({params.short_spades})
						.map({meta,fqs -> [meta,fqs,[]]})
				)
				SHORT_UNICYCLER(
					fqs_ch
						.filter({params.short_unicycler})
						.map({meta,fqs -> [meta,fqs,[]]})
				)
				
				// Long reads only assemblies
				LONG_FLYE_MEDAKA(
					fql_ch
						.filter({params.long_flye_medaka|params.hybrid_flye_medaka_pilon})
				)
				LONG_HYBRACTER(
					fql_ch
						.filter({params.long_hybracter})
						.map({meta,fql -> [meta,[],fql]})
				)
				LONG_UNICYCLER(
					fql_ch
						.filter({params.long_unicycler})
						.map({meta,fql -> [meta,[],fql]})
				)

				// Hybrid assemblies
				HYBRID_HYBRACTER(
					fql_ch.join(fqs_ch)
						.filter({params.hybrid_hybracter})
						.map({meta,fql,fqs -> [meta,fqs,fql]})
				)
				HYBRID_UNICYCLER(
					fql_ch.join(fqs_ch)
						.filter({params.hybrid_unicycler})
						.map({meta,fql,fqs -> [meta,fqs,fql]})
				)
				HYBRID_FLYE_MEDAKA_PILON(
					LONG_FLYE_MEDAKA.out,
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
		    hybrid_flye_medaka_pilon = HYBRID_FLYE_MEDAKA_PILON.out
}

