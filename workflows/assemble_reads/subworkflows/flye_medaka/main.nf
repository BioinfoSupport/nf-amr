
include { FLYE             } from './modules/flye'
include { MEDAKA_CONSENSUS } from './modules/medaka/consensus'

process OUTPUT_FOLDER {
    input:
    	tuple val(meta),path('flye_medaka/01_flye'),path('flye_medaka/02_medaka')
    output:
    	tuple val(meta),path("flye_medaka",type: 'dir')
    script:
    """
    """
}

workflow FLYE_MEDAKA {
	take:
		fql_ch
	main:
		FLYE(fql_ch)
			.join(fql_ch)
			.map({meta,flye,fql -> [meta,flye / 'assembly.fasta',fql]})
			| MEDAKA_CONSENSUS
		
		OUTPUT_FOLDER(FLYE.out.join(MEDAKA_CONSENSUS.out))
	emit:
		OUTPUT_FOLDER.out
}



