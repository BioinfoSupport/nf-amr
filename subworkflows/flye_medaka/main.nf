
include { FLYE   } from '../../modules/assembly/flye'
include { MEDAKA_CONSENSUS } from '../../modules/medaka/consensus'

workflow FLYE_MEDAKA {
	take:
		fql_ch
	main:
		FLYE(fql_ch)
			.join(fql_ch)
			.map({meta,flye,fql -> [meta,flye/ 'assembly.fasta',fql]})
			| MEDAKA_CONSENSUS
	emit:
		medaka_consensus = Channel.empty()
}



