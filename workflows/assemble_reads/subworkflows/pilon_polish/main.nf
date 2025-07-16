
include { BWA_INDEX } from './modules/bwa/index'
include { BWA_MEM   } from './modules/bwa/mem'
include { PILON     } from './modules/pilon'

workflow PILON_POLISH {
	take:
		fa_ch
		fqs_ch
	main:
		BWA_MEM(BWA_INDEX(fa_ch).join(fqs_ch))
		PILON(fa_ch.join(BWA_MEM.out.bam).join(BWA_MEM.out.bai))
	emit:
		PILON.out
}



