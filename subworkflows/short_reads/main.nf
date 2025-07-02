
include { FASTQC             } from '../../modules/fastqc'
//include { RESFINDER      } from '../../modules/cgetools/resfinder'
//include { PLASMIDFINDER  } from '../../modules/cgetools/plasmidfinder'

workflow ONT_READS {
	take:
		fql_ch

	main:
			FASTQC(fqs_ch)
			//RESFINDER(fql_ch,"nanopore")
			//PLASMIDFINDER(fql_ch)
			
	emit:
			fastqc         = FASTQC.out.nanostat
			
			//resfinder     = RESFINDER.out
			//plasmidfinder = PLASMIDFINDER.out
}



