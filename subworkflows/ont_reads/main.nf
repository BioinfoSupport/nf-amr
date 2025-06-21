
include { NANOPLOT       } from '../../modules/local/nanoplot'
include { RESFINDER      } from '../../modules/local/cgetools/resfinder'
include { PLASMIDFINDER  } from '../../modules/local/cgetools/plasmidfinder'

workflow ONT_READS {
	take:
		fql_ch

	main:
			NANOPLOT(fql_ch)
			RESFINDER(fql_ch,"nanopore")
			PLASMIDFINDER(fql_ch)
			
	emit:
			nanostat      = NANOPLOT.out.nanostat
			nanoplot      = NANOPLOT.out.nanoplot
			resfinder     = RESFINDER.out
			plasmidfinder = PLASMIDFINDER.out
}



