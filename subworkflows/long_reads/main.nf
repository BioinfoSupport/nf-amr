
include { NANOPLOT       } from './modules/nanoplot'
include { RESFINDER      } from './modules/cgetools/resfinder'
include { PLASMIDFINDER  } from './modules/cgetools/plasmidfinder'

workflow LONG_READS {
	take:
		fql_ch

	main:
			NANOPLOT(fql_ch)
			RESFINDER(fql_ch)
			PLASMIDFINDER(fql_ch)
			
	emit:
			nanostat      = NANOPLOT.out.nanostat
			nanoplot      = NANOPLOT.out.nanoplot
			resfinder     = RESFINDER.out
			plasmidfinder = PLASMIDFINDER.out
}



