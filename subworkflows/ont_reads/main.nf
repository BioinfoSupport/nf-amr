
include { NANOPLOT           } from '../../modules/local/nanoplot'
include { RESFINDER_RUN      } from '../../modules/local/cgetools/resfinder'
include { PLASMIDFINDER_RUN  } from '../../modules/local/cgetools/plasmidfinder'

workflow ONT_READS {
	take:
		fql_ch

	main:
			NANOPLOT(fql_ch)
			RESFINDER_RUN(fql_ch,"nanopore")
			PLASMIDFINDER_RUN(fql_ch)
			
	emit:
			nanostat      = NANOPLOT.out.nanostat
			nanoplot      = NANOPLOT.out.nanoplot
			resfinder     = RESFINDER_RUN.out
			plasmidfinder = PLASMIDFINDER_RUN.out
}



