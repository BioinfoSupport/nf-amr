
include { FASTQC             } from '../../modules/fastqc'
include { RESFINDER      } from '../../modules/cgetools/resfinder'
include { PLASMIDFINDER  } from '../../modules/cgetools/plasmidfinder'

workflow SHORT_READS {
	take:
		fqs_ch
	main:
			FASTQC(fqs_ch)
			RESFINDER(fqs_ch,'illumina')
			PLASMIDFINDER(fqs_ch)
	emit:
			fastqc_html     = FASTQC.out.html
			fastqc_zip      = FASTQC.out.zip
			resfinder     = RESFINDER.out
			plasmidfinder = PLASMIDFINDER.out
}



