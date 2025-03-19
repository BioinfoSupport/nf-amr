#!/usr/bin/env nextflow

nextflow.preview.output = true


include { AMR_REPORT } from './subworkflows/local/amr'

workflow {
	main:
			println("""
			 8888b.  88888b.d88b.  888d888 
			    "88b 888 "888 "88b 888P"   
			.d888888 888  888  888 888     
			888  888 888  888  888 888     
			"Y888888 888  888  888 888     
			""")
			
			fa_ch = Channel.fromPath(params.input)
					.map({x -> tuple(["id":x.baseName],x)})
			amr_ch = AMR_REPORT(fa_ch)
			
	publish:
			amr_ch.resfinder >> 'resfinder'
			amr_ch.org_ani >> 'org_ani'
			amr_ch.org_db >> '.'
			amr_ch.plasmidfinder >> 'plasmidfinder'
			amr_ch.mlst >> 'mlst'
			amr_ch.report_rds >> 'report_rds'
			amr_ch.report_txt >> 'report_txt'
			amr_ch.report_html >> 'report_hrml'
}


output {
	resfinder {
		path({x -> "results/${x[0].id}/"})
		mode 'copy'
	}
	plasmidfinder {
		path({x -> "results/${x[0].id}/"})
		mode 'copy'
	}
	mlst {
		path({x -> "results/${x[0].id}/"})
		mode 'copy'
	}
	org_ani {
		path({x -> "results/${x[0].id}/"})
		mode 'copy'
	}
	report_rds {
		path({x -> "results/${x[0].id}/"})
		mode 'copy'
	}
	report_txt {
		path({x -> "results/${x[0].id}/"})
		mode 'copy'
	}
	report_html {
		path({x -> "results/${x[0].id}/"})
		mode 'copy'
	}	
} // required to publish the output !



