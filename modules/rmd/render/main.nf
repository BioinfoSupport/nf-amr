
process RMD_RENDER {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
    memory '8 GB'
    cpus 2
    time '30 min'
    input:
    		tuple(val(meta),path(files),val(render_params))
    		each path('report_template.Rmd')
    		each path(extra)
    output:
        tuple(val(meta),path('report.html'),emit:'html')
    script:
				"""
				#!/usr/bin/env Rscript
				rmarkdown::render(
				  knit_root_dir = getwd(),
				  'report_template.Rmd',
				  ${if (render_params==null) {""} else {"params = list(${render_params}),"}}
					output_dir = getwd(),
					output_file = "report.html"
				)
				"""
}
