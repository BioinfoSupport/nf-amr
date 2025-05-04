
process RMD_RENDER {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
    memory '8 GB'
    cpus 1
    input:
    		tuple(val(meta),path(files))
    		each path("report_template.Rmd")
    output:
        tuple(val(meta),path("report.html"))
    script:
				"""
				#!/usr/bin/env Rscript
				p <- list(isolate_dir = getwd())
				rmarkdown::render(
				  knit_root_dir = getwd(),
				  'report_template.Rmd',
					params = p,
					output_dir = getwd(),
					output_file = "report.html"
				)
				"""
}
