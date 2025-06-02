#!/usr/bin/env nextflow

include { RMD_RENDER        } from '../modules/local/rmd/render'


workflow MULTI_REPORT {
	  take:
	  	amr_ch
		main:
			//RMD_RENDER([["id":"test"],file("assets/rmd")],file("assets/rmd/multireport.Rmd"),[file("assets/rmd/lib_typing.R")])
			println("DONE")
}
	



