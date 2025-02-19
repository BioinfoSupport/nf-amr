#!/usr/bin/env nextflow

include { CGE_RESFINDER_RUN } from '../../../modules/local/cgetools/resfinder'
include { CGE_RESFINDER_FORMAT } from '../../../modules/local/cgetools/resfinder'

workflow RESFINDER {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				res_ch = fa_ch 
				  | CGE_RESFINDER_RUN 
					| CGE_RESFINDER_FORMAT
		emit:
				res_ch    // channel: [ val(meta), path(resfinder_json) ]
}



