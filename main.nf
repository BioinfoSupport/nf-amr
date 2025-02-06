#!/usr/bin/env nextflow

include { AMR_ANNOTATE } from './subworkflows/local/amr_annotate'
//include {FATOOLS_ROTATE}    from './modules/local/fatools/rotate'
//include {FATOOLS_REHEADER}  from './modules/local/fatools/reheader'

workflow {
	fa_ch = Channel.fromPath(params.input)
			.map({x -> tuple(["id":x.baseName],x)})

	AMR_ANNOTATE(fa_ch)
}





