#!/usr/bin/env nextflow

include { KRAKEN2_DB         } from './modules/local/kraken2/db'
include { KRAKEN2_CLASSIFY   } from './modules/local/kraken2/classify'

params.kraken2_db = null

workflow TEST_KRAKEN2 {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]

		main:
			k2_db = Channel.empty()
			if (params.kraken2_db=="download") {
				k2_db = KRAKEN2_DB()
			} else if (params.kraken2_db) {
				k2_db = Channel.fromPath(params.kraken2_db)
			}
			k2_ch = KRAKEN2_CLASSIFY(k2_db,fa_ch)
			
		emit:
				kraken2 = k2_ch
}
	



