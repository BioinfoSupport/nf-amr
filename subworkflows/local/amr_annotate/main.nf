#!/usr/bin/env nextflow

include { ORGANISM_MAP      } from '../../../modules/local/org_map'
include { CGE_MLST          } from '../../../modules/local/cgetools/mlst'
include { CGE_PLASMIDFINDER } from '../../../modules/local/cgetools/plasmidfinder'
include { CGE_RESFINDER     } from '../../../modules/local/cgetools/resfinder'

workflow AMR_ANNOTATE {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				res_ch = CGE_RESFINDER(fa_ch)

				// Add org_name metadata
				org_ch = ORGANISM_MAP(fa_ch)
						.map({meta,org_name,x,y -> [meta,org_name]})
				mlst_ch = fa_ch
						.join(org_ch)
						.filter({meta,fasta,org_name -> 
								params.organisms.containsKey(org_name) 
								&& params.organisms[org_name].containsKey("mlst_flags")
								&& params.organisms[org_name]["mlst_flags"]
						})
						| CGE_MLST
				plf_ch = fa_ch
				    .join(org_ch)
						.filter({meta,fasta,org_name -> 
								params.organisms.containsKey(org_name) 
								&& params.organisms[org_name].containsKey("plasmidfinder_flags")
								&& params.organisms[org_name]["plasmidfinder_flags"]
						})
						| CGE_PLASMIDFINDER
		emit:
		    org_map       = org_ch    // channel: [ val(meta), val(org_name) ]
				resfinder     = res_ch    // channel: [ val(meta), path(resfinder_json) ]
				plasmidfinder = plf_ch    // channel: [ val(meta), path(plasmidfinder_json) ]
				mlst          = mlst_ch   // channel: [ val(meta), path(mlst_json) ]
}
	



