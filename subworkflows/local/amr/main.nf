#!/usr/bin/env nextflow

include { RESFINDER } from '../resfinder'
include { PLASMIDFINDER } from '../plasmidfinder'
include { MLST } from '../mlst'
include { ORGANISM_MAP      } from '../../../modules/local/organism/org_map'
include { ORGANISM_DB      } from '../../../modules/local/organism/org_db'
//include { REPORTING_AMR     } from '../../../modules/local/reporting/amr'

workflow AMR_ANNOTATE {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				res_ch = RESFINDER(fa_ch)
				
				// Determine organism by mapping the assembly on organism database
				ORGANISM_DB()
				org_ch = ORGANISM_MAP(fa_ch)
				
				// Determine MLST
				mlst_ch = fa_ch
						.join(org_ch.map({meta,org_name,map,env -> [meta,org_name]}))
						.filter({meta,fasta,org_name -> 
								params.organisms.containsKey(org_name) 
								&& params.organisms[org_name].containsKey("mlst_flags")
								&& params.organisms[org_name]["mlst_flags"]
						})
						| MLST

				// Perform plasmid type
				plf_ch = fa_ch
				    .join(org_ch.map({meta,org_name,map,env -> [meta,org_name]}))
						.filter({meta,fasta,org_name -> 
								params.organisms.containsKey(org_name) 
								&& params.organisms[org_name].containsKey("plasmidfinder_flags")
								&& params.organisms[org_name]["plasmidfinder_flags"]
						})
						| PLASMIDFINDER

				res_ch.view()
				mlst_ch.view()
				plf_ch.view()
	
				//REPORTING_AMR(aggr_results_ch.collect(flat:false))

		emit:
				resfinder     = res_ch    // channel: [ path(resfinder_rds) ]
        org_map       = org_ch    // channel: [ val(meta), val(org_name) ]
        org_db        = ORGANISM_DB.out // channel: path(org_db) ]
				plasmidfinder = plf_ch    // channel: [ path(plasmidfinder_rds) ]
				mlst          = mlst_ch   // channel: [ path(mlst_rds) ]
				//report        = REPORTING_AMR.out.report
}
	



