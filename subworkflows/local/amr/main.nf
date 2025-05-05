#!/usr/bin/env nextflow

include { ORG_DB        } from '../../../modules/local/org/db'
include { ORG_DETECT       } from '../../../modules/local/org/detect'

include { AMRFINDERPLUS_UPDATE } from '../../../modules/local/amrfinderplus/update'
include { AMRFINDERPLUS_RUN } from '../../../modules/local/amrfinderplus/run'
include { PROKKA_RUN        } from '../../../modules/local/prokka'

include { RESFINDER_FA_RUN  } from '../../../modules/local/cgetools/resfinder'
include { PLASMIDFINDER_RUN } from '../../../modules/local/cgetools/plasmidfinder'
include { MLST_RUN          } from '../../../modules/local/cgetools/mlst'
include { MOBTYPER_RUN      } from '../../../modules/local/mobsuite/mobtyper'

include { TO_JSON           } from '../../../modules/local/tojson'
//include { RMD_RENDER        } from '../../../modules/local/rmd/render'


params.default_args = [
	'resfinder_args'	: '',
	'mobtyper_args' : '',
	'amrfinderplus_args' : '',
	'plasmidfinder_args' : '',
	'mlst_args' : null,
	'prokka_args'	: '--kingdom Bacteria'
]
params.skip_prokka = true


def get_org(meta) {
		def org_key = 'org_name'
		if (params.containsKey(org_key)) return params[org_key]
    if (meta.containsKey(org_key)) return meta[org_key]
		return null
}

def get_key(meta,key,org_name=null) {
		if (params.containsKey(key)) return params[key]
    if (meta.containsKey(key)) return meta[key]
    if (org_name==null) return params.default_args[key]
    def org_args = params.organisms.containsKey(org_name)?params.organisms[org_name] : [:]
    if (org_args.containsKey(key)) return org_args[key]
    return params.default_args[key]
}

def get_tool_args(tool_name, meta, org_name=null) {
		return get_key(meta,tool_name + "_args", org_name)
}



workflow AMR_REPORT {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
	
	      // ---------------------------------------------------------------------
	      // Tools that can run directly on a FASTA witout specifying an organism
	      // ---------------------------------------------------------------------
				// CGE - RESFINDER
				resfinder_ch = fa_ch
					.map({meta,fasta -> [meta,fasta,get_tool_args('resfinder',meta)]})
					.filter({meta,fasta,args -> args!=null})
	        | RESFINDER_FA_RUN
	
				// NCBI AMRfinder+
				amrfinderplus_db = AMRFINDERPLUS_UPDATE()
				amrfinderplus_ch = AMRFINDERPLUS_RUN(
						fa_ch
						  .map({meta,fasta -> [meta,fasta,get_tool_args('amrfinderplus',meta)]})
							.filter({meta,fasta,args -> args!=null}),
						amrfinderplus_db
				)
	
				// MOBsuite - MOBtyper
				mobtyper_ch = fa_ch
					.map({meta,fasta -> [meta,fasta,get_tool_args('mobtyper',meta)]})
					.filter({meta,fasta,args -> args!=null})
	        | MOBTYPER_RUN


		
	      // ----------------------------------------------------
	      // Tools specific to an organism
	      // ----------------------------------------------------
				// Run organism detection when organism is unknown
				ORG_DB()
				detected_org_ch = ORG_DETECT(fa_ch.filter({meta,fa -> get_org(meta)==null}),ORG_DB.out)
				
				// Update fa_ch to add detected organism
				fa_org_ch = fa_ch
					.join(detected_org_ch.org_name,remainder:true)
					.map({meta,fa,detected_org_name -> [meta,fa,get_org(meta)?:detected_org_name]})
				
				// Plasmid typing
				plf_ch = fa_org_ch
					.map({meta,fa,org_name -> [meta, fa, get_tool_args('plasmidfinder',meta, org_name)]})
					.filter({meta,fasta,args -> args!=null})
					| PLASMIDFINDER_RUN

				// MLST typing
				mlst_ch = fa_org_ch
					.map({meta,fa,org_name -> [meta, fa, get_tool_args('mlst',meta,org_name)]})
					.filter({meta,fasta,args -> args!=null})
					| MLST_RUN
					
				// PROKKA annotations
				prokka_ch = Channel.empty()
				if (!params.skip_prokka) {
						prokka_ch = fa_org_ch
						  .map({meta,fa,org_name -> [meta, fa, get_tool_args('prokka',meta,org_name)]})
						  .filter({meta,fasta,args -> args!=null})
							| PROKKA_RUN
				}

/*
				meta_json_ch = fa_ch
					.join(org_ch.org_name,remainder:true)
					.join(org_ch.org_ani,remainder:true)
					.join(org_ch.org_acc,remainder:true)
					.map({meta,fa,org_name,org_ani,org_acc -> [meta, [meta:[assembly:meta,org:[org_name:org_name,org_ani:org_ani,org_acc:org_acc]]] ]})
					| TO_JSON
*/
				meta_json_ch = Channel.empty()

/*
				// Aggregate isolate annotations
				isolate_ch = fa_ch
					.join(res_ch,remainder:true)
					.join(mlst_ch,remainder:true)
					.join(plf_ch,remainder:true)
					.join(org_ch,remainder:true)
					.join(meta_json_ch,remainder:true)
					.map({meta,fa,res,mlst,plf,meta_org,ani,meta_json -> 
							[meta,fa,meta_json,[ani,res,mlst,plf].findAll({x->x!=null})]
					})
					//report_ch = RMD_RENDER(isolate_ch,file("${moduleDir}/isolate_report.Rmd"))
*/
				report_ch = Channel.empty()
				
		emit:
		    meta_json        = meta_json_ch     // channel: [ val(meta), path(resfinder) ]
		    prokka           = prokka_ch        // channel: [ val(meta), path(prokka) ]
		    amrfinderplus_db = amrfinderplus_db // channel: path(amrfinderplus_db) ]
		    amrfinderplus    = amrfinderplus_ch // channel: val(meta), path(amrfinderplus) ]
				resfinder        = resfinder_ch     // channel: [ val(meta), path(resfinder) ]
				mobtyper         = mobtyper_ch     // channel: [ val(meta), path(mobtyper) ]
        org_ani          = detected_org_ch.all_ani     // channel: [ val(meta), val(org_name) ]
        org_db           = ORG_DB.out // channel: path(org_db) ]
				plasmidfinder    = plf_ch     // channel: [ val(meta), path(plasmidfinder) ]
				mlst             = mlst_ch    // channel: [ val(meta), path(mlst) ]
				report_html      = report_ch  // channel: [ val(meta), path(html) ]
}
	



