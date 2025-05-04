#!/usr/bin/env nextflow

include { get_tool_args } from '../../../modules/local/functions.nf'

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
	'prokka_args'	: '--kingdom Bacteria',
	'amrfinderplus_args' : ''
]

params.skip_prokka = true



workflow AMR_REPORT {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				// Determine organism by mapping the assembly on organism database
				ORG_DB()
				org_ch = ORG_DETECT(fa_ch,ORG_DB.out)
			
				// CGE - RESFINDER
				res_ch = RESFINDER_FA_RUN(fa_ch)
				
				// NCBI AMRfinder+
				amrfinderplus_db = AMRFINDERPLUS_UPDATE()
				amrfinderplus_ch = AMRFINDERPLUS_RUN(
						fa_ch.map({meta,fasta -> [meta,fasta,get_tool_args('amrfinderplus',meta,params.organisms,params.default_args)]}),
						amrfinderplus_db
				)
				
				// Plasmid typing
				plf_ch = fa_ch
					.join(org_ch.org_name,remainder:true)
					.map({meta,fa,org_name -> [meta, fa, get_tool_args('plasmidfinder',meta,params.organisms,params.default_args)]})
					.filter({meta,fasta,args -> args!=null})
					| PLASMIDFINDER_RUN

				// MOBsuite - MOBtyper
				mob_ch = fa_ch
					.map({meta,fasta -> [meta,fasta,get_tool_args('mobtyper',meta)]})
					.filter({meta,fasta,args -> args!=null})
          | MOBTYPER_RUN

				// MLST typing
				mlst_ch = fa_ch
					.join(org_ch.org_name,remainder:true)
					.map({meta,fa,org_name -> [meta, fa, get_tool_args('mlst',meta,params.organisms,params.default_args,null)]})
					.filter({meta,fasta,args -> args!=null})
					| MLST_RUN
					
				// PROKKA annotations
				prokka_ch = Channel.empty()
				if (!params.skip_prokka) {
						prokka_ch = fa_ch
						  .join(org_ch.org_name,remainder:true)
						  .map({meta,fa,org_name -> [meta, fa, get_tool_args('prokka',meta,params.organisms,params.default_args,null)]})
						  .filter({meta,fasta,args -> args!=null})
							| PROKKA_RUN
				}

				meta_json_ch = fa_ch
					.join(org_ch.org_name,remainder:true)
					.join(org_ch.org_ani,remainder:true)
					.join(org_ch.org_acc,remainder:true)
					.map({meta,fa,org_name,org_ani,org_acc -> [meta, [meta:[assembly:meta,org:[org_name:org_name,org_ani:org_ani,org_acc:org_acc]]] ]})
					| TO_JSON

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
				resfinder        = res_ch     // channel: [ val(meta), path(resfinder) ]
				mobtyper         = mob_ch     // channel: [ val(meta), path(mobtyper) ]
        org_ani          = org_ch.all_ani     // channel: [ val(meta), val(org_name) ]
        org_db           = ORG_DB.out // channel: path(org_db) ]
				plasmidfinder    = plf_ch     // channel: [ val(meta), path(plasmidfinder) ]
				mlst             = mlst_ch    // channel: [ val(meta), path(mlst) ]
				report_html      = report_ch  // channel: [ val(meta), path(html) ]
}
	



