#!/usr/bin/env nextflow

include { RESFINDER     } from '../../../modules/local/resfinder'
include { PLASMIDFINDER } from '../../../modules/local/plasmidfinder'
include { MLST          } from '../../../modules/local/mlst'
include { ORG_MAP       } from '../../../modules/local/org/map'
include { ORG_DB        } from '../../../modules/local/org/db'


process AMR_REPORT {
	  //container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '8 GB'
    cpus 4
    input:
    		val(meta)
    		path(fastas)
    		path(org_maps)
    		path(org_envs)
    		path(resfinders)
    		path(mlsts)
    		path(plasmidfinders)
    output:
        path("amr_report.json"), emit: json
        path("amr_report.html"), emit: html
    script:
		    def args = task.ext.args ?: ''
		    builder = new groovy.json.JsonBuilder(["version":"v2"])
		    builder.content.meta = meta.collect({x->x.id})
		    builder.content.fasta = fastas.collect({it.toString()})
		    builder.content.org_maps = org_maps.collect({it.toString()})
				builder.content.resfinders = resfinders.collect({it.toString()})
		    builder.content.mlsts = mlsts.collect({it.toString()})
		    builder.content.plasmidfinders = plasmidfinders.collect({it.toString()})
		    println(builder.toPrettyString())
		    
		    """
			    cat <<EOF > amr_report.html
			    ${meta.inject('',{x0,x -> x0 + ' "' + x.id + '"'})}
			    -------------
			    EOF
			    #ls -lha -R ./ >> amr_report.html
			    cp amr_report.html amr_report.json
		    """
}


workflow AMR_ANNOTATE {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				res_ch = RESFINDER(fa_ch)
				
				// Determine organism by mapping the assembly on organism database
				ORG_DB()
				org_ch = ORG_MAP(fa_ch)
				
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
	
				// Aggregate results joining on assembly id
				aggr_ch = fa_ch
				.join(org_ch.map({meta,org_name,map,env -> [meta,["name": org_name,"map": map,"env": env]]}),remainder:true)
				.join(res_ch,remainder:true)
				.join(mlst_ch,remainder:true)
				.join(plf_ch,remainder:true)
				.multiMap({meta,fa,org,res,mlst,plf -> 
					meta: meta
					fasta: fa
					org_map: org.map
					org_env: org.env
					resfinder: res
					mlst: mlst
					plasmidfinder: plf
				})
				
				// Collect all results and call reporting
				/*
				AMR_REPORT(
					aggr_ch.meta.collect(flat:false),
					aggr_ch.fasta.collect(flat:false),
					aggr_ch.org_map.collect(flat:false),
					aggr_ch.org_env.collect(flat:false),
					aggr_ch.resfinder.collect(flat:false),
					aggr_ch.mlst.collect(flat:false),
					aggr_ch.plasmidfinder.collect(flat:false)
				)
				*/

		emit:
				resfinder     = res_ch     // channel: [ path(resfinder_rds) ]
        org_map       = org_ch     // channel: [ val(meta), val(org_name) ]
        org_db        = ORG_DB.out // channel: path(org_db) ]
				plasmidfinder = plf_ch     // channel: [ path(plasmidfinder_rds) ]
				mlst          = mlst_ch    // channel: [ path(mlst_rds) ]
				//report        = AMR_REPORT.out.html
}
	



