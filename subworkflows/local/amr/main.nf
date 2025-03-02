#!/usr/bin/env nextflow

include { RESFINDER     } from '../resfinder'
include { PLASMIDFINDER } from '../plasmidfinder'
include { MLST          } from '../mlst'
include { ORG_MAP       } from '../../../modules/local/org/map'
include { ORG_DB        } from '../../../modules/local/org/db'


process AMR_REPORT {
	  //container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '8 GB'
    cpus 4
    input:
    		val(meta)
    		path(fastas)
    		val(org_names)
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
		    
		    //meta.collect({x -> x.id}).collectFile(name:'meta.txt',newLine:true)
		    
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
	

				aggr_ch = fa_ch
				.join(org_ch.map({meta,org_name,map,env -> [meta,[org_name,map,env]]}),remainder:true)
				.join(res_ch,remainder:true)
				.join(mlst_ch,remainder:true)
				.join(plf_ch,remainder:true)
				.multiMap({meta,fa,org,res,mlst,plf -> 
					meta: meta
					fasta: fa
					org_name: org[0]
					org_map: org[1]
					org_env: org[2]
					resfinder: res
					mlst: mlst
					plasmidfinder: plf
				})
				AMR_REPORT(
					aggr_ch.meta.collect(flat:false),
					aggr_ch.fasta.collect(flat:false),
					aggr_ch.org_name.collect(flat:false),
					aggr_ch.org_map.collect(flat:false),
					aggr_ch.org_env.collect(flat:false),
					aggr_ch.resfinder.collect(flat:false),
					aggr_ch.mlst.collect(flat:false),
					aggr_ch.plasmidfinder.collect(flat:false)
				)

		emit:
				resfinder     = res_ch     // channel: [ path(resfinder_rds) ]
        org_map       = org_ch     // channel: [ val(meta), val(org_name) ]
        org_db        = ORG_DB.out // channel: path(org_db) ]
				plasmidfinder = plf_ch     // channel: [ path(plasmidfinder_rds) ]
				mlst          = mlst_ch    // channel: [ path(mlst_rds) ]
				report        = AMR_REPORT.out.html
}
	



