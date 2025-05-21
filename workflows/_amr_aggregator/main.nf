#!/usr/bin/env nextflow

include { RESFINDER     } from '../resfinder'
include { PLASMIDFINDER } from '../plasmidfinder'
include { MLST          } from '../mlst'
include { ORG_MAP       } from '../org/map'
include { ORG_DB        } from '../org/db'


process RUN_ISOLATE_REPORT {
	  //container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '8 GB'
    cpus 4
    input:
    		tuple(val(meta),path(fasta),val(meta_org),path(ani),path(resfinder),path(mlst),path(plasmidfinder))
    output:
        path("amr_report.json"), emit: json
        path("amr_report.html"), emit: html
    script:
    /*
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
		 */
		 """
		 echo Hello > amr_report.html
		 echo Hello > amr_report.json
		 """
}






process RUN_AMR_REPORT {
	  //container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '8 GB'
    cpus 4
    input:
    		val(meta)
    		path("????/fasta")
    		val(meta_org)
    		path("????/ani")
    		path("????/resfinder")
    		path("????/mlst")
    		path("????/plasmidfinder")
    output:
        path("amr_report.json"), emit: json
        path("amr_report.html"), emit: html
    script:
    /*
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
		 */
		 """
		 echo Hello > amr_report.html
		 echo Hello > amr_report.json
		 """
}


workflow AMR_REPORT {
	take:
		fa_ch
		res_ch
		plf_ch
		mlst_ch
		org_ch
	main:
			// Aggregate results joining on assembly id
			aggr_ch = fa_ch
				.join(res_ch,remainder:true)
				.join(mlst_ch,remainder:true)
				.join(plf_ch,remainder:true)
				.join(org_ch,remainder:true)
				.multiMap({meta,fa,res,mlst,plf,meta_org,ani -> 
					meta: meta
					meta_org: meta_org
					fasta: fa
					ani: ani
					resfinder: res
					mlst: mlst
					plasmidfinder: plf
				})
				
			// Collect all results and call reporting
			RUN_AMR_REPORT(
				aggr_ch.meta.collect(flat:false),
				aggr_ch.fasta.collect(flat:false),
				aggr_ch.meta_org.collect(flat:false),
				aggr_ch.ani.collect(flat:false),
				aggr_ch.resfinder.collect(flat:false),
				aggr_ch.mlst.collect(flat:false),
				aggr_ch.plasmidfinder.collect(flat:false)
			)
	emit:
		html_report = RUN_AMR_REPORT.out.html
}


workflow AMR_ANNOTATE {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				// Determine organism by mapping the assembly on organism database
				ORG_DB()
				org_ch = ORG_MAP(fa_ch)
				res_ch = RESFINDER(fa_ch)
				
				// Plasmid typing
				plf_ch = fa_ch
					.join(org_ch,remainder:true)
					.map({meta,fa,meta_org,ani -> [meta, meta_org, fa]})
					| PLASMIDFINDER
					
				// MLST typing
				mlst_ch = fa_ch
					.join(org_ch,remainder:true)
					.map({meta,fa,meta_org,ani -> [meta, meta_org, fa]})
					| MLST

				AMR_REPORT(fa_ch,res_ch,plf_ch,mlst_ch,org_ch)

		emit:
				resfinder     = res_ch     // channel: [ val(meta), path(resfinder) ]
        org_ani       = org_ch     // channel: [ val(meta), val(org_name) ]
        org_db        = ORG_DB.out // channel: path(org_db) ]
				plasmidfinder = plf_ch     // channel: [ val(meta), path(plasmidfinder) ]
				mlst          = mlst_ch    // channel: [ val(meta), path(mlst) ]
				//report        = AMR_REPORT.out.html_report
}
	



