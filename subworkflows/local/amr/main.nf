#!/usr/bin/env nextflow

include { get_tool_args } from '../../../modules/local/functions.nf'

include { AMRFINDERPLUS_UPDATE } from '../../../modules/local/amrfinderplus/update'
include { AMRFINDERPLUS_RUN } from '../../../modules/local/amrfinderplus/run'

include { RESFINDER     } from '../cgetools/resfinder'
include { PLASMIDFINDER } from '../cgetools/plasmidfinder'
include { MLST          } from '../cgetools/mlst'
include { MOBTYPER      } from '../mobsuite/mobtyper'
include { PROKKA        } from '../prokka'
include { ORG_MAP       } from '../org/map'
include { ORG_DB        } from '../org/db'










params.skip_amr_report = true // do not run by default as we are still in debugging phase

process META_TO_JSON {
	input:
		tuple(val(meta),val(json_content))
	output:
		tuple(val(meta),path("meta.json"))
	script:
		def builder = new groovy.json.JsonBuilder(json_content)
"""
cat >> meta.json << __EOF_META_JSON__
${builder.toPrettyString()}
__EOF_META_JSON__
"""
}

process BUILD_RMD_REPORT {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
    memory '8 GB'
    cpus 1
    input:
    		tuple(val(meta),path('assembly.fasta'),path('meta.json'),path(files))
    		each path(report_template)
    output:
        tuple(val(meta),path("report.html"))
    script:
				"""
				#!/usr/bin/env Rscript
				p <- list(isolate_dir = getwd())
				rmarkdown::render(
				  knit_root_dir = getwd(),
				  '${report_template}',
					params = p,
					output_dir = getwd(),
					output_file = "report.html"
				)
				"""
}

workflow AMR_REPORT {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
				// Determine organism by mapping the assembly on organism database
				ORG_DB()
				org_ch = ORG_MAP(fa_ch)
				res_ch = RESFINDER(fa_ch)
				mob_ch = MOBTYPER(fa_ch)
				
				// NCBI AMRfinder+
				amrfinderplus_db = AMRFINDERPLUS_UPDATE()
				amrfinderplus_ch = AMRFINDERPLUS_RUN(
						fa_ch.map({meta,fasta -> [meta,fasta,get_tool_args('amrfinderplus',meta)]}),
						amrfinderplus_db
				)
				
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

				// PROKKA annotations
				prokka_ch = fa_ch
					.join(org_ch,remainder:true)
					.map({meta,fa,meta_org,ani -> [meta, meta_org, fa]})
					| PROKKA

				meta_json_ch = fa_ch
					.join(org_ch,remainder:true)
					.map({meta,fa,meta_org,ani -> [meta, [meta:[assembly:meta,org:meta_org]] ]})
					| META_TO_JSON

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

				if (params.skip_amr_report) {
					report_ch = Channel.empty()
				} else {
					report_ch = BUILD_RMD_REPORT(isolate_ch,file("${moduleDir}/isolate_report.Rmd"))
				}

		emit:
		    meta_json        = META_TO_JSON.out // channel: [ val(meta), path(resfinder) ]
		    prokka           = prokka_ch  // channel: [ val(meta), path(prokka) ]
		    amrfinderplus_db = amrfinderplus_db // channel: path(amrfinderplus_db) ]
		    amrfinderplus    = amrfinderplus_ch // channel: val(meta), path(amrfinderplus) ]
				resfinder        = res_ch     // channel: [ val(meta), path(resfinder) ]
				mobtyper         = mob_ch     // channel: [ val(meta), path(mobtyper) ]
        org_ani          = org_ch     // channel: [ val(meta), val(org_name) ]
        org_db           = ORG_DB.out // channel: path(org_db) ]
				plasmidfinder    = plf_ch     // channel: [ val(meta), path(plasmidfinder) ]
				mlst             = mlst_ch    // channel: [ val(meta), path(mlst) ]
				report_html      = report_ch  // channel: [ val(meta), path(html) ]
}
	



