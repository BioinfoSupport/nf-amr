#!/usr/bin/env nextflow

include { RESFINDER     } from '../resfinder'
include { PLASMIDFINDER } from '../plasmidfinder'
include { MLST          } from '../mlst'
include { ORG_MAP       } from '../org/map'
include { ORG_DB        } from '../org/db'


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

process AGGREGATE_ISOLATE_DATA {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
    memory '8 GB'
    cpus 4
    input:
    		tuple(val(meta),path("meta.json"),path("assembly.fasta"),path("ani.tsv"),path("resfinder"),path("mlst"),path("plasmidfinder"))
    output:
        tuple(val(meta),path("isolate_data.rds"))
    script:
				"""
				aggregate_isolate_data
				"""
}

process ISOLATE_HTML_REPORT {
	  container "registry.gitlab.unige.ch/amr-genomics/rscript:main"
    memory '8 GB'
    cpus 4
    input:
    		tuple(val(meta),path(rds_data))
    output:
        tuple(val(meta),path("isolate_report.html"))
    script:
				"""
				#!/usr/bin/env Rscript
				rmarkdown::render(
					'${moduleDir}/resources/isolate_report.Rmd',
					params = list(rds = '${rds_data}'),
					output_file = "./isolate_report.html"
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
							[meta,meta_json,fa,ani,res,mlst,plf]
					})
					| AGGREGATE_ISOLATE_DATA

				ISOLATE_HTML_REPORT(isolate_ch)

		emit:
				resfinder     = res_ch     // channel: [ val(meta), path(resfinder) ]
        org_ani       = org_ch     // channel: [ val(meta), val(org_name) ]
        org_db        = ORG_DB.out // channel: path(org_db) ]
				plasmidfinder = plf_ch     // channel: [ val(meta), path(plasmidfinder) ]
				mlst          = mlst_ch    // channel: [ val(meta), path(mlst) ]
				report_rds    = AGGREGATE_ISOLATE_DATA.out // channel: [val(meta), path(rds)]
				report_html   = ISOLATE_HTML_REPORT.out    // channel: [val(meta), path(html)]
}
	



