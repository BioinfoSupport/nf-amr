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

process ISOLATE_AGGREGATE_DATA {
	  container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '8 GB'
    cpus 4
    input:
    		tuple(val(meta),path("meta.json"),path("assembly.fasta"),path("ani.tsv"),path("resfinder"),path("mlst"),path("plasmidfinder"))
    output:
        tuple(val(meta),path("isolate_data.rds"))
    script:
				"""
				isolate_aggregate_data
				"""
}

process ISOLATE_TEXT_REPORT {
	  container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '8 GB'
    cpus 4
    input:
    		tuple(val(meta),path(rds_data))
    output:
        tuple(val(meta),path("isolate_report.txt"))
    script:
				"""
				isolate_text_report
				"""
}

process ISOLATE_HTML_REPORT {
	  container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '8 GB'
    cpus 4
    input:
    		tuple(val(meta),path(rds_data))
    output:
        tuple(val(meta),path("isolate_report.html"))
    script:
				"""
				isolate_html_report
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
					| ISOLATE_AGGREGATE_DATA

				ISOLATE_TEXT_REPORT(isolate_ch)
				ISOLATE_HTML_REPORT(isolate_ch)

		emit:
				resfinder     = res_ch     // channel: [ val(meta), path(resfinder) ]
        org_ani       = org_ch     // channel: [ val(meta), val(org_name) ]
        org_db        = ORG_DB.out // channel: path(org_db) ]
				plasmidfinder = plf_ch     // channel: [ val(meta), path(plasmidfinder) ]
				mlst          = mlst_ch    // channel: [ val(meta), path(mlst) ]
				report_rds    = ISOLATE_AGGREGATE_DATA.out // channel: [val(meta), path(rds)]
				report_txt    = ISOLATE_TEXT_REPORT.out    // channel: [val(meta), path(json)]
				report_html   = ISOLATE_HTML_REPORT.out    // channel: [val(meta), path(html)]
}
	



