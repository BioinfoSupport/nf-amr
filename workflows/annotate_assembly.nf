#!/usr/bin/env nextflow

include { ORGFINDER_DETECT    } from '../modules/local/orgfinder/detect'

include { AMRFINDERPLUS_UPDATE } from '../modules/local/amrfinderplus/update'
include { AMRFINDERPLUS_RUN } from '../modules/local/amrfinderplus/run'
include { PROKKA_RUN        } from '../modules/local/prokka'

include { RESFINDER_FA_RUN  } from '../modules/local/cgetools/resfinder'
include { PLASMIDFINDER_RUN } from '../modules/local/cgetools/plasmidfinder'
include { MLST_RUN          } from '../modules/local/cgetools/mlst'
include { MOBTYPER_RUN      } from '../modules/local/mobsuite/mobtyper'

include { TO_JSON           } from '../modules/local/tojson'
include { SAMTOOLS_FAIDX    } from '../modules/local/samtools/faidx'
//include { RMD_RENDER        } from '../modules/local/rmd/render'

params.skip_prokka = true

params.resfinder_default_args = ''
params.mobtyper_default_args = ''
params.amrfinderplus_default_args = ''
params.plasmidfinder_default_args = ''
params.mlst_default_args = null // do not run by default
params.prokka_default_args = '--kingdom Bacteria'


// Return a boolean of whether the tool should be skip or not
def skip_tool(tool_name) {
	def key = 'skip_' + tool_name
	if (params.containsKey(key)) return params[key]
	return false
}

// Get name of the organism to use for given sample
def org_name(meta) {
		def org_key = 'org_name'
		if (params.containsKey(org_key)) return params[org_key]
    if (meta.containsKey(org_key)) return meta[org_key]
		return null
}

def get_key(meta,tool_name,org_name=null) {
		def default_args_key = tool_name + "_default_args"
		def key = tool_name + "_args"
		if (params.containsKey(key)) return params[key]
    if (meta.containsKey(key)) return meta[key]
    
    if (org_name==null) return params[default_args_key]
    def org_args = params.organisms.containsKey(org_name)?params.organisms[org_name] : [:]
    if (org_args.containsKey(key)) return org_args[key]
    return params[default_args_key]
}

// Retreive arguments to use for the given tool on the given sample
def tool_args(tool_name, meta, org_name=null) {
		return get_key(meta,tool_name, org_name)
}


workflow ANNOTATE_ASSEMBLY {
		take:
	    	fa_ch    // channel: [ val(meta), path(assembly_fna) ]
		main:
	
	      // ---------------------------------------------------------------------
	      // Tools that can run directly on a FASTA witout specifying an organism
	      // ---------------------------------------------------------------------
				SAMTOOLS_FAIDX(fa_ch)
	      
				// CGE - RESFINDER
				if (skip_tool('resfinder')) {
						resfinder_ch = Channel.empty()
				} else {
						resfinder_ch = fa_ch
							.map({meta,fasta -> [meta,fasta,tool_args('resfinder',meta)]})
							.filter({meta,fasta,args -> args!=null})
			        | RESFINDER_FA_RUN
				}
				
				// NCBI AMRfinder+
				amrfinderplus_db = AMRFINDERPLUS_UPDATE()
				if (skip_tool('amrfinderplus')) {
						amrfinderplus_ch = Channel.empty()
				} else {
						amrfinderplus_ch = AMRFINDERPLUS_RUN(
								fa_ch
								  .map({meta,fasta -> [meta,fasta,tool_args('amrfinderplus',meta)]})
									.filter({meta,fasta,args -> args!=null}),
								amrfinderplus_db
						)
				}
				
				// MOBsuite - MOBtyper
				if (skip_tool('mobtyper')) {
						mobtyper_ch = Channel.empty()
				} else {
						mobtyper_ch = fa_ch
							.map({meta,fasta -> [meta,fasta,tool_args('mobtyper',meta)]})
							.filter({meta,fasta,args -> args!=null})
			        | MOBTYPER_RUN
				}

		
	      // ----------------------------------------------------
	      // Tools specific to an organism
	      // ----------------------------------------------------
				//detected_org_ch = ORG_DETECT(fa_ch.filter({meta,fa -> org_name(meta)==null}),ORG_DB.out)
				detected_org_ch = ORGFINDER_DETECT(fa_ch)
				
				// Update fa_ch to add detected organism
				fa_org_ch = fa_ch
					.join(detected_org_ch.org_name,remainder:true)
					.map({meta,fa,detected_org_name -> [meta,fa,org_name(meta)?:detected_org_name]})
				
				// Plasmid typing
				if (skip_tool('plasmidfinder')) {
					plf_ch = Channel.empty()
				} else {
						plf_ch = fa_org_ch
							.map({meta,fa,org_name -> [meta, fa, tool_args('plasmidfinder',meta, org_name)]})
							.filter({meta,fasta,args -> args!=null})
							| PLASMIDFINDER_RUN
				}
				
				// MLST typing
				if (skip_tool('mlst')) {
						mlst_ch = Channel.empty()
				} else {
						mlst_ch = fa_org_ch
							.map({meta,fa,org_name -> [meta, fa, tool_args('mlst',meta,org_name)]})
							.filter({meta,fasta,args -> args!=null})
							| MLST_RUN
				}
				
				// PROKKA annotations
				if (skip_tool('prokka')) {
						prokka_ch = Channel.empty()
				} else {
						prokka_ch = fa_org_ch
						  .map({meta,fa,org_name -> [meta, fa, tool_args('prokka',meta,org_name)]})
						  .filter({meta,fasta,args -> args!=null})
							| PROKKA_RUN
				}

				runinfo_json_ch = fa_org_ch
					.join(detected_org_ch.org_name,remainder:true)
					.join(detected_org_ch.org_ani,remainder:true)
					.join(detected_org_ch.org_acc,remainder:true)
					.map({meta,fa,org_name,detected_org_name,detected_org_ani,detected_org_acc -> [
						meta, 
						[runinfo: [
							meta: meta, 
							org_name: org_name,
							org_detection: [org_name: detected_org_name, org_ani: detected_org_ani, org_acc: detected_org_acc]
						]]]})
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
		    runinfo          = runinfo_json_ch    // channel: [ val(meta), path(resfinder) ]
		    faidx            = SAMTOOLS_FAIDX.out // channel: [ val(meta), path(fai) ]
		    prokka           = prokka_ch          // channel: [ val(meta), path(prokka) ]
		    amrfinderplus_db = amrfinderplus_db   // channel: path(amrfinderplus_db) ]
		    amrfinderplus    = amrfinderplus_ch   // channel: val(meta), path(amrfinderplus) ]
				resfinder        = resfinder_ch       // channel: [ val(meta), path(resfinder) ]
				mobtyper         = mobtyper_ch        // channel: [ val(meta), path(mobtyper) ]
				plasmidfinder    = plf_ch             // channel: [ val(meta), path(plasmidfinder) ]
				mlst             = mlst_ch            // channel: [ val(meta), path(mlst) ]
				report_html      = report_ch          // channel: [ val(meta), path(html) ]
				orgfinder        = detected_org_ch.orgfinder // channel: [ val(meta), val(org_name) ]
}
	



