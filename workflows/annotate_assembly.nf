#!/usr/bin/env nextflow

include { ORGFINDER_DETECT    } from '../modules/local/orgfinder/detect'

include { AMRFINDERPLUS_UPDATE } from '../modules/local/amrfinderplus/update'
include { AMRFINDERPLUS_RUN } from '../modules/local/amrfinderplus/run'
include { PROKKA_RUN        } from '../modules/local/tseemann/prokka'
include { MLST_RUN          } from '../modules/local/tseemann/mlst'

include { RESFINDER         } from '../modules/local/cgetools/resfinder'
include { PLASMIDFINDER     } from '../modules/local/cgetools/plasmidfinder'
include { CGEMLST_RUN       } from '../modules/local/cgetools/cgemlst'
include { MOBTYPER_RUN      } from '../modules/local/mobsuite/mobtyper'

include { TO_JSON           } from '../modules/local/tojson'
include { SAMTOOLS_FAIDX    } from '../modules/local/samtools/faidx'

params.skip_prokka = true
params.mobtyper_default_args = ''
params.amrfinderplus_default_args = ''
params.plasmidfinder_default_args = ''
params.cgemlst_default_args = null // do not run by default
params.MLST_default_args = ''      // autodetect species by default
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
				fai_ch = SAMTOOLS_FAIDX(fa_ch)
	      
				// CGE - RESFINDER
				if (skip_tool('resfinder')) {
						resfinder_ch = Channel.empty()
				} else {
						resfinder_ch = RESFINDER(fa_ch,"fasta")
				}

				// Plasmid typing
				if (skip_tool('plasmidfinder')) {
						plasmidfinder_ch = Channel.empty()
				} else {
						plasmidfinder_ch = PLASMIDFINDER(fa_ch)
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
				//orgfinder_ch = ORG_DETECT(fa_ch.filter({meta,fa -> org_name(meta)==null}),ORG_DB.out)
				orgfinder_ch = ORGFINDER_DETECT(fa_ch)
				
				// Update fa_ch to add detected organism
				fa_org_ch = fa_ch
					.join(orgfinder_ch.org_name,remainder:true)
					.map({meta,fa,detected_org_name -> [meta,fa,org_name(meta)?:detected_org_name]})
				
				// MLST typing
				if (skip_tool('cgemlst')) {
						cgemlst_ch = Channel.empty()
				} else {
						cgemlst_ch = fa_org_ch
							.map({meta,fa,org_name -> [meta, fa, tool_args('cgemlst',meta,org_name)]})
							.filter({meta,fasta,args -> args!=null})
							| CGEMLST_RUN
				}
				
				if (skip_tool('MLST')) {
						MLST_ch = Channel.empty()
				} else {
						MLST_ch = fa_org_ch
						  .map({meta,fa,org_name -> [meta, fa, tool_args('MLST',meta,org_name)]})
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

				runinfo_ch = fa_org_ch
					.join(orgfinder_ch.org_name,remainder:true)
					.join(orgfinder_ch.ani,remainder:true)
					.join(orgfinder_ch.org_acc,remainder:true)
					.map({meta,fa,org_name,orgfinder_org_name,orgfinder_ani,orgfinder_org_acc -> [
						meta, 
						[runinfo: [
							meta: meta, 
							org_name: org_name,
							orgfinder: [org_name: orgfinder_org_name, ani: orgfinder_ani, org_acc: orgfinder_org_acc]
						]]]})
				runinfo_ch = TO_JSON(runinfo_ch)

		emit:
				fai = fai_ch
				runinfo = runinfo_ch
				orgfinder = orgfinder_ch.orgfinder
				amrfinderplus = amrfinderplus_ch
				resfinder = resfinder_ch
				mobtyper = mobtyper_ch
				plasmidfinder = plasmidfinder_ch
				cgemlst = cgemlst_ch
				MLST = MLST_ch
				prokka = prokka_ch
}
	



