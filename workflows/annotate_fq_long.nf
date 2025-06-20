#!/usr/bin/env nextflow


include { RESFINDER_FA_RUN  } from '../modules/local/cgetools/resfinder'
include { PLASMIDFINDER_RUN } from '../modules/local/cgetools/plasmidfinder'
include { CGEMLST_RUN       } from '../modules/local/cgetools/cgemlst'


workflow ANNOTATE_FQ_LONG {
		take:
	    	fq_ch    // channel: [ val(meta), path(long_reads_fq) ]
		main:
	      
				// CGE - RESFINDER
				if (skip_tool('resfinder')) {
						resfinder_ch = Channel.empty()
				} else {
						resfinder_ch = fa_ch
							.map({meta,fasta -> [meta,fasta,tool_args('resfinder',meta)]})
							.filter({meta,fasta,args -> args!=null})
			        | RESFINDER_FA_RUN
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
				
				// Plasmid typing
				if (skip_tool('plasmidfinder')) {
						plasmidfinder_ch = Channel.empty()
				} else {
						plasmidfinder_ch = fa_org_ch
							.map({meta,fa,org_name -> [meta, fa, tool_args('plasmidfinder',meta, org_name)]})
							.filter({meta,fasta,args -> args!=null})
							| PLASMIDFINDER_RUN
				}
				
				// MLST typing
				if (skip_tool('cgemlst')) {
						cgemlst_ch = Channel.empty()
				} else {
						cgemlst_ch = fa_org_ch
							.map({meta,fa,org_name -> [meta, fa, tool_args('cgemlst',meta,org_name)]})
							.filter({meta,fasta,args -> args!=null})
							| CGEMLST_RUN
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
				runinfo = runinfo_ch
				orgfinder = orgfinder_ch.orgfinder
				resfinder = resfinder_ch
				plasmidfinder = plasmidfinder_ch
				cgemlst = cgemlst_ch
}
	



