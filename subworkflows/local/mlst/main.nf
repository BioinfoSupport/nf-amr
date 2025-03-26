#!/usr/bin/env nextflow

params.skip_mlst = false
params.mlst_default_args = null // do not run by default

process RUN_MLST {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna), val(args)
    output:
    		tuple val(meta), path("mlst/", type: 'dir')
    script:
		    """
				mkdir -p mlst/ tmp/
				mlst.py -p ${task.ext.args?:''} '/db/mlst_db' -i '${assembly_fna}' -o 'mlst/' ${args} --tmp_dir 'tmp/' -x
		    """
}


workflow MLST {
		take:
	    	fa_ch  // [meta,meta_org,fasta]
		main:
			if (params.skip_mlst) {
				out_ch = Channel.empty()
			} else {
				out_ch = fa_ch
					.map({meta,meta_org,fasta -> 
							args = meta.mlst_args
							if ((args==null) && (meta_org.org_name) && (meta_org.org_name in params.organisms)) {
								args = params.organisms[meta_org.org_name].mlst_args
							}
							if (args==null) {
								args = params.mlst_default_args
							}					
							[meta,fasta,args]
					})
					.filter({meta,fasta,args -> args!=null})
					| RUN_MLST
			}
		emit:
				out_ch  // [meta,mlst]
}

