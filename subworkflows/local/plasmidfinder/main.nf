#!/usr/bin/env nextflow

params.skip_plasmidfinder = false
params.plasmidfinder_default_args = ''

process RUN_PLASMIDFINDER {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna), val(args)
    output:
				tuple val(meta), path("plasmidfinder/", type: 'dir')
    script:
		    """
	  		mkdir -p plasmidfinder
		  	plasmidfinder.py ${task.ext.args?:''} -q -x -p /db/plasmidfinder_db ${args} -i '${assembly_fna}' -o 'plasmidfinder/'
		    """    
}


workflow PLASMIDFINDER {
		take:
	    	fa_ch // [meta,meta_org,fasta]
		main:
			if (params.skip_plasmidfinder) {
				out_ch = Channel.empty()
			} else {
				out_ch = fa_ch
					.map({meta,meta_org,fasta -> 
							args = meta.plasmidfinder_args
							args = args!=null?args:params.organisms[meta_org.org_name].plasmidfinder_args
							args = args!=null?args:params.plasmidfinder_default_args
							[meta,fasta,args]
					})
					.filter({meta,fasta,args -> args!=null})
					| RUN_PLASMIDFINDER
			}
		emit:
				out_ch // [meta,plasmidfinder]
}

