#!/usr/bin/env nextflow

params.skip_prokka = false
params.prokka_default_args = "--kingdom Bacteria"

process RUN_PROKKA {
	  container "docker.io/staphb/prokka:1.14.6"
    memory '12 GB'
    cpus 4
    input:
        tuple val(meta), path('assembly.fasta'), val(args)
    output:
				tuple val(meta), path("prokka/", type: 'dir')
    script:
		    """
				prokka ${task.ext.args?:''} --outdir 'prokka/' --cpus ${task.cpus} --prefix prokka ${args} assembly.fasta
		    """    
}

workflow PROKKA {
		take:
	    	fa_ch
		main:
			if (params.skip_prokka) {
				out_ch = Channel.empty()
			} else {
				out_ch = fa_ch
					.map({meta,meta_org,fasta -> 
							args = meta.prokka_args
							if ((args==null) && (meta_org.org_name) && (meta_org.org_name in params.organisms)) {
								args = params.organisms[meta_org.org_name].prokka_args
							}
							if (args==null) {
								args = params.prokka_default_args
							}
							[meta,fasta,args]
					})
					.filter({meta,fasta,args -> args!=null})
					| RUN_PROKKA
			}
		emit:
				out_ch
}

