#!/usr/bin/env nextflow

params.skip_mobtyper = false
params.mobtyper_default_args = ''

process RUN_MOBTYPER {
	  container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna), val(args)
    output:
				tuple val(meta), path('mobtyper.tsv')
    script:
		    """
		    mob_typer ${task.ext.args?:''} ${args} --multi -i '${assembly_fna}' -o mobtyper.tsv
		    """    
}


workflow MOBTYPER {
		take:
	    	fa_ch // [meta,fasta]
		main:
			if (params.skip_mobtyper) {
				out_ch = Channel.empty()
			} else {
				out_ch = fa_ch
					.map({meta,fasta -> 
							args = meta.mobtyper_args
							if (args==null) {
								args = params.mobtyper_default_args
							}
							[meta,fasta,args]
					})
					.filter({meta,fasta,args -> args!=null})
					| RUN_MOBTYPER
			}
		emit:
				out_ch // [meta,plasmidfinder]
}
