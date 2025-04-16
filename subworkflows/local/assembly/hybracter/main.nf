#!/usr/bin/env nextflow

params.hybracter_default_args = ''

process RUN_HYBRACTER {
	  container "quay.io/gbouras13/hybracter:0.11.2"
    memory '12 GB'
    cpus 4
    input:
        tuple val(meta), path('ont_reads.fastq.gz'), path('R1.fastq.gz'), path('R2.fastq.gz'), val(args)
    output:
				tuple val(meta), path("hybracter/", type: 'dir')
    script:
		    """
				hybracter ${task.ext.args?:''} --outdir 'hybracter/' --cpus ${task.cpus} ont_reads.fastq.gz
		    """    
}

workflow HYBRACTER {
		take:
	    	ont_fq_ch
		main:
			out_ch = RUN_HYBRACTER(ont_fq_ch)
		emit:
				out_ch
}

