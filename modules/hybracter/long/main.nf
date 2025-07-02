#!/usr/bin/env nextflow

process HYBRACTER_LONG {
	  container "quay.io/gbouras13/hybracter:0.11.2"
    memory '12 GB'
    cpus 8
    input:
        tuple val(meta), path('ont_reads.fastq.gz')
    output:
				tuple val(meta), path("hybracter_out/", type: 'dir')
    script:
		    """
				hybracter long-single ${task.ext.args?:''} -l ont_reads.fastq.gz --auto --threads ${task.cpus}
		    """
}

