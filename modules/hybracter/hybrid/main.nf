#!/usr/bin/env nextflow


process HYBRACTER_HYBRID {
	  container "quay.io/gbouras13/hybracter:0.11.2"
    memory '12 GB'
    cpus 8
    input:
        tuple val(meta), path('ont_reads.fastq.gz'), path('R1.fastq.gz'), path('R2.fastq.gz')
    output:
				tuple val(meta), path("hybracter_out/FINAL_OUTPUT/complete", type: 'dir')
    script:
		    """
				hybracter hybrid-single ${task.ext.args?:''} -l ont_reads.fastq.gz -1 R1.fastq.gz -2 R2.fastq.gz --auto --threads ${task.cpus}
		    """
}

