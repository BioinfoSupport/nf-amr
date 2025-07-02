process HYBRACTER_HYBRID {
	  container "quay.io/gbouras13/hybracter:0.11.2"
    memory '16 GB'
    cpus 8
    time '4h'
    input:
        tuple val(meta), path('ont_reads.fastq.gz'), path('R1.fastq.gz'), path('R2.fastq.gz')
    output:
				tuple val(meta), path('hybracter_out/', type: 'dir')
    script:
		    """
				hybracter hybrid-single ${task.ext.args?:''} -l ont_reads.fastq.gz -1 R1.fastq.gz -2 R2.fastq.gz --auto --threads ${task.cpus}
		    """
}

