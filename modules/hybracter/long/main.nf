process HYBRACTER_LONG {
	  container "quay.io/gbouras13/hybracter:0.11.2"
    memory '16 GB'
    cpus 8
    time '4h'
    input:
        tuple val(meta), path('ont_reads.fastq.gz')
    output:
				tuple val(meta), path('hybracter_out/', type: 'dir')
    script:
		    """
				hybracter long-single ${task.ext.args?:''} -l ont_reads.fastq.gz --auto --threads ${task.cpus}
		    """
}

