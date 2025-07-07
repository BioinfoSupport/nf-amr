process CHOPPER {
    container 'quay.io/biocontainers/chopper:0.10.0--hcdda2d0_0'
    cpus 4
    memory '4 GB'
    time '30 min'
    input:
	    tuple val(meta), path("reads.fastq.gz")
    output:
	    tuple val(meta), path("filtered_reads.fastq.gz"), emit: fastq
    script:
    """
    zcat $fastq \\
    | chopper --threads ${task.cpus} ${task.ext.args?:'-q 15 -l 500'} \\
    | gzip \\
    > filtered_reads.fastq.gz
    """
}

