process FLYE {
    container 'quay.io/biocontainers/flye:2.9.6--py311h2de2dd3_0'
    memory '20 GB'
    cpus 8
    time '4h'
    input:
        tuple val(meta), path('long_reads.fastq.gz')
    output:
        tuple val(meta), path('flye',type:'dir')
    script:
		    """
		    flye \\
		      ${task.ext.args?:''} \\
		      --threads ${task.cpus} \\
		      --out-dir ./flye \\
		      --nano-hq long_reads.fastq.gz
		    """
}
