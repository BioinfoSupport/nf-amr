
process SAMTOOLS_FASTQ {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '2 GB'
    cpus 1
    time '15 min'
    input:
    		tuple val(meta), path("input.bam")
    output:
    		tuple val(meta), path("output.fastq.gz")
    script:
				"""
				samtools fastq -@ ${task.cpus} ${task.ext.args?:'-n -T "*" -0 -'} input.bam \\
				| bgzip -@ ${task.cpus} \\
				> output.fastq.gz
				"""
}


