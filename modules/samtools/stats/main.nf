
process SAMTOOLS_STATS {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '12 GB'
    cpus 1
    time '30 min'
    input:
    		tuple val(meta), path("input.bam")
    output:
    		tuple val(meta), path("input.bam.stats")
    script:
				"""
				samtools stats -@ ${task.cpus} ${task.ext.args?:''} input.bam > input.bam.stats
				"""
}
