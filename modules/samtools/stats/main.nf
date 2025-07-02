
process SAMTOOLS_STATS {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '2 GB'
    cpus 1
    time '15 min'
    input:
    		tuple val(meta), path("input.cram")
    output:
    		tuple val(meta), path("input.cram.stats")
    script:
				"""
				samtools stats -@ ${task.cpus} ${task.ext.args?:''} input.cram > input.cram.stats
				"""
}
