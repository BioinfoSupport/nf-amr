
process SAMTOOLS_INDEX {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '4 GB'
    cpus 2
    time '15 min'
    input:
    		tuple val(meta), path("align.bam")
    output:
    		tuple val(meta), path("align.bam.bai")
    script:
				"""
				samtools index -@ ${task.cpus} ${task.ext.args?:''} align.bam
				"""
}
