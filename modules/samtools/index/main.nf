
process SAMTOOLS_INDEX {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '4 GB'
    cpus 2
    input:
    		tuple val(meta), path("align.cram")
    output:
    		tuple val(meta), path("align.cram.crai")
    script:
				"""
				samtools index -@ ${task.cpus} ${task.ext.args?:''} align.cram
				"""
}
