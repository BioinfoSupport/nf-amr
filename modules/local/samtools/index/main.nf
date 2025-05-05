
process SAMTOOLS_INDEX {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    memory '2 GB'
    cpus 1
    input:
    		tuple val(meta), path(bam)
    output:
    		tuple val(meta), path ("*.bam.bai")
    script:
		    """
		    samtools index -@ ${task.cpus} ${task.ext.args?:''} bam
		    """
}
