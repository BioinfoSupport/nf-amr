
process BWA_MEM {
    container 'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4'
    memory '5 GB'
    cpus 4    
    input:
	    tuple val(meta), path('index'), path(reads)
    output:
	    tuple val(meta), path("out.cram"), emit: cram
	    tuple val(meta), path("out.cram.crai"), emit: crai
    script:
	    """
	    bwa mem ${task.ext.args?:''} -t ${task.cpus} index/index ${reads.join(' ')} | samtools sort -@ ${task.cpus} --write-index -O CRAM -o out.cram
	    """
}
