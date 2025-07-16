process BWA_MEM {
    container 'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4'
    memory '5 GB'
    cpus 4
    time '2h'
    input:
	    tuple val(meta), path('index'), path(reads)
    output:
	    tuple val(meta), path("out.bam"), emit: bam
	    tuple val(meta), path("out.bam.bai"), emit: bai
    script:
	    """
	    bwa mem ${task.ext.args?:''} -t ${task.cpus} index/index ${reads.join(' ')} | samtools sort -@ ${task.cpus} --write-index -O BAM -o out.bam##idx##out.bam.bai
	    """
}
