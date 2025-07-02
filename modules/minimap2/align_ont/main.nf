process MINIMAP2_ALIGN_ONT {
    container 'community.wave.seqera.io/library/minimap2_samtools:33bb43c18d22e29c'
    memory '10 GB'
    cpus 4
    time '1h'
    input:
	    tuple val(meta), path('ref.fasta'), path('reads.fastq.gz')
    output:
	    tuple val(meta), path("out.cram"), emit: cram
	    tuple val(meta), path("out.cram.crai"), emit: crai
    script:
	    """
	    minimap2 -x map-ont ${task.ext.args?:''} -t ${task.cpus} ref.fasta reads.fastq.gz -a | samtools sort -@ ${task.cpus} --write-index -O CRAM -o out.cram
	    """
}

