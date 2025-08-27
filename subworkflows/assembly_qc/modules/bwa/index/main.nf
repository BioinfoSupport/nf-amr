process BWA_INDEX {
    container 'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4'
    memory '5 GB'
    cpus 3
    time '1h'
    input:
	    tuple val(meta), path('ref.fasta')
    output:
	    tuple val(meta), path('index',type:'dir'), emit: index
    script:
	    """
	    mkdir -p "index" && bwa index ${task.ext.args?:''} -p "index/index" ref.fasta
	    """
}
