process MEDAKA_CONSENSUS {
    container 'docker.io/ontresearch/medaka:shac4e11bfa4e65668b28739ba32edc3af12baf7574-amd64'
    memory '20 GB'
    cpus 8
    time '4h'
    input:
        tuple val(meta), path('assembly.fasta'), path('long_reads.fastq.gz')
    output:
        tuple val(meta), path('medaka',type:'dir')
    script:
		    """
		    medaka_consensus \\
		      ${task.ext.args?:'--bacteria'} \\
			    -f \\
			    -d assembly.fasta \\
			    -i long_reads.fastq.gz \\
			    -o medaka
		    """
}
