process PILON {
    container 'quay.io/biocontainers/pilon:1.24--hdfd78af_0'
    memory '20 GB'
    cpus 8
    time '4h'
    input:
        tuple val(meta), path('assembly.fasta'), path('short_reads.bam'), path('short_reads.bam.bai')
    output:
        tuple val(meta), path('pilon',type:'dir')
    script:
		    """
		    pilon ${task.ext.args?:'--fix all --changes --frags'} \\
		      --threads ${args.cpus} \\
		      --genome assembly.fasta \\
		      --bam short_reads.bam \\
		      --output pilon
		    """
}

