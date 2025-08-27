process PILON {
    container 'quay.io/biocontainers/pilon:1.24--hdfd78af_0'
    memory '10 GB'
    cpus 6
    time '4h'
    input:
        tuple val(meta), path('assembly.fasta'), path('short_reads.bam'), path('short_reads.bam.bai')
    output:
        tuple val(meta), path('pilon',type:'dir')
    script:
    """
    mkdir pilon \
    && java -Xmx8000M -jar /usr/local/share/pilon*/pilon.jar \
       --genome assembly.fasta \
       --bam short_reads.bam \
       ${task.ext.args?:'--fix all --changes'} \
       --threads ${task.cpus} \
       --output pilon/pilon
    """
}

