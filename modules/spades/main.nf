process SPADES {
    container 'quay.io/biocontainers/spades:4.2.0--h8d6e82b_1'
    memory '16 GB'
    cpus 8
    time '4h'
    input:
        tuple val(meta), path(illumina), path(nanopore)
    output:
        tuple val(meta), path('spades',type='dir')
    script:
	      def illumina_reads = illumina?( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ):""
	      def nanopore_reads = nanopore?"--nanopore $nanopore":""
		    """
		    mkdir -p spades && spades.py \\
		        ${task.ext.args?:''} \\
		        --threads ${task.cpus} \\
		        --memory ${task.memory.toGiga()} \\
		        $illumina_reads \\
		        $nanopore_reads \\
		        -o ./spades/
		    """
}