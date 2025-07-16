process UNICYCLER {
    container 'quay.io/biocontainers/unicycler:0.5.1--py312hdcc493e_4'
    memory '20 GB'
    cpus 8
    time '4h'
    input:
        tuple val(meta), path(illumina), path(nanopore)
    output:
        tuple val(meta), path('unicycler',type:'dir')
    script:
	      def nanopore_reads = nanopore?"-l $nanopore":''
	      illumina = illumina instanceof List?illumina:[illumina]
	      def short_reads = illumina?(illumina.size()==1?"-s ${illumina[0]}":"-1 ${illumina[0]} -2 ${illumina[1]}"):''
		    """
		    #${illumina.size()}
		    #${illumina.size()}
		    #${illumina.class}
		    #${illumina.join(" ")}
		    unicycler \\
		        ${task.ext.args?:''} \\
		        --threads ${task.cpus} \\
		        ${short_reads} ${nanopore_reads} \\
		        --out ./unicycler/
		    """
}
