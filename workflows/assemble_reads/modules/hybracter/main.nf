process HYBRACTER {
	  container "quay.io/gbouras13/hybracter:0.11.2"
    memory '20 GB'
    cpus 8
    time '4h'
    input:
        tuple val(meta), path(illumina), path(nanopore)
    output:
				tuple val(meta), path('hybracter/', type: 'dir')
    script:
    		def cmd = illumina?'hybrid-single':'long-single'
	      def nanopore_reads = nanopore?"-l $nanopore":''
	      illumina = illumina instanceof List?illumina:[illumina]
	      def short_reads = illumina?(illumina.size()==1?"-s ${illumina[0]}":"-1 ${illumina[0]} -2 ${illumina[1]}"):''
		    """
				hybracter ${cmd} \\
				    ${task.ext.args?:''} \\
				    --threads ${task.cpus} \\
				    --auto \\
				    ${short_reads} ${nanopore_reads} \\
				    -o hybracter
		    """
}


