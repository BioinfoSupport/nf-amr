#!/usr/bin/env nextflow

include { get_tool_args } from '../../../modules/local/functions.nf'



process AMRFINDERPLUS_RUN {
	  container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.22--hf69ffd2_0'
    memory '6 GB'
    cpus 4
    input:
		    tuple val(meta), path(fasta), val(args)
		    path('amrfinder_db')
    output:
		    tuple val(meta), path("amrfinderplus", type: 'dir')
    script:
		    """
		    mkdir amrfinderplus
			  amrfinder -n $fasta $args --mutation_all amrfinderplus/mutations.tsv --database amrfinder_db --threads $task.cpus > amrfinderplus/report.tsv
		    """
}



process AMRFINDERPLUS_UPDATE {
    container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.22--hf69ffd2_0'
    memory '4 GB'
    cpus 1
    output:
		    path('amrfinder_db', type: 'dir')
    script:
    """
    amrfinder_update -d db && mv \$(realpath db/latest) amrfinder_db
    """
}

workflow AMRFINDERPLUS {
		take:
	    	fa_ch  // [meta,fasta]
		main:
			amrfinderplus_db = AMRFINDERPLUS_UPDATE()
			amrfinderplus_ch = AMRFINDERPLUS_RUN(
					fa_ch.map({meta,fasta -> [meta,fasta,get_tool_args('amrfinderplus',meta)]}),
					amrfinderplus_db
			)
			amrfinderplus_ch = Channel.empty()
		emit:
		    db = amrfinderplus_db
				results = amrfinderplus_ch  // channel: [meta,amrfinderplus]
}


