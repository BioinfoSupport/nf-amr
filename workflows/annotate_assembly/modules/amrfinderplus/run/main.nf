
process AMRFINDERPLUS_RUN {
	  container 'quay.io/biocontainers/ncbi-amrfinderplus:4.0.22--hf69ffd2_0'
    memory '6 GB'
    cpus 4
    time '30 min'
    input:
		    tuple val(meta), path('assembly.fasta'), val(args)
		    path('amrfinder_db')
    output:
		    tuple val(meta), path("amrfinderplus", type: 'dir')
    script:
		    """
		    mkdir amrfinderplus
			  amrfinder ${task.ext.args?:''} -n assembly.fasta ${args} --mutation_all amrfinderplus/mutations.tsv --database amrfinder_db --threads $task.cpus > amrfinderplus/report.tsv
		    """
}

