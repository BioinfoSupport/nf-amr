

process MOBTYPER_RUN {
	  container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0"
    memory '10 GB'
    cpus 4
    time '30 min'
    input:
        tuple val(meta), path('assembly.fna'), val(args)
    output:
				tuple val(meta), path('mobtyper.tsv')
    script:
		    """
		    mob_typer ${task.ext.args?:''} ${args} --num_threads=${task.cpus} --multi -i 'assembly.fna' -o mobtyper.tsv
		    """    
}
