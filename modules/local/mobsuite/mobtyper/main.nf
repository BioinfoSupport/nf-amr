

process MOBTYPER_RUN {
	  container "quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_0"
	  errorStrategy 'ignore'
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path('assembly.fna'), val(args)
    output:
				tuple val(meta), path('mobtyper.tsv')
    script:
		    """
		    mob_typer ${task.ext.args?:''} ${args} --multi -i 'assembly.fna' -o mobtyper.tsv
		    """    
}
