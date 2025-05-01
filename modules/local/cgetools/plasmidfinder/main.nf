

process PLASMIDFINDER_RUN {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assembly_fna), val(args)
    output:
				tuple val(meta), path("plasmidfinder/", type: 'dir')
    script:
		    """
	  		mkdir -p plasmidfinder
		  	plasmidfinder.py ${task.ext.args?:''} -q -x -p /db/plasmidfinder_db ${args} -i '${assembly_fna}' -o 'plasmidfinder/'
		    """    
}


