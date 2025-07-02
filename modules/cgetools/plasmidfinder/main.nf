process PLASMIDFINDER {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    time '1h'
    input:
        tuple val(meta), path(seq)
    output:
				tuple val(meta), path("plasmidfinder/", type: 'dir')
    script:
		    """
	  		mkdir -p plasmidfinder
		  	plasmidfinder.py ${task.ext.args?:''} -q -x -p /db/plasmidfinder_db -i ${seq.join(' ')} -o 'plasmidfinder/'
		    """
}


