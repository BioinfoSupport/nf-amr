
process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0'
    cpus = 2
    memory = '2 GB'
    time '1 h'
    input:
	    path('db')
	    path('config.yml')
    output:
	    path 'multiqc.html', emit: html
    script:
	    """
	    multiqc ${task.ext.args?:''} --force --config config.yml --filename multiqc.html db/
	    """
}

