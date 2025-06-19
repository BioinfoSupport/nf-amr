process MULTIQC {
    container 'quay.io/biocontainers/multiqc:1.28--pyhdfd78af_0'
    cpus = 2
    memory = '2 GB'
    input:
	    path('?/*')
	    path(multiqc_config)
    output:
	    path 'multiqc.html', emit: html
    script:
	    def config = multiqc_config ? "" : ''
	    """
	    multiqc ${task.ext.args?:''} --force --config ${multiqc_config} --filename multiqc.html ./
	    """
}