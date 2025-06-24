process FASTQC {
    container 'biocontainers/fastqc:0.12.1--hdfd78af_0'
    cpus 3
    memory '5 GB'
    input:
	    tuple val(meta), path(reads)
    output:
	    tuple val(meta), path("report.html"), emit: html
    script:
	    """
	    fastqc \\
	        ${task.ext.args?:''} \\
	        --threads ${task.cpus} \\
	        --memory 5000 \\
	        ${reads}
	    """
}