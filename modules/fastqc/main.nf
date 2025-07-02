process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    cpus 3
    memory '5 GB'
    time '1h'
    input:
	    tuple val(meta), path(reads)
    output:
	    tuple val(meta), path("*_fastqc.html"), emit: html
	    tuple val(meta), path("*_fastqc.zip"), emit: zip
    script:
	    """
	    fastqc \\
	        ${task.ext.args?:''} \\
	        --threads ${task.cpus} \\
	        --memory 5000 \\
	        ${reads}
	    """
}