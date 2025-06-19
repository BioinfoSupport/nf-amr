
process NANOPLOT {
    container 'quay.io/biocontainers/nanoplot:1.41.0--pyhdfd78af_0'
    cpus = 2
    memory = '2 GB'
    input:
    	tuple val(meta), path("reads.fastq.gz")
    output:
    	tuple val(meta), path("nanoplot",type:'dir')
    script:
    """
    NanoPlot --threads ${task.cpus} --fastq reads.fastq.gz -o nanoplot
    """
}
