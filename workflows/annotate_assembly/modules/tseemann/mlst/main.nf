process MLST {
	  container "docker.io/staphb/mlst:2.23.0-2025-02-01"
    memory '4 GB'
    cpus 1
    time '30 min'
    input:
        tuple val(meta), path("assembly.fasta"), val(args)
    output:
    		tuple val(meta), path("MLST.txt")
    script:
		    """
		    mlst --threads=${task.cpus} ${task.ext.args?:''} ${args} assembly.fasta > MLST.txt
		    """
}

