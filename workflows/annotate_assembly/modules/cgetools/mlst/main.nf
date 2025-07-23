process MLST {
		label 'cgetools'
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    time '1h'
    input:
        tuple val(meta), path(assembly_fna), val(args)
    output:
    		tuple val(meta), path("mlst/", type: 'dir')
    script:
		    """
				mkdir -p mlst/ tmp/
				mlst.py -p ${task.ext.args?:''} '/db/mlst_db' -i '${assembly_fna}' -o 'mlst/' ${args} --tmp_dir 'tmp/' -x
		    """
}

