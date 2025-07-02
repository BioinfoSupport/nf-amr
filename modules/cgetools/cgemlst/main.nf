process CGEMLST_RUN {
	  container "registry.gitlab.unige.ch/amr-genomics/cgetools:main"
    memory '4 GB'
    cpus 1
    time '1h'
    input:
        tuple val(meta), path(assembly_fna), val(args)
    output:
    		tuple val(meta), path("cge_mlst/", type: 'dir')
    script:
		    """
				mkdir -p cge_mlst/ tmp/
				mlst.py -p ${task.ext.args?:''} '/db/mlst_db' -i '${assembly_fna}' -o 'cge_mlst/' ${args} --tmp_dir 'tmp/' -x
		    """
}

