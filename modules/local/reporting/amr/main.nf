
process REPORTING_AMR {
	  //container "registry.gitlab.unige.ch/amr-genomics/fatools:main"
    memory '8 GB'
    cpus 4
    input:
    		tuple val(meta), path(fasta), path(resfinder), path(org_map), path(org_env), path(mlst), path(plasmidfinder)
    output:
        path("*.amr_report.json"), emit: report
    script:
		    def args = task.ext.args ?: ''
		    """
			    cat <<EOF > amr_report.html
			    ${fna_meta}
			    ${res_meta}
			    -------------
			    EOF
			    find ./ >> amr_report.html
		    """
}


