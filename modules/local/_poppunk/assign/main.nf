

process POPPUNK_RUN {
	  container "docker.io/staphb/mlst:2.23.0-2025-02-01"
    memory '4 GB'
    cpus 1
    input:
        tuple val(meta), path(assemblies.fasta), val(args)
    output:
    		tuple val(meta), path("poppunk",	type:'dir')
    script:
		    """
		    mkdir -p poppunk

# Create input file
cat <<EOF > poppunk/qfile.txt
MS1 ms1_assembled.fa
MS2 ms2_assembled.fa
SM14 SM14_1.fq.gz SM14_2.fq.gz
EOF

		    poppunk_assign --threads ${task.cpus} ${task.ext.args?:''} ${args} --db database --query poppunk/qfile.txt --output poppunk/clusters
		    """
}

