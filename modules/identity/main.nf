
process IDENTITY {
	input:
		tuple val(meta),path(x)
	output:
		tuple val(meta),path(x)
	script:
		"ls"
}

