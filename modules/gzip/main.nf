
process GZIP_DECOMPRESS {
    container 'ubuntu:24.04'
    memory '2 GB'
    cpus 1
    time '30 min'
    input:
    		tuple val(meta), path(gzfile)
    output:
    		tuple val(meta), path(gzfile.baseName)
    script:
				"""
				gzip -d ${gzfile}
				"""
}


