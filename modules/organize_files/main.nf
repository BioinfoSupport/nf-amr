
process ORGANIZE_FILES {
    input:
    	tuple(val(meta),val(file_pairs)) // liste de [[path(src), target]]
    output:
    	tuple(val(meta),path("output",type:'dir'))
    script:
    	def cmd = file_pairs
    		.findAll({src,target -> src})
    		.collect({src,target ->
			    def dest_path = "output/${target}"
			    return "mkdir -p \$(dirname ${dest_path}) && ln -s ${src} ${dest_path}".stripIndent()
    	})
    """
    mkdir -p output
    ${cmd.join('\n')}
    """
}



