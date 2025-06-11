
process COLLECT_FILES {
	input:
		val(all_files)
		val(outdir)
	output:
		path("${outdir}")
	script:
	"""
  mkdir -p ${outdir}
  ${ 
  	all_files
  		.findAll({it.size()==2})
	  	.collect({
	  			meta, file -> "mkdir -p '${outdir}/${meta.id}' && ln -s '${file}' '${outdir}/${meta.id}/'"
	  	})
	  	.join('\n')
  }
  ${ 
  	all_files
  		.findAll({it.size()==3})
	  	.collect({
	  			meta, file, target -> "mkdir -p '${outdir}/${meta.id}' && ln -s '${file}' '${outdir}/${meta.id}/${target}'"
	  	})
	  	.join('\n')
  }
  
	"""
}
