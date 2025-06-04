
process COLLECT_FILES {
	input:
		val all_files
	output:
		path('results')
	script:
	"""
  mkdir -p results
  ${ 
  	all_files
  		.findAll({it.size()==2})
	  	.collect({
	  			meta, file -> "mkdir -p 'results/${meta.id}' && ln -s '${file}' 'results/${meta.id}/'"
	  	})
	  	.join('\n')
  }
  ${ 
  	all_files
  		.findAll({it.size()==3})
	  	.collect({
	  			meta, file, target -> "mkdir -p 'results/${meta.id}' && ln -s '${file}' 'results/${meta.id}/${target}'"
	  	})
	  	.join('\n')
  }
  
	"""
}
