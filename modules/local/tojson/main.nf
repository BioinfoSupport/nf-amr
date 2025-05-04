
process TO_JSON {
	input:
		tuple(val(meta),val(json_content))
	output:
		tuple(val(meta),path("meta.json"))
	script:
		def builder = new groovy.json.JsonBuilder(json_content)
"""
cat >> meta.json << __EOF_META_JSON__
${builder.toPrettyString()}
__EOF_META_JSON__
"""
}

