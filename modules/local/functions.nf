


def get_tool_args(tool_name, meta, orgs_args=null, default_args=null, unknown_value="") {
		if (orgs_args==null) orgs_args = [:]
		if (default_args==null) default_args = [:]
    def tool_key = "${tool_name}_args"
    if (meta.containsKey(tool_key)) return meta[tool_key]
    if (meta.containsKey('organism') && orgs_args.containsKey(meta.organism) && orgs_args[meta.organism].containsKey(tool_key)) {
    		default_args = orgs_args[meta.organism][tool_key]
    } else if (default_args.containsKey(tool_key)) {
    		return default_args[tool_key]
    } else {
    		return unknown_value
    }
}



