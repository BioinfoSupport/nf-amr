#!/usr/bin/env nextflow

include {ORGANISM_MAP} from './modules/local/org_map'
include {CGE_MLST} from './modules/local/cgetools/mlst'




workflow {
	fa_ch = Channel.fromPath(params.input)
			.map({x -> tuple(["id":x.baseName],x)})
	fa_ch = ORGANISM_MAP(fa_ch)
		.map({meta,assembly_fna,org_map,org_env,org_name -> 
			meta.org_name = org_name
			[meta,assembly_fna]
		})
	
	
	fa_ch.filter({meta,x -> 
		meta.containsKey("org_name")
		&& params.organisms.containsKey(meta.org_name) 
		&& params.organisms[meta.org_name]
	})
	| CGE_MLST
}
	



