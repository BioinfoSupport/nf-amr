
process ORG_DB {
		container 'registry.gitlab.unige.ch/amr-genomics/species_profiler:main'
		memory '2 GB'
		cpus 1
		output:
		    path("org_db", type: 'dir')
		script:
				"""
				cp -r /app/db ./org_db
				"""
}
