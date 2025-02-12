
process ORGANISM_DB {
		container 'registry.gitlab.unige.ch/amr-genomics/species_profiler:main'
		memory '2 GB'
		cpus 1
		output:
		    path("org_db.tsv")
		script:
				"""
				cp /app/db/db.tsv org_db.tsv
				"""
}
