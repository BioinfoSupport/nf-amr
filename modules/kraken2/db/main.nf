process NCBI_DOWNLOAD_GENOMES {
		container 'docker.io/staphb/ncbi-datasets:18.0.2'
		memory '8 GB'
		cpus 2
		time '1h'
		input:
		    tuple val(taxon),val(args)
		output:
		    path('genomes/ncbi_dataset/data/*/*_genomic.fna')
		script:
		"""
		datasets download genome taxon '${taxon}' --reference ${task.ext.args} ${args} --filename archive.zip
		unzip archive.zip -d genomes
		"""
}

process KRAKEN2_DB_DOWNLOAD_TAX {
    container 'docker.io/staphb/kraken2:2.1.5'
    memory '10 GB'
    cpus 4
    time '1h'
    output:
    		path("db/taxonomy",type: 'dir')
    script:
		    """
				kraken2-build --download-taxonomy --db db
		    """
}

process KRAKEN2_DB_BUILD {
    container 'docker.io/staphb/kraken2:2.1.5'
    memory '10 GB'
    cpus 4
    time '1h'
    input:
    		each path("taxonomy")
        path('genomes/*')
    output:
    		path("db",type: 'dir')
    script:
		    """
		    mkdir -p db && ln -s ../taxonomy db/taxonomy
				find genomes/ -name '*.fna' -print0 | xargs -0 -P ${task.cpus} -I{} -n1 kraken2-build --add-to-library {} --db db
		    k2 build --db db ${task.ext.args?:''} --threads ${task.cpus}
		    rm -f db/taxonomy
		    kraken2-build --clean --db db
		    """
}



workflow KRAKEN2_DB {
  main:
		fa_ch = NCBI_DOWNLOAD_GENOMES(Channel.of(
			['Enterobacterales','--assembly-level complete'],
			['Pseudomonas aeruginosa','--assembly-level complete'],
			['Acinetobacter baumannii','--assembly-level complete'],
			['Enterococcus','--assembly-level complete'],
			['Staphylococcus','--assembly-level complete'],
			['Streptococcus','--assembly-level complete'],
			['Enterococcus faecalis',''],
			['Citrobacter murliniae','']
		))
		.collect()
		KRAKEN2_DB_DOWNLOAD_TAX()
		KRAKEN2_DB_BUILD(KRAKEN2_DB_DOWNLOAD_TAX.out,fa_ch)
  emit:
    db = KRAKEN2_DB_BUILD.out
}


