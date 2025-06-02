
process ORGFINDER_DB_DOWNLOAD_GENOMES {
		container 'docker.io/staphb/ncbi-datasets:16.30.0'
		memory '6 GB'
		cpus 1
		output:
		    path("genomes", type: 'dir')
		script:
"""
mkdir -p genomes/

# Download some complete reference genomes
datasets download genome taxon 'Pseudomonas aeruginosa' --reference --assembly-level complete --filename Paeruginosa.zip && unzip Paeruginosa.zip -d genomes/Paeruginosa
datasets download genome taxon 'Acinetobacter baumannii' --reference --assembly-level complete --filename Abaumannii.zip && unzip Abaumannii.zip -d genomes/Abaumannii
datasets download genome taxon 'Enterococcus' --reference --assembly-level complete --filename Enterococcus.zip && unzip Enterococcus.zip -d genomes/Enterococcus
datasets download genome taxon 'Staphylococcus' --reference --assembly-level complete --filename Staphylococcus.zip && unzip Staphylococcus.zip -d genomes/Staphylococcus
datasets download genome taxon 'Streptococcus' --reference --assembly-level complete --filename Streptococcus.zip && unzip Streptococcus.zip -d genomes/Streptococcus

# Download additional uncomplete reference genomes
datasets download genome taxon 'Enterococcus faecalis' --reference --filename Efaecalis.zip && unzip Efaecalis.zip -d genomes/Efaecalis
datasets download genome taxon 'Citrobacter murliniae' --reference --filename Cmurliniae.zip && unzip Cmurliniae.zip -d genomes/Cmurliniae

# Download Enterobacterales genomes
datasets download genome taxon 'Enterobacterales' --reference --assembly-level complete --filename Enterobacterales.zip && unzip Enterobacterales.zip -d genomes/Enterobacterales

# Make the tsv file with all accession numbers
cat genomes/*/ncbi_dataset/data/assembly_data_report.jsonl | dataformat tsv genome --force --fields accession,organism-name,organism-tax-id,assmstats-total-sequence-len,assmstats-total-number-of-chromosomes > genomes/db_accession.tsv

# Retreive taxonomy and untar
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && mkdir -p genomes/taxdump && tar -C genomes/taxdump -zxf taxdump.tar.gz
"""
}



process ORGFINDER_DB_BUILD {
		memory '6 GB'
		cpus 1
		input:
		    path("genomes", type: 'dir')
		output:
		    path("db", type: 'dir')
		script:
		"""
		orgfinder_db_build.R
		"""
}

workflow ORGFINDER_DB {
	ORGFINDER_DB_DOWNLOAD_GENOMES 
	| ORGFINDER_DB_BUILD
}




