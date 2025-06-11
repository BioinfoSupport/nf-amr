
library(tidyverse)
library(GenomicRanges)
library(Biostrings)

read_runinfo_json <- function(json_file) {
	#json_file <- "results/samples/r62b17.hdr/runinfo.json"
	json <- jsonlite::fromJSON(json_file,simplifyVector=FALSE)
	tibble(json) |> 
		unnest_wider(1) |> 
		unnest_wider(c(meta,orgfinder),names_sep=".") |>
		mutate(across(any_of("orgfinder.ani"),as.numeric))
}


read_resfinder_json <- function(json_file) {
	json <- jsonlite::fromJSON(json_file,simplifyVector=FALSE)
	expected_structure <- tibble(
		query_id = character(0),
		name = character(0),
		query_start_pos = numeric(0),
		query_end_pos = numeric(0),
		coverage = numeric(0),
		identity = numeric(0)
	)
	json$seq_regions |>
		tibble() |> 
		unnest_wider(1) |>
		mutate(contig_id = str_replace(query_id," .*","")) |>
		bind_rows(expected_structure) |>
		mutate(position=str_glue("{query_start_pos}-{query_end_pos}")) |>
		relocate(contig_id,resistance_name=name,query_start_pos,query_end_pos,coverage,identity,position)
}
#fs::dir_ls("results/samples/",glob = "*/resfinder",recurse = 1,type = "dir") |> fs::path("data.json") |> head(1) |> read_resfinder_json()


read_amrfinderplus_tsv <- function(tsv_file) {
	#tsv_file <- "results/samples/r62b17.hdr/amrfinderplus/report.tsv"
	read_tsv(tsv_file,show_col_types = FALSE) |>
		dplyr::rename(contig_id=`Contig id`,resistance_name=`Element symbol`,coverage=`% Coverage of reference`,identity=`% Identity to reference`) |>
		mutate(position=str_glue("{Start}-{Stop}:{Strand}")) |>
		relocate(contig_id,coverage,identity,resistance_name,position)
}


read_plasmidfinder_json <- function(json_file) {
	#json_file <- "results/samples/r62b17.hdr/plasmidfinder/data.json"
	#json_file <- "results/samples/RS1357_contig_1/plasmidfinder/data.json"
	json <- jsonlite::fromJSON(json_file,simplifyVector=FALSE)
	expected_structure <- tibble(
		contig_name = character(0),
		plasmid = character(0),
		coverage = numeric(0),
		identity = numeric(0)
	)
	enframe(json$plasmidfinder$results,name = "db_lev1") |>
		unnest_longer(value,indices_to = "db_lev2") |>
		filter(value != "No hit found") |>		
		unnest_longer(value,indices_include = FALSE) |>
		unnest_wider(value) |>
		bind_rows(expected_structure) |>
		mutate(contig_id = str_replace(contig_name," .*","")) |>
		relocate(contig_id,plasmid_type=plasmid,coverage,identity)
}
#fs::dir_ls("results/samples/",glob = "*/plasmidfinder",recurse = 1,type = "dir") |> fs::path("data.json") |> tail(1) |> read_plasmidfinder_json()



read_cgemlst_json <- function(json_file) {
	if (!file.exists(json_file)) return(tibble(mlst_type="?"))
	#json_file <- "results/samples/r62b17.hdr/cge_mlst/data.json"
	json <- jsonlite::fromJSON(json_file,simplifyVector=FALSE)
	json$mlst$results$sequence_type |>
		enframe(value = "mlst_type",name = NULL) |>
		mutate(mlst_type=if_else(mlst_type %in% c("Unknown",NA),"?",mlst_type))
}
#fs::dir_ls("results/samples/",glob = "*/mlst",recurse = 1,type = "dir") |> fs::path("data.json") |> tail(1) |> read_mlst_json()


read_mobtyper_tsv <- function(tsv_file) {
		expected_structure <- tibble(
			contig_id = character(0),
			contig_name = character(0),
			relaxase_types = character(0),
			rep_types = character(0)
		)	
		if (!file.exists(tsv_file)) return(expected_structure)
		#tsv_file <- "results/samples/r62b17.hdr/mobtyper.tsv"
		read_tsv(tsv_file,show_col_types = FALSE) |>
			mutate(contig_id=str_replace(sample_id," .*","")) |>
			dplyr::rename(contig_name=sample_id,relaxase_types=`relaxase_type(s)`,rep_types=`rep_type(s)`) |>
			bind_rows(expected_structure) |>
			relocate(contig_id)
}


read_orgfinder_tax <- function(tsv_file) {
	#tsv_file <- "results/samples/r62b17.hdr/orgfinder/tax.tsv"
	read_tsv(tsv_file,show_col_types = FALSE,col_types = cols(.default = "c"))	
}

read_orgfinder_tsv <- function(tsv_file) {
	#tsv_file <- "results/samples/r62b17.hdr/orgfinder/ani.tsv"
	read_tsv(tsv_file,show_col_types = FALSE,col_types = cols(ANI = "n",bi_frag="i",query_frag="i",.default = "c"))
}


contig_meta <- function(fasta_filename) {
	fasta_filename <- as.character(fasta_filename)
	fa <- Biostrings::readDNAStringSet(fasta_filename)
	meta <- tibble(
		contig_id = str_replace(names(fa)," .*",""),
		contig_length = lengths(fa),
		header = names(fa),
		GC = as.vector(letterFrequency(fa,"GC",as.prob=TRUE))
	)
	stopifnot("FASTA contains duplicated contig_id" = all(!duplicated(meta$contig_id)))
	
	# regex-pattern for modifier and header line
	tag_pat <- "((\\[([^=]+)=([^\\]]+)\\])|(([^=]+)=([^\\s]+)))"
	hdr_pat <- str_glue("^([^\\s]+)((\\s*{tag_pat})*)\\s*(.*)$")
	if (!all(str_detect(meta$header,hdr_pat))) rlang::abort("Some header lines are malformed")
	
	meta <- meta |>
		mutate(title = str_extract(header,hdr_pat,group=11)) |>
		mutate(tags_str=str_trim(str_extract(header,hdr_pat,group=2)))
	
	
	# Extract key-value pairs from modifiers
	tags <- tibble(select(meta,contig_id,tags_str)) |>
		mutate(tags = str_extract_all(tags_str,tag_pat),tags_str = NULL) |>
		unnest(tags) |>
		mutate(tags = str_trim(tags)) |>
		# TODO: use coalesce() here ?
		mutate(
			tag_name = if_else(is.na(str_extract(tags,tag_pat,3)),str_extract(tags,tag_pat,6),str_extract(tags,tag_pat,3)),
			tag_value = if_else(is.na(str_extract(tags,tag_pat,4)),str_extract(tags,tag_pat,7),str_extract(tags,tag_pat,4)),
			tags = NULL
		)
	
	stopifnot("Some modifiers are set several times in the same header" = !(summarize(tags,.by=c(contig_id,tag_name),any(n()>1)) |> pull() |> any()))
	tags <- pivot_wider(tags,id_cols="contig_id",names_from = "tag_name",values_from = "tag_value",names_prefix = "tag_")
	
	left_join(meta,tags,by="contig_id",relationship = "one-to-one") |>
		select(!c(tags_str,header))
}



db_load <- function(amr_dir) {
	fs::dir_ls(amr_dir,recurse = 1,glob = "*/assembly.fasta") |>
			fs::path_dir() |>
			enframe(name = NULL,value = "basepath") |>
			mutate(assembly_id = basename(basepath)) |>
			mutate(runinfo = map(fs::path(basepath,"runinfo.json"),read_runinfo_json)) |>
			mutate(orgfinder = map(fs::path(basepath,"orgfinder","tax.tsv"),read_orgfinder_tax)) |>
			mutate(contigs = map(fs::path(basepath,"assembly.fasta"),contig_meta)) |>
			mutate(mlst = map(fs::path(basepath,"cge_mlst","data.json"),read_cgemlst_json)) |>
			mutate(plasmidfinder = map(fs::path(basepath,"plasmidfinder","data.json"),read_plasmidfinder_json)) |>
			mutate(resfinder = map(fs::path(basepath,"resfinder","data.json"),read_resfinder_json)) |>
			mutate(amrfinderplus = map(fs::path(basepath,"amrfinderplus","report.tsv"),read_amrfinderplus_tsv)) |>		
			mutate(mobtyper = map(fs::path(basepath,"mobtyper.tsv"),read_mobtyper_tsv))
}


summarise_assembly <- function(db) {
	#db <- db_load("results") 
	assemlbies <- db |> 
		select(assembly_id,contigs) |> 
		unnest(contigs) |> 
		group_by(assembly_id) |> 
		summarise(num_contig=n(),assembly_length=sum(contig_length),GC=weighted.mean(GC,contig_length),N50=N50(contig_length)) 
	mlst <- db |> select(assembly_id,mlst) |> unnest(mlst)
	runinfo <- db |> select(assembly_id,runinfo) |> unnest(runinfo) 
	orgfinder <- db |> select(assembly_id,orgfinder) |> unnest(orgfinder) |> select(assembly_id,org_name,species_name,genus_name)
	assemlbies |>
		left_join(runinfo,by="assembly_id",relationship = "one-to-one")	|>
		left_join(mlst,by="assembly_id",relationship = "one-to-one") |>
		left_join(orgfinder,by=c("assembly_id","org_name"),relationship = "many-to-one") |>
		left_join(rename_with(orgfinder,.cols=!assembly_id,~str_c("orgfinder.",.)),by=c("assembly_id","orgfinder.org_name"),relationship = "one-to-one")
}

summarise_resistances <- function(db) {
	#db <- db_load("results")
	bind_rows(
		resfinder = db |> 
			select(assembly_id,resfinder) |> 
			unnest(resfinder) |>
			select(assembly_id,contig_id,resistance_name,coverage,identity,position),
		amrfinderplus = db |> 
			select(assembly_id,amrfinderplus) |> 
			unnest(amrfinderplus) |>
			select(assembly_id,contig_id,resistance_name,coverage,identity,position),
		.id = "source"
	) |>
		arrange(assembly_id,contig_id,resistance_name,desc(coverage),desc(identity)) |>
		group_by(assembly_id,contig_id,resistance_name) |>
		slice_head(n = 1) |>
		ungroup()
}

summarise_contigs <- function(db) {
	#db <- db_load("results") 
	contigs <- db |> 
		select(assembly_id,contigs) |> 
		unnest(contigs)
	res <- summarise_resistances(db) |>
		group_by(assembly_id,contig_id) |>
		summarise(resistance_names = list(resistance_name))
	plf <- db |> 
		select(assembly_id,plasmidfinder) |> 
		unnest(plasmidfinder) |>
		group_by(assembly_id,contig_id) |>
		summarise(plasmid_types=list(plasmid_type))
	mob <- db |> 
		select(assembly_id,mobtyper) |> 
		unnest(mobtyper) |>
		mutate(relaxase_types=str_split(relaxase_types,",")) |>
		select(assembly_id,contig_id,relaxase_types)
	contigs |>
		select(assembly_id,contig_id,contig_length,GC,topology=tag_topology) |>
		left_join(res,by=c("assembly_id","contig_id"),relationship = "one-to-one") |>
		left_join(plf,by=c("assembly_id","contig_id"),relationship = "one-to-one") |>
		left_join(mob,by=c("assembly_id","contig_id"),relationship = "one-to-one")
}






