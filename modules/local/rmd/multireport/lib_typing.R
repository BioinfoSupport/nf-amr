
library(tidyverse)
library(GenomicRanges)
library(Biostrings)


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
		relocate(contig_id,resistance_name=name,query_start_pos,query_end_pos,coverage,identity)
}
#fs::dir_ls("results/samples/",glob = "*/resfinder",recurse = 1,type = "dir") |> fs::path("data.json") |> head(1) |> read_resfinder_json()

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



################### END OF FILE ##########################


read_mlst_jsons <- function(json_files) {
	x <- enframe(json_files,value = "json_filename",name = "assembly_id") |>
		mutate(json_content = map(json_filename,~if (file.exists(.)) jsonlite::fromJSON(.) else NULL)) |>
		mutate(mlst_type=map(json_content,~.$mlst$results$sequence_type)) |> 
		select(assembly_id,mlst_type) |>
		unnest(mlst_type)
	x
}

read_mobtyper <- function(tsv_files) {
	enframe(tsv_files,value = "tsv_filename",name = "assembly_id") |>
		mutate(content = map(tsv_files,~read_tsv(.,show_col_types = FALSE))) |>
		unnest(content) |>
		mutate(contig_id=str_replace(sample_id," .*",""))
}

read_species_tsv <- function(tsv_file) {
	readr::read_tsv(tsv_file,col_types = cols(ANI = "n",bi_frag="i",query_frag="i",.default = "c"))
}

contig_meta <- function(f) {
	f <- as.character(f)
	fa <- Biostrings::readDNAStringSet(f)
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




db_load <- function(ss,dir) {
	message("Load contigs metadata (name, length, %GC, tags)")
	contigs <- ss |>
		mutate(assembly_name = file.path(dir,str_glue("{assembly_id}.hdr/assembly.fasta"))) |>
		pull("assembly_name",name = "assembly_id") |>
		BiocParallel::bplapply(\(f) if (file.exists(f)) contig_meta(f) else NULL) |>
		bind_rows(.id = "assembly_id")
	
	message("Load species informations")
	species <- ss |>
		mutate(species_tsv = file.path(dir,str_glue("{assembly_id}.hdr/org.ani"))) |>
		pull("species_tsv",name = "assembly_id") |>
		BiocParallel::bplapply(\(f) if (file.exists(f)) read_species_tsv(f) else NULL) |>
		bind_rows(.id = "assembly_id")
	
	message("Load resfinder results")
	resfinder <- ss |>
		mutate(species_tsv = file.path(dir,str_glue("{assembly_id}.hdr/resfinder/data.json"))) |>
		pull("species_tsv",name = "assembly_id") |>
		read_resfinder_jsons()
	
	message("Load plasmidfinder results")
	plasmidfinder <- ss |>
		mutate(plasmidfinder_json = file.path(dir,str_glue("{assembly_id}.hdr/plasmidfinder/data.json"))) |>
		pull("plasmidfinder_json",name = "assembly_id") |>
		read_plasmidfinder_jsons()

	message("Load MLST results")
	mlst <- ss |>
		mutate(mlst_json = file.path(dir,str_glue("{assembly_id}.hdr/mlst/data.json"))) |>
		pull("mlst_json",name = "assembly_id") |>
		read_mlst_jsons()
	
	message("Load Mobtyper results")
	mobtyper <- ss |>
		mutate(mobtyper_tsv = file.path(dir,str_glue("{assembly_id}.hdr/mobtyper.tsv"))) |>
		pull("mobtyper_tsv",name = "assembly_id") |>
		read_mobtyper()
	
	list(contigs=contigs,species=species,resfinder=resfinder,plasmidfinder=plasmidfinder,mlst=mlst,mobtyper=mobtyper)
}

