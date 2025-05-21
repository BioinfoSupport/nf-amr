

library(tidyverse)


read_resfinder_json <- function() {
	jsonlite::fromJSON("resfinder.json")\$seq_regions |>
		tibble() |>
		unnest_wider(1) |>
		mutate(contig_id = str_replace(query_id," .*","")) |>
		mutate(assembly_id = "${meta.id}")
}

read_plasmidfinder_json <- function() {
	jsonlite::fromJSON("plasmidfinder.json") |>
		pluck("plasmidfinder","results") |>
		enframe("plasmidfinder_Species") |> 
		unnest_longer(value,indices_to = "plasmidfinder_species") |>
		unnest_longer(value) |>
		select(!value_id) |>
		unnest_wider(value) |>
		mutate(contig_id = str_replace(contig_name," .*",""))	|>
		mutate(assembly_id = "${meta.id}") |>
		relocate(assembly_id,contig_id)
}


read_mlst_json <- function() {
	jsonlite::fromJSON("mlst.json")\$mlst\$results\$sequence_type |>
		enframe(value="mlst_type") |>
		mutate(assembly_id = "${meta.id}")
}

