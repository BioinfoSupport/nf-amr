
# nf-amrseq

**This Page is under construction and not up to date**

`nf-amrseq` is a nextflow pipeline to process sequencing data of Antimicrobial Multi-Resistance bacteria.

The pipeline:

 1) Detect the organism by comparing the assembly to a reference species database with [`orgfinder`](https://gitlab.unige.ch/amr-genomics/orgfinder).
 
 2) Run `ResFinder`, `amrfinder+`, `mobtyper`, `PlasmidFinder`.
 
 3) Run `MLST` with the schema automatically selected from detected species. 

 4) Generate HTML reports.
 
 
## Usage

### Local computer

Running the pipeline on a local computer requires [`docker`](https://www.docker.com) 
(to run containerized software) and [`nextflow`](https://www.nextflow.io).
`nextflow` is however optional as containerized version exists (see section **Nextflow container** below).

**Note:** The pipeline cannot be run directly on a NAS but only on a local folder of your hard-drive.

If the FASTA files to process are in subfolder 'data/' of your working directory

```bash
nextflow run BioinfoSupport/nf-amr -resume --input_assembly=data/*.fasta
```

#### Nextflow container

Nextflow installation is optional as containerized versions of nextflow exists. 
For example, on MacOS, to run a nextflow container with access to docker daemon
of the host, we use:

```bash
docker run --rm -it \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd):$(pwd) \
  --platform linux/amd64 \
  --workdir $(pwd) \
  --env NXF_HOME=$(pwd)/.nextflow_home \
  nextflow/nextflow:25.04.2 bash
```

### HPC

To run the pipeline on a HPC cluster with `slurm` and `singularity` use `-profile=hpc`:

```bash
nextflow run BioinfoSupport/nf-amr -profile hpc -bg -resume --input_assembly=data/*.fasta
```

If `nextflow` is not installed on your HPC, it can be installed with:

```bash
curl -s https://get.nextflow.io | bash && chmod +x nextflow
```

And if you need to update the pipeline latest version, use:

```bash
nextflow pull BioinfoSupport/nf-amr
```

Or run a specific version:

```bash
nextflow run BioinfoSupport/nf-amr -r v0.4.12 -profile hpc -bg -resume --input_assembly=data/*.fasta
```


## Input

 - `--input_assembly`: FASTA files with assembled genomes to annotated 

## Output
By default, pipeline outputs are stored in `results/` but it can controlled with 
option `-output-dir` 


