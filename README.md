
# nf-amr

`nf-amr` is a nextflow pipeline to annotate FASTA assemblies for Antimicrobial Multi-Resistance.

The pipeline:

 1) Detect the organism by comparing the assembly to a reference species database with [`orgfinder`](https://gitlab.unige.ch/amr-genomics/orgfinder).
 
 2) Run `ResFinder`, `amrfinder+`, `mobtyper`, `PlasmidFinder`.
 
 3) Run `MLST` with the schema automatically selected from detected species. 

 4) Generate a HTML report.
 
 
## Usage

### Local computer

Running the pipeline on a local computer requires [`docker`](https://www.docker.com) 
(to run containerized software) and [`nextflow`](https://www.nextflow.io).
`nextflow` is however optional as containerized version exists (see section **Nextflow container** below).

**Note:** The pipeline cannot be run directly on a NAS but only on a local folder of your hard-drive.

If the FASTA files to process are in subfolder 'data/' of your working directory

```bash
nextflow run BioinfoSupport/nf-amr -resume --input=data/*.fasta
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
ml Nextflow/24.04.2
nextflow run BioinfoSupport/nf-amr -profile hpc -resume
```

If `nextflow` is not installed on your HPC, it can be installed with:

```bash
ml Java/17.0.2
curl -s https://get.nextflow.io | bash
```



## Input

If `--input` argument is missing, the pipeline process all FASTA files located 
in subfolder `data/` (`--input='data/*.fasta'` by default).



# TODO:
 - Setup all tools so they use external DB given in param: orgfinder, plasmidfinder, resfinder, amrfinder, mlst
 - add sample sheet from CSV file (with id, org_name)
 - add multi-reporting / aggregator (which include independant report)

 - add ISfinder annotations (mobile elements)
 - add generation of CGView file
 - add hybracter assembly ?
 - add cgMLST analysis
 - add plasmid cgMLST
 - après assemblage ajouter les mapping des reads sur l'assemblage 
 
 



# Parking 


### Cleaning

```bash
rm -rf results/ .nextflow* work/
```

```bash
nextflow run . -resume
nextflow run . -resume --org_name="Citrobacter freundii"
```

### Test ANI speed
```bash
nextflow run . --input=data/r62b14.hdr.fasta --skip_plasmidfinder --skip_resfinder --skip_mlst --skip_prokka
```


### nf-core
```bash
docker run --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd):$(pwd) \
  --platform linux/amd64 \
  --workdir $(pwd) \
  --env NXF_HOME=$(pwd)/.nextflow_home \
  nfcore/tools \
  pipelines schema build
```

```bash
docker run --rm nfcore/tools --help
docker run --rm nfcore/tools modules list remote
docker run --rm nfcore/tools pipelines list

docker run -it --rm --entrypoint bash nfcore/tools
nf-core pipelines create
nf-core pipelines create --name americ -a "J. Prados" --organisation "amr-genomics" --description "anti-microbial resistance infection control pipeline"
```





