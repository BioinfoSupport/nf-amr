
# nf-amr

`nf-amr` is a nextflow pipeline to annotate FASTA assemblies for Antimicrobial Multi-Resistance.

The pipeline:

 1) Detect the species by comparing the assembly to a reference species database.
 
 2) Run `ResFinder` and `PlasmidFinder`.
 
 3) Run `MLST` with the schema automatically selected from detected species.

 4) Generate a HTML report.
 
## Usage

### Local computer

To run the pipeline on your local computer (with `docker`  and `nextflow` already installed):

```bash
nextflow run -r main BioinfoSupport/nf-amr -resume --input=data/r62b14.hdr.fasta
```


### HPC

To run the pipeline on a HPC cluster with `slurm` and `singularity` use `-profile=hpc`:

```bash
nextflow run -r main BioinfoSupport/nf-amr -profile hpc -resume
```

If `nextflow` is not installed on your HPC, it can be installed with:

```bash
ml Java/17.0.2
curl -s https://get.nextflow.io | bash
```





If `--input` argument is missing, the pipeline process all FASTA files located 
in subfolder `data/` (`--input='data/*.fasta'` by default).

If only `docker` is installed on your computer, but not `nextflow`, 
you can run a docker environment with `nextflow` inside:

```bash
docker run --rm -it \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd):$(pwd) \
  --platform linux/amd64 \
  --workdir $(pwd) \
  --env NXF_HOME=$(pwd)/.nextflow_home \
  nextflow/nextflow:24.10.4 bash
```

If in addition you need `nf-core` toolbox:

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


# Parking 

### Cleaning

```bash
rm -rf results/ .nextflow* work/
```

```bash
nextflow run . -resume
nextflow run . -resume --org_name="Citrobacter freundii"
```

# Test ANI speed
```bash
nextflow run . --input=data/r62b14.hdr.fasta --skip_plasmidfinder --skip_resfinder --skip_mlst --skip_prokka
```



```bash
docker run --rm nfcore/tools --help
docker run --rm nfcore/tools modules list remote
docker run --rm nfcore/tools pipelines list

docker run -it --rm --entrypoint bash nfcore/tools
nf-core pipelines create
nf-core pipelines create --name americ -a "J. Prados" --organisation "amr-genomics" --description "anti-microbial resistance infection control pipeline"
```





# TODO:
 - add sample sheet from CSV file (with id, org_name)

 - OrgDB and Org detection
   * custom ord_db is not functional yet because it uses absolute path
   * add module to generate db_tax.tsv
 
 - add multi-reporting / aggregator (which include independant report)

 - add ISfinder annotations (mobile elements)
 - add generation of CGView file
 - add hybracter assembly ?
 - add cgMLST analysis
 - add plasmid cgMLST
 
 - après assemblage ajouter les mapping des reads sur l'assemblage 
 
 


