
# nf

Nextflow pipeline for AMR detection


## Usage

```bash
nextflow run -r main BioinfoSupport/nf-amr -resume
```

By default the pipeline process all FASTA files in subfolder `data` (`data/*.fasta`).



### Docker

If `nextflow``` is not installed on your system, you can run a docker environement
running it:

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
rm -rf results/
rm -rf .nextflow*
rm -rf work/
```

```bash
docker run --rm nfcore/tools --help
docker run --rm nfcore/tools modules list remote
docker run --rm nfcore/tools pipelines list

docker run -it --rm --entrypoint bash nfcore/tools
nf-core pipelines create
nf-core pipelines create --name americ -a "J. Prados" --organisation "amr-genomics" --description "anti-microbial resistance infection control pipeline"
```




