
# nf

Nextflow for AMR detection

```bash
docker run --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd):$(pwd) \
  --platform linux/amd64 \
  --workdir $(pwd) \
  nextflow/nextflow:24.10.1 \
  nextflow run ./ -resume
```

```bash
docker run --rm \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd):$(pwd) \
  --platform linux/amd64 \
  --workdir $(pwd) \
  nextflow/nextflow:24.10.1 \
  nextflow run -resume ./ -params-file params.json -profile local
```



```bash
docker run --rm nfcore/tools --help
docker run --rm nfcore/tools modules list remote
docker run --rm nfcore/tools pipelines list

docker run -it --rm --entrypoint bash nfcore/tools
nf-core pipelines create
nf-core pipelines create --name americ -a "J. Prados" --organisation "amr-genomics" --description "anti-microbial resistance infection control pipeline"
```




