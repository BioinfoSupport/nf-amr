
# Introduction

This is a tutorial to run EPI2ME/wf-bacterial-genomes pipeline on UNIGE cluster.


# Step1: connect to cluster

```bash
ssh login1.yggrasil.hpc.unige.ch
cd ~/scratch/run49
```

# Step2: create config file

To run on UNIGE HPC, the pipeline requires the follwing `nextflow.config` in your working directory:

```bash
# Create this nextflow.config file in the working directory
cat <<EOF > nextflow.config
profiles {
	hpc {
		singularity {
			enabled = true
			autoMounts = true
		}
		process.queue = 'shared-cpu'
		process."withLabel:gpu".containerOptions = "--nv"
		process."withLabel:gpu".queue = 'shared-gpu'
		process."withLabel:gpu".clusterOptions = '--ntasks=1 --gpus-per-task=1'
		process.executor = "slurm"
		process.time = '4h'
	}
}
EOF
```


# Step3: Set nextflow environment variables

```bash
export NXF_SINGULARITY_CACHEDIR=~/scratch/singularity
export NXF_SINGULARITY_TMPDIR=~/scratch/singularity_tmp
```

**Note**: To avoid this step, this can be added in your `~/.bashrc`



# Step4: Run bacterial-genomes pipeline

```bash
nextflow run -r v1.4.2 -bg -resume \
  epi2me-labs/wf-bacterial-genomes \
  --isolates \
  --bam 'output_wf-basecalling_sup_v5.0.0/demuxed/demuxed/' \
  --out_dir 'output_wf-bacterial-genomes_sup_v5.0.0' \
  -profile hpc
```



