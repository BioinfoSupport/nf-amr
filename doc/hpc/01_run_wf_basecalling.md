
# Introduction

This is a tutorial to run EPI2ME/wf-basecalling pipeline on UNIGE cluster.


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



# Step4: Run basecalling pipeline

```bash
nextflow run -r v1.5.3 -bg -resume \
  epi2me-labs/wf-basecalling \
  -profile hpc,discrete_gpus \
  --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_sup@v5.0.0' \
  --barcode_kit 'SQK-NBD114-24' \
  --out_dir 'wf-basecalling_sup_v5.0.0' \
  --input 'pod5_pass' \
  --basecaller_chunk_size 2
```



