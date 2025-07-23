

docker run --rm -it \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd):$(pwd) \
  --platform linux/amd64 \
  --workdir $(pwd) \
  --env NXF_HOME=$(pwd)/.nextflow_home \
  nextflow/nextflow:25.04.6 bash




nextflow run BioinfoSupport/nf-amr -r v0.4.4 -resume --input_assembly='data/urgent/*.fasta' 

nextflow run . -resume -profile arm64 -output-dir 'results_run64_sup_v5.0.0' --input_assembly='data/run64_wf-bacterial/sup_v5.0.0/*.fasta*' 

nextflow run . -resume -profile arm64 --long_reads=data/run64_sup_v5.2.0/demuxed/demuxed/*_barcode17/reads.bam --flye_long

# Test hybrid assembly 
nextflow run . -resume --long_reads=data/test/long/RH1.fastq.gz --short_reads='data/test/short/RH1_*.fastq.gz' --flye_long --flye_hybrid


nextflow run . -resume --input data/r62b14.hdr.fasta


nextflow run . -resume --samplesheet=data/samplesheet_test.csv --unicycler_hybrid  --unicycler_short  --unicycler_long

nextflow run . -resume --long_reads=data/test/long/*.fastq.gz --short_reads=data/test/short/*.fastq.gz --input_assembly=data/test/assembly/*.fasta

nextflow run . -resume --fastq_long=data/test/long/RH2*.fastq.gz --fastq_short=data/test/short/RH2*.fastq.gz --input=data/test/assembly/RH2*.fasta

nextflow run . -resume --input='data/run99/*.fasta'

nextflow run . -resume --input data/oxa48_prd_ricai/assemblies/r2b3.fasta 


nextflow run . -resume --fastq_long=data/*.fastq.gz --fastq_short=data/*_R{1,2}.fastq.gz

nextflow run -r dev BioinfoSupport/nf-amr -resume --input='data/oxa48_prd_ricai/assemblies/*.fasta'

nextflow run . -resume --input='data/oxa48_prd_ricai/assemblies/r2b3.fasta'



docker run --rm -it registry.gitlab.unige.ch/amr-genomics/cgetools:main bash

# Merge dev branch with master branch
git checkout dev
git merge master
git checkout master
git merge dev
git tag -a v0.4.5 -m "This version include new Aeromonas"
git push origin --tags
git push origin
git checkout dev











ssh yggdrasil
cd ~/gvfs/*/BioinfoSupport/colombie
rsync -av 2025-07-09_sierra_triplecarba_run64 ~/scratch/
cd  ~/scratch/2025-07-09_sierra_triplecarba_run64

cd ~/gvfs/*/BioinfoSupport/colombie
rsync -av --no-perms ~/scratch/2025-07-09_sierra_triplecarba_run64/results ./2025-07-09_sierra_triplecarba_run64/



export NXF_SINGULARITY_CACHEDIR=~/scratch/singularity
export NXF_SINGULARITY_TMPDIR=~/scratch/singularity_tmp
nextflow run -r v0.4.5 BioinfoSupport/nf-amr -profile hpc -resume -bg --samplesheet=data/SampleSheet_pipeline.csv



scp yggdrasil:~/scratch/2025-07-09_sierra_triplecarba_run64/results/*.html ~/Downloads/

./nextflow run BioinfoSupport/nf-amr -profile hpc -bg -resume --input=data/*.fasta
