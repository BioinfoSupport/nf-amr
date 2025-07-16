

docker run --rm -it \
  -v /var/run/docker.sock:/var/run/docker.sock \
  -v $(pwd):$(pwd) \
  --platform linux/amd64 \
  --workdir $(pwd) \
  --env NXF_HOME=$(pwd)/.nextflow_home \
  nextflow/nextflow:25.04.6 bash




nextflow run . -resume -profile arm64 -output-dir 'results_run64_sup_v5.0.0' --input_assembly='data/run64_wf-bacterial/sup_v5.0.0/*.fasta' 

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
git tag -a v0.3 -m "This version include assembly options"
git push origin --tags
git checkout dev

