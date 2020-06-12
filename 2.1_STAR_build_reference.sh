#!/bin/bash
# Define your base directory
basedir="/home/vagrant/RNASeq/"

# Define your Tools directory and your Data directory (where you have your fastq files)
ToolsDir=${basedir}"Tools/"

# Create a folder for your reference, download and decompress genome and GTF files in it
mkdir -p ${ToolsDir}"/ref_ucsc_hg38"
cd ${ToolsDir}"/ref_ucsc_hg38/"
curl -LOJ "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" # hg38 fasta file with reference genome
gzip -d "hg38.fa.gz"
curl -LOJ "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz" # 'hg38.ensGene.gtf.gz' with Ensembl IDs is probably more complete, but 'refGene' contains directly gene names we are familiar with.
gzip -d "hg38.refGene.gtf.gz"

# Indexing with STAR. To use multiple threads with the option '--runThreadN'.
${ToolsDir}STAR-2.7.3a/bin/Linux_x86_64/STAR --runMode genomeGenerate --genomeDir ${ToolsDir}"/ref_ucsc_hg38/" --genomeFastaFiles ${ToolsDir}"/ref_ucsc_hg38/hg38.fa" --sjdbGTFfile ${ToolsDir}"/ref_ucsc_hg38/hg38.refGene.gtf"
