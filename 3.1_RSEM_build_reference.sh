#!/bin/bash

# Define your base directory
basedir="/home/vagrant/RNASeq/"

# Define your Tools, Data and STAR directories
ToolsDir=${basedir}"Tools/"

# Add RSEM path to the PATH environment variable
export PATH=$PATH:${ToolsDir}RSEM-1.3.3/

# Create directory where RSEM reference files will be stored
mkdir -p ${ToolsDir}/RSEM-1.3.3/ref/

# Run rsem-prepare-reference
${ToolsDir}/RSEM-1.3.3/rsem-prepare-reference --gtf ${ToolsDir}"/ref_ucsc_hg38/hg38.refGene.gtf" ${ToolsDir}"/ref_ucsc_hg38/hg38.fa" ${ToolsDir}"RSEM-1.3.3/ref/ucsc_hg38"
