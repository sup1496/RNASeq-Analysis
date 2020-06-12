#!/bin/bash

# Define your base directory
basedir="/home/vagrant/RNASeq/"

# Define your Tools directory (where you unzipped FastQC) and your Data directory (where you have your fastq files)
ToolsDir=${basedir}"Tools/"
DataDir=${basedir}"Data/"

# Define and create your output directory, where results of FASTQC will be saved
FastqcDir=${basedir}"Processed/FastQC/"
mkdir -p ${FastqcDir} # Flag '-p' creates the directory only if that directory does not exist already

# Define sample file name. Note that I don't include "_R1.fastq.gz" or "_R2.fastq.gz", I will concatenate them to the file name when I will run FastQC
SampleName="Condition1_Rep1"

# Run FastQC and save results into the output directory
${ToolsDir}/FastQC/fastqc ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --extract --outdir ${FastqcDir}

# Run FastQC on all the other samples. To keep things simple, here I am just copy-pasting and changing the 'SampleName'.
# If you have a lot of samples and you are becoming skilled with bash scripting, you can instead implement a for loop over all your samples.
SampleName="Condition1_Rep2"
${ToolsDir}/FastQC/fastqc ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --extract --outdir ${FastqcDir}
SampleName="Condition1_Rep3"
${ToolsDir}/FastQC/fastqc ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --extract --outdir ${FastqcDir}
SampleName="Condition2_Rep1"
${ToolsDir}/FastQC/fastqc ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --extract --outdir ${FastqcDir}
SampleName="Condition2_Rep2"
${ToolsDir}/FastQC/fastqc ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --extract --outdir ${FastqcDir}
SampleName="Condition2_Rep3"
${ToolsDir}/FastQC/fastqc ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --extract --outdir ${FastqcDir}
