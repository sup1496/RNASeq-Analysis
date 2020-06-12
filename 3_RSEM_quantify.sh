#!/bin/bash

# Define your base directory
basedir="/home/vagrant/RNASeq/"

# Define your Tools, Data and STAR directories
ToolsDir=${basedir}"Tools/"
DataDir=${basedir}"Data/"
StarOutDir=${basedir}"Processed/STAR_alignments/"

# Define and create your output directory, where results of RSEM will be saved, and cd into it
RsemOutDir=${basedir}"Processed/RSEM_quantification/"
mkdir -p ${RsemOutDir}
cd ${RsemOutDir}

# Add RSEM path to the PATH environment variable
export PATH=$PATH:${ToolsDir}RSEM-1.3.3/

# Run RSEM. We give it as input the transcriptome-aligned files from STAR, the RSEM reference and we specify that it was paired-end sequencing.
# Results will be saved in the current directory, and we specify 'SampleName' as prefix for output files.
SampleName="Condition1_Rep1"
${ToolsDir}/RSEM-1.3.3/rsem-calculate-expression --bam --paired-end --num-threads 16 ${StarOutDir}${SampleName}"_Aln_Aligned.toTranscriptome.out.bam" ${ToolsDir}"RSEM-1.3.3/ref/ucsc_hg38" ${SampleName}
SampleName="Condition1_Rep2"
${ToolsDir}/RSEM-1.3.3/rsem-calculate-expression --bam --paired-end --num-threads 16 ${StarOutDir}${SampleName}"_Aln_Aligned.toTranscriptome.out.bam" ${ToolsDir}"RSEM-1.3.3/ref/ucsc_hg38" ${SampleName}
SampleName="Condition1_Rep3"
${ToolsDir}/RSEM-1.3.3/rsem-calculate-expression --bam --paired-end --num-threads 16 ${StarOutDir}${SampleName}"_Aln_Aligned.toTranscriptome.out.bam" ${ToolsDir}"RSEM-1.3.3/ref/ucsc_hg38" ${SampleName}
SampleName="Condition2_Rep1"
${ToolsDir}/RSEM-1.3.3/rsem-calculate-expression --bam --paired-end --num-threads 16 ${StarOutDir}${SampleName}"_Aln_Aligned.toTranscriptome.out.bam" ${ToolsDir}"RSEM-1.3.3/ref/ucsc_hg38" ${SampleName}
SampleName="Condition2_Rep2"
${ToolsDir}/RSEM-1.3.3/rsem-calculate-expression --bam --paired-end --num-threads 16 ${StarOutDir}${SampleName}"_Aln_Aligned.toTranscriptome.out.bam" ${ToolsDir}"RSEM-1.3.3/ref/ucsc_hg38" ${SampleName}
SampleName="Condition2_Rep3"
${ToolsDir}/RSEM-1.3.3/rsem-calculate-expression --bam --paired-end --num-threads 16 ${StarOutDir}${SampleName}"_Aln_Aligned.toTranscriptome.out.bam" ${ToolsDir}"RSEM-1.3.3/ref/ucsc_hg38" ${SampleName}
