#!/bin/bash

# Define your base directory
basedir="/home/vagrant/RNASeq/"

# Define your Tools directory and your Data directory (where you have your fastq files)
ToolsDir=${basedir}"Tools/"
DataDir=${basedir}"Data/"

# Define and create your output directory, where results of STAR will be saved, and cd into it
StarOutDir=${basedir}"Processed/STAR_alignments/"
mkdir -p ${StarOutDir}
cd ${StarOutDir}

# Align reads with STAR. If they are paired-end, the two R1 and R2 files should be provided one after the other.
# By default, results will be saved in the current directory. We specify the file name root with '--outFileNamePrefix'
# We are providing also the reference built previously.
# We can give as input to STAR compressed fastq files. The option "--readFilesCommand gunzip -c" tells STAR to decompress the .gz fastq files with zcat "on the fly", without having to occupy disk space by uncompressing the fastq files beforehand.
# The option '--outSAMtype BAM SortedByCoordinate' allows to store output alignments in (sorted) .bam files, which are lighter than SAM files.
# The option '--quantMode TranscriptomeSAM' will output alignments translated into transcript coordinates, that will be used by RSEM for quantification of expression at the gene and transcript level.
# Clearly, there are several other options that can be set and tuned. See STAR reference manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
SampleName="Condition1_Rep1"
${ToolsDir}STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir ${ToolsDir}"/ref_ucsc_hg38/" --readFilesIn ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --readFilesCommand 'zcat' --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${SampleName}"_Aln_"
SampleName="Condition1_Rep2"
${ToolsDir}STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir ${ToolsDir}"/ref_ucsc_hg38/" --readFilesIn ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --readFilesCommand 'zcat' --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${SampleName}"_Aln_"
SampleName="Condition1_Rep3"
${ToolsDir}STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir ${ToolsDir}"/ref_ucsc_hg38/" --readFilesIn ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --readFilesCommand 'zcat' --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${SampleName}"_Aln_"
SampleName="Condition2_Rep1"
${ToolsDir}STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir ${ToolsDir}"/ref_ucsc_hg38/" --readFilesIn ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --readFilesCommand 'zcat' --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${SampleName}"_Aln_"
SampleName="Condition2_Rep2"
${ToolsDir}STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir ${ToolsDir}"/ref_ucsc_hg38/" --readFilesIn ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --readFilesCommand 'zcat' --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${SampleName}"_Aln_"
SampleName="Condition2_Rep3"
${ToolsDir}STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 16 --genomeDir ${ToolsDir}"/ref_ucsc_hg38/" --readFilesIn ${DataDir}${SampleName}"_R1.fastq.gz" ${DataDir}${SampleName}"_R2.fastq.gz" --readFilesCommand 'zcat' --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${SampleName}"_Aln_"
