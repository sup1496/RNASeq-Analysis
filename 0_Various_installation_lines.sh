#!/bin/bash

#### This script is part of "RNA-seq data analysis: from raw reads to differential expression" tutorial of keepsciencegoing.org initiative. 
#### You can find explanations and references in the related videos under keepsciencegoing.org
#### The script is provided "as is" for educational purposes. 
#### Author: Daniele Tavernari ( daniele.tavernari@unil.ch ), Date: 12th April 2020

### You will not need admin rights

basedir="/path_to_your_base_directory/"

# Move to the Tools directory
cd ${basedir}"/Tools/"

# Download FastQC with curl - you can also download it with your browser
curl -LOJ "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip"
unzip "fastqc_v0.11.9.zip"
chmod 777 "FastQC/fastqc"

# Download STAR
curl -LOJ "https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz"
tar -xzf "STAR-2.7.3a.tar.gz"

# Download RSEM
curl -LOJ "https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz"
tar -xzf "RSEM-1.3.3.tar.gz"
cd "RSEM-1.3.3"
make
