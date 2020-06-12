# clear workspace
rm(list = ls())

# Define your base directory
basedir = "/home/vagrant/RNASeq/"

# Define your RSEM directory
RsemOutDir = paste0(basedir, "Processed/RSEM_quantification/")

# Define list of samples
All_samples = c("Condition1_Rep1","Condition1_Rep2","Condition1_Rep3","Condition2_Rep1","Condition2_Rep2","Condition2_Rep3")

# Reading RSEM tables and creating two genes X samples tables, one with TPMs and one with expected counts
# We will just stich together the appropriate columns of the RSEM gene.results tables for all samples
# When doing this, always make sure that the ordering of the genes is consistent across all samples (see 'merge')!
is_first = TRUE
for (s in All_samples)
{
	this_table_path = paste0(RsemOutDir, s, ".genes.results")
	this_table = read.table(file = this_table_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
	if (is_first)
	{
		tpm_table = this_table[,c("gene_id","TPM")]
		colnames(tpm_table) = c("gene_id",s)

		counts_table = this_table[,c("gene_id","expected_count")]
		colnames(counts_table) = c("gene_id",s)

		is_first = FALSE
	} else {
		this_tpm_table = this_table[,c("gene_id","TPM")]
		colnames(this_tpm_table) = c("gene_id",s)
		tpm_table = merge(tpm_table,this_tpm_table, by = "gene_id")

		this_counts_table = this_table[,c("gene_id","expected_count")]
		colnames(this_counts_table) = c("gene_id",s)
		counts_table = merge(counts_table,this_counts_table, by = "gene_id")
	}
}
rownames(counts_table) = counts_table$gene_id
counts_table$gene_id = NULL
rownames(tpm_table) = tpm_table$gene_id
tpm_table$gene_id = NULL

# counts table will be used as input for differential expression analysis
write.table(counts_table, file = paste0(RsemOutDir, "tutorial_counts_table.txt"), col.names = TRUE, quote = FALSE, row.names = TRUE, sep = "\t" )

# we also save a table with TPMs. Sometimes such table is used as input for other downstream analyses (gene set enrichment analysis, bulk RNA-seq deconvolution, ...)
write.table(tpm_table, file = paste0(RsemOutDir, "tutorial_tpm_table.txt"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t" )
