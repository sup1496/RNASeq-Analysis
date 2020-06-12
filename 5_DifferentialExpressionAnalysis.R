
# clear workspace
rm(list = ls())

# loading libraries
library(limma)
library(edgeR)

# Define your base directory
basedir = "/home/vagrant/RNASeq/"

# Define your RSEM directory
RsemOutDir = paste0(basedir, "Processed/RSEM_quantification/")

# Reading the counts table
counts_table = read.table(file = paste0(RsemOutDir, "tutorial_counts_table.txt"), header = TRUE, row.names = 1, quote = '', sep = "\t", stringsAsFactors = F )

# Setting up annotation table
All_samples = colnames(counts_table)
annot = data.frame(row.names = All_samples, Sample = All_samples, Treatment = substr(All_samples,1,10), stringsAsFactors = F)

x = DGEList(counts = counts_table, group = annot[colnames(counts_table),"Treatment"], remove.zeros = TRUE , genes = rownames(counts_table))
group = as.factor(annot[colnames(counts_table),"Treatment"])

# Scale transformation
cpm = cpm(x)
lcpm = cpm(x, log=TRUE)

# Removing genes that are lowly expressed
keep.exprs = filterByExpr(x, group = group) 
x = x[keep.exprs, keep.lib.sizes=FALSE]
nrow(x)
plot(density(lcpm))
plot(density(cpm(x, log = TRUE)))

# Normalising for library size
x = calcNormFactors(x, method = "TMM")

# Unsupervised clustering (can be explored by group or other variable, also using different components)
plotMDS(x, labels = group)

# Setting up design matrix and contrasts
design = model.matrix(~0+group)
colnames(design) = substr(colnames(design), 6, nchar(colnames(design)))
design
contr.matrix = makeContrasts(
   Condition1_vs_Condition2 = Condition1-Condition2,
   levels = colnames(design))

# Removing heteroscedasticity
xvoom = voom(x, design, plot=TRUE)

# Fitting linear model
vfit = lmFit(xvoom, design)
vfit = contrasts.fit(vfit, contrasts=contr.matrix) 
efit = eBayes(vfit)

# Extracting differentially expressed genes
top = topTable(efit, number = Inf, adjust.method="BH", sort.by = "p")

adj_pval_thresh = 0.05
abs_logFC_thresh = 1

top$DiffExpr_status = "not_significant"
top[ (top$logFC > abs_logFC_thresh) & (top$adj.P.Val < adj_pval_thresh),"DiffExpr_status"] = "up_in_cond1"
top[ (top$logFC < (-abs_logFC_thresh) ) & (top$adj.P.Val < adj_pval_thresh),"DiffExpr_status"] = "up_in_cond2"

# Heatmap
library(gplots)
de_genes = top[top$DiffExpr_status!="not_significant","genes"]

heatmap.2(lcpm[de_genes,], scale = 'row', col = colorRampPalette(c("blue","white","red"))(100), dendrogram = "column", trace = "none", cexCol = 0.5 )

# Volcano plot
plot(x = top$logFC, y = -log10(top$P.Value) )
plot(x = top$logFC, y = -log10(top$P.Value), pch = 19, cex = 0.5, 
     xlab = "logFC (1 vs 2)", ylab = "-log10(p-value)")

top$color = "gray"
top[top$DiffExpr_status=="up_in_cond1","color"] = "red"
top[top$DiffExpr_status=="up_in_cond2","color"] = "blue"
plot(x = top$logFC, y = -log10(top$P.Value), pch = 19, cex = 0.5, 
     xlab = "logFC (1 vs 2)", ylab = "-log10(p-value)", col = top$color)

# Write top table to a file
write.table(top, file = paste0(RsemOutDir, "top_tab.txt"), row.names = T, col.names = T, sep = "\t", quote = F)

