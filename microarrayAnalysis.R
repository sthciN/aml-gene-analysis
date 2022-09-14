
R version 4.1.2 (2021-11-01) -- "Bird Hippie"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin17.0 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[R.app GUI 1.77 (8007) x86_64-apple-darwin17.0]

# Set working directory
setwd("/Users/ati/bio")

# Install packages
# install.packages(c("Biobase", 
#                   "GEOquery", 
#                   "limma", 
#                   "pheatmap",
#                   "reshape2",
#                   "plyr",
#                   "gplots",
#                   "ggplot2"
#                   ))

# Import packages
library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(reshape2)
library(plyr)
library(gplots)
library(ggplot2)

## aml <- read.delim("./data/GSE9476.top.table.tsv")
## head(aml)
## subset(aml, logFC > 1 & adj.P.Value < 0.05)
## aml.up <- subset(aml, logFC > 1 & adj.P.Val < 0.05)
## dim(aml.up)
## aml <- read.delim("GSE9476.top.table.tsv")
## patient <- subset(aml, logFC > 1 & adj.P.Val < 0.05)
## colnames(patient)
## dim(patient)
## length(patient$Gene.symbol)
## length(unique(patient$Gene.symbol))
## head(unique(patient$Gene.symbol))
## patient.gene = unique(patient$Gene.symbol)
## length(patient)
## length(patient.gene)


### Load Data

## Define series and platform
series <- "GSE9476"
platform <- "GPL96"

## Load the data with the series ID
gset <- getGEO(series,GSEMatrix=TRUE, AnnotGPL=TRUE, destdir="./data/")

## Or load the data from downloaded file 
# aml <- read.delim("GSE9476.top.table.tsv")

# You can see the length/class/names of your gset
length(gset) # 1
class(gset) # "list"
names(gset) # "GSE9476_series_matrix.txt.gz"

if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
gset


class(gset)

## Grouping the samples
gr <- c("CD34", 
        rep("BM", 10), 
        rep("CD34", 7), 
        rep("AML", 26), 
        rep("PB", 10), 
        rep("CD34", 10))
gr
length(gr)

### Expressing gset
ex = exprs(gset)
dim(ex) # 22283    64

## Plot the data
pdf("./result/boxplot.pdf", width=64)
boxplot(ex)
dev.off()


max(ex)
min(ex)

### Log2 scale, if required
# ex <- log2(ex + 1)
# exprs(gset) <- ex


#### Normalize, if required
# ex <- normalizeQuantiles(ex)
# exprs(gset) <- ex

### Correlation Heatmap
pdf("./result/core_heatmap.pdf", width=15, height=15)
pheatmap(cor(ex))
pheatmap(cor(ex), labels_row=gr, labels_col=gr)
dev.off()

### Principal Component Analysis
pc <- prcomp(ex)
pdf("./result/pc.pdf", width=15, height=15)
plot(pc)
dev.off()

names(pc) # [1] "sdev"     "rotation" "center"   "scale"    "x"
colnames(pc$x)
#  [1] "PC1"  "PC2"  "PC3"  "PC4"  "PC5"  "PC6"  "PC7"  "PC8"  "PC9"  "PC10"
# [11] "PC11" "PC12" "PC13" "PC14" "PC15" "PC16" "PC17" "PC18" "PC19" "PC20"
# [21] "PC21" "PC22" "PC23" "PC24" "PC25" "PC26" "PC27" "PC28" "PC29" "PC30"
# [31] "PC31" "PC32" "PC33" "PC34" "PC35" "PC36" "PC37" "PC38" "PC39" "PC40"
# [41] "PC41" "PC42" "PC43" "PC44" "PC45" "PC46" "PC47" "PC48" "PC49" "PC50"
# [51] "PC51" "PC52" "PC53" "PC54" "PC55" "PC56" "PC57" "PC58" "PC59" "PC60"
# [61] "PC61" "PC62" "PC63" "PC64"

pdf("./result/pc$x.pdf", width=15, height=15)
plot(pc$x[,1:2])
dev.off()


?scale
ex.scale <- t(scale(t(ex))) 
mean(ex.scale[1,]) # -1.671772e-16
head(ex)
# With scale FALSE
ex.scale <- t(scale(t(ex), scale=FALSE))

# Do the principal component analysis again
PCsc <- prcomp(ex.scale) 
pdf("./result/pcsc.pdf")
plot(PCsc)
plot(PCsc$x[,1:2])
dev.off()


pcrotation <- data.frame(PCsc$rotation[,1:3], Group=gr)
pdf("./result/pcrotation.pdf", width=15, height=15)
ggplot(pcrotation, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()


### Differential Expression Analysis

# Factoring the gr
gr <- factor(gr)
gr

# Create a design matrix that determine what is the group of each sample
gset$description <- gr
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gr)
head(design)

# Fit a linear model to gset based on the design matrix
fit <- lmFit(gset, design)

# Find the difference between AML and CD34 groups
cont.matrix <- makeContrasts (AML-CD34, levels=design)

# Find the slope of the fit based on the contrast
fit2 <- contrasts.fit(fit, cont.matrix)

# Apply bayesian model to get the P value
fit2 <- eBayes (fit2, 0.01)

# Benjamini-Hochberg: False Discovery Rate (fdr)
# The adjust can be changed to Bonferroni correction
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

# The topTable will provide 
tT <- subset (tT, select=c ("Gene.symbol", 
                            "Gene.ID", 
                            "adj.P.Val", 
                            "P.Value", 
                            "logFC"))
write.table(tT, "result/AML-CD34.txt", row.names=F, sep="\t", quote=F)


### Gene Antology

# Genes with P value < 0.05 and 
aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)

# Clean data , remove ///
# aml.up.genes <- sub("///.*", "", aml.up.genes)
# or split 
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))

dim(aml.up)
# [1] 566   5

write.table(aml.up.genes, 
            file="./result/AML_CD34_UP.txt", 
            quote=F, 
            row.names=F, 
            col.names=F)



aml.down <- subset(tT, logFC <-1 & adj.P.Val < 0.05)
# aml.down.genes = unique(aml.down$Gene.symbol)

## remove ///
# aml.down.genes <- sub("///.*", "", aml.down.genes)
## or split
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))

dim(aml.down)
write.table(aml.down.genes, 
            file="./result/AML_CD34_DOWN.txt", 
            quote=F,
            row.names=F,
            col.names=F)






