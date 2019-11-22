#NoAb_Sham - light blue
#Ab_Sham - dark grey
#NoAB_TBI - blue
#Ab_TBI - orange

#Do:
#DONE - FDR of 0.05
#DONE - heatmap top 40 up and down defined by noAb_TBI v. Ab_TBI
#GO term / BART analysis for same top 40 up/down
#square PCA plot and change colors, no difference in transparency
#x-ref inflammatory genes
#vignette!

library(lattice)
library(DESeq2)
library(pheatmap)

source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

setwd('~/Desktop/Ashley_RNAseq')

categorize.deseq.df <- function(df, fdr = 0.05, log2fold = 0.0, treat
= 'Auxin') {

     df.activated = data.frame(matrix(nrow = 0, ncol = 0))
     df.repressed = data.frame(matrix(nrow = 0, ncol = 0))

     if (nrow(df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold,]) != 0) {
     	df.activated = df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold,]
	df.activated$arfauxin = paste(treat, 'Activated')
	}

     if (nrow(df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold,]) != 0) {
     	df.repressed = df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold,]
	df.repressed$arfauxin = paste(treat, 'Repressed')
	}
    
    df.unchanged = df[df$padj > 0.5 & !is.na(df$padj) & abs(df$log2FoldChange) < 0.25,]
    df.unchanged$arfauxin = paste(treat, 'Unchanged')

    df.dregs = df[!(df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold) &
                  !(df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold) &
                  !(df$padj > 0.5 & !is.na(df$padj) &
    		  abs(df$log2FoldChange) < 0.25), ]
    df.dregs$arfauxin = paste(treat, 'All Other Genes')
    
    df.effects.lattice = 
    rbind(df.activated, 
          df.unchanged, 
          df.repressed, 
          df.dregs)
    df.effects.lattice$arfauxin = factor(df.effects.lattice$arfauxin)
    df.effects.lattice$arfauxin = relevel(df.effects.lattice$arfauxin, ref = paste(treat, 'Unchanged'))
    df.effects.lattice$arfauxin = relevel(df.effects.lattice$arfauxin, ref = paste(treat, 'All Other Genes'))
    return(df.effects.lattice)
}

run.deseq.list.any <- function(mat, unt = 2, trt =2) {
    sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))        
    deseq.counts.table = DESeqDataSetFromMatrix(mat, DataFrame(sample.conditions), ~ sample.conditions);
    colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
    dds = DESeq(deseq.counts.table)
    return(dds)
}

#These two codes create a data frame for all the libraries, x is the gene names
x <- read.table("Ab_Sham_rep1.gene.counts.txt",header=FALSE)
all.counts <- data.frame(row.names = x$V1)

#Putting name of all *.combined.gene.counts.txt file into a vector. For loop then append individual gene counts to merged.counts. 
samples <- c("Ab_Sham_rep1","Ab_Sham_rep2","Ab_Sham_rep3","Ab_Sham_rep4","Ab_TBI_rep1","Ab_TBI_rep2","Ab_TBI_rep3","Ab_TBI_rep4","noAb_Sham_rep1",
             "noAb_Sham_rep2","noAb_Sham_rep3","noAb_Sham_rep4","noAb_TBI_rep1","noAb_TBI_rep2","noAb_TBI_rep3","noAb_TBI_rep4")

for (i in samples){
  j<-read.table(print(paste0(i,".gene.counts.txt")))
  all.counts <- cbind(all.counts,data.frame(i = j[2]))
}

colnames(all.counts) <- samples

#Delete Last 5 lines of non-feature counts
merged.counts <- all.counts[-(55422:55426),]

#Replace Ensembl ID with Gene Name as rownames

#Generating Gene ID alongside EnsembleID (CAUTION: rm() these data frames when you're done because they take up a lot of memory to store)
ensembl.all = read.table('Mus_musculus.GRCm38.98.corrected.gtf', sep='\t', header =F);
ensembl.gene.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_name '), "[", 2), ";"), "[", 1));
ensembl.id.names = data.frame(sapply(strsplit(sapply(strsplit(as.character(ensembl.all[,9]),'gene_id '), "[", 2), ";"), "[", 1));
ensembl.code = cbind(ensembl.gene.names, ensembl.id.names);
ensembl.code = ensembl.code[!duplicated(ensembl.code[,2]),];
rownames(ensembl.code) = ensembl.code[,2];
colnames(ensembl.code) = c('gene', 'id');

save(ensembl.code, file = "ensembl.code.Rdata")

rm(ensembl.all)
rm(ensembl.gene.names)
rm(ensembl.id.names)

merged.counts = merge(merged.counts, ensembl.code, by="row.names", all.x=F)
rownames(merged.counts) <- make.names(merged.counts$gene, unique=TRUE)
merged.counts <- merged.counts[,-c(1,18,19)]

save(merged.counts, file="merged.counts.Rdata")

rm(ensembl.code)

#PCA

merged_dds = run.deseq.list.any(merged.counts, unt=4, trt=12)
rld_HH = rlogTransformation(merged_dds)
pca.dat = plotPCA(rld_HH, intgroup="condition", returnData=TRUE)
percentVar = round(100 * attr(pca.dat, "percentVar"))
plotPCAlattice(pca.dat, file = 'PCA_RNAseq.pdf')

###### Ab_Sham v. Ab_TBI ######

merged.counts.small = merged.counts[,1:8]

# number of replicates per condition
unt = 4
trt = 4

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)

mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)

ab.sham.v.ab.tbi.lattice = 
    categorize.deseq.df(res.mm.atac, 
                        fdr = 0.05, log2fold = 0.0, treat = 'TBI')

pdf("ab.sham.v.ab.tbi.MAplot.pdf", useDingbats = FALSE, width=4, height=3.33);
print(xyplot(ab.sham.v.ab.tbi.lattice$log2FoldChange ~ 
                 log(ab.sham.v.ab.tbi.lattice$baseMean, base=10), 
             groups=ab.sham.v.ab.tbi.lattice$arfauxin,
             col=c("grey90", "grey60", "red", "blue"),
             scales="free",
             aspect=1,
             ylim=c(-6.5, 6.5),
#             xlim=c(-1,4.2),
             par.strip.text=list(cex=1.0, font = 1),
             pch=20,
             cex=0.5,
             ylab=expression("+TBI log"[2]~"Expression fold change"), 
             xlab=expression("log"[10]~"Mean of Normalized Counts"),
             main = 'Ab_Sham v. Ab_TBI',
            par.settings=list(par.xlab.text=list(cex=1.1,font=2), 
                              par.ylab.text=list(cex=1.1,font=2), 
                              strip.background=list(col="grey85"))))
dev.off()

save(ab.sham.v.ab.tbi.lattice, file='ab.sham.v.ab.tbi.lattice.Rdata')

###### noAb_TBI v. Ab_TBI ######

merged.counts.small = merged.counts[,c(13:16,5:8)]

# number of replicates per condition
unt = 4
trt = 4

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)

mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)

noAb.tbi.v.ab.tbi.lattice = 
    categorize.deseq.df(res.mm.atac, 
                        fdr = 0.05, log2fold = 0.0, treat = 'TBI')

pdf("noAb.tbi.v.ab.tbi.MAplot.pdf", useDingbats = FALSE, width=4, height=3.33);
print(xyplot(noAb.tbi.v.ab.tbi.lattice$log2FoldChange ~ 
                 log(noAb.tbi.v.ab.tbi.lattice$baseMean, base=10), 
             groups=noAb.tbi.v.ab.tbi.lattice$arfauxin,
             col=c("grey90", "grey60", "red", "blue"),
             scales="free",
             aspect=1,
             ylim=c(-6.5, 6.5),
#             xlim=c(-1,4.2),
             par.strip.text=list(cex=1.0, font = 1),
             pch=20,
             cex=0.5,
             ylab=expression("+Ablation log"[2]~"Expression fold change"), 
             xlab=expression("log"[10]~"Mean of Normalized Counts"),
             main = 'noAb_TBI v. Ab_TBI',
            par.settings=list(par.xlab.text=list(cex=1.1,font=2), 
                              par.ylab.text=list(cex=1.1,font=2), 
                              strip.background=list(col="grey85"))))
dev.off()

save(noAb.tbi.v.ab.tbi.lattice, file='noAb.tbi.v.ab.tbi.lattice.Rdata')

###### noAb_Sham v. noAb_TBI ######

merged.counts.small = merged.counts[,c(9:12,13:16)]

# number of replicates per condition
unt = 4
trt = 4

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)

mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)

noAb.sham.v.noAb.tbi.lattice = 
    categorize.deseq.df(res.mm.atac, 
                        fdr = 0.05, log2fold = 0.0, treat = 'TBI')

pdf("noAb.sham.v.noAb.tbi.MAplot.pdf", useDingbats = FALSE, width=4, height=3.33);
print(xyplot(noAb.sham.v.noAb.tbi.lattice$log2FoldChange ~ 
                 log(noAb.sham.v.noAb.tbi.lattice$baseMean, base=10), 
             groups=noAb.sham.v.noAb.tbi.lattice$arfauxin,
             col=c("grey90", "grey60", "red", "blue"),
             scales="free",
             aspect=1,
             ylim=c(-6.5, 6.5),
#             xlim=c(-1,4.2),
             par.strip.text=list(cex=1.0, font = 1),
             pch=20,
             cex=0.5,
             ylab=expression("+TBI log"[2]~"Expression fold change"), 
             xlab=expression("log"[10]~"Mean of Normalized Counts"),
             main = 'noAb_sham v. noAb_TBI',
            par.settings=list(par.xlab.text=list(cex=1.1,font=2), 
                              par.ylab.text=list(cex=1.1,font=2), 
                              strip.background=list(col="grey85"))))
dev.off()

save(noAb.sham.v.noAb.tbi.lattice, file='noAb.sham.v.noAb.tbi.lattice.Rdata')

###### noAb_Sham v. Ab_Sham ######

merged.counts.small = merged.counts[,c(9:12,1:4)]

# number of replicates per condition
unt = 4
trt = 4

sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))
mm.deseq.counts.table = DESeqDataSetFromMatrix(merged.counts.small, DataFrame(sample.conditions), ~ sample.conditions)

mm.atac = mm.deseq.counts.table
atac.size.factors = estimateSizeFactorsForMatrix(merged.counts.small)

sizeFactors(mm.atac) = atac.size.factors
mm.atac = estimateDispersions(mm.atac)
mm.atac = nbinomWaldTest(mm.atac)
res.mm.atac = results(mm.atac)

noAb.sham.v.Ab.sham.lattice = 
    categorize.deseq.df(res.mm.atac, 
                        fdr = 0.05, log2fold = 0.0, treat = 'TBI')

pdf("noAb.sham.v.Ab.sham.MAplot.pdf", useDingbats = FALSE, width=4, height=3.33);
print(xyplot(noAb.sham.v.Ab.sham.lattice$log2FoldChange ~ 
                 log(noAb.sham.v.Ab.sham.lattice$baseMean, base=10), 
             groups=noAb.sham.v.Ab.sham.lattice$arfauxin,
             col=c("grey90", "grey60", "red", "blue"),
             scales="free",
             aspect=1,
             ylim=c(-6.5, 6.5),
#             xlim=c(-1,4.2),
             par.strip.text=list(cex=1.0, font = 1),
             pch=20,
             cex=0.5,
             ylab=expression("+Ablation log"[2]~"Expression fold change"), 
             xlab=expression("log"[10]~"Mean of Normalized Counts"),
             main = 'noAb_sham v. Ab_sham',
            par.settings=list(par.xlab.text=list(cex=1.1,font=2), 
                              par.ylab.text=list(cex=1.1,font=2), 
                              strip.background=list(col="grey85"))))
dev.off()

save(noAb.sham.v.Ab.sham.lattice, file='noAb.sham.v.Ab.sham.lattice.Rdata')

#for GO term analysis:
x = noAb.tbi.v.ab.tbi.lattice
x = na.omit(x)

y = x[x$log2FoldChange > 0,]
y = y[order(y$padj),]
inc = rownames(y)[1:40]
write(inc,file='top40_increased_genes.txt',sep='\t')

y = x[x$log2FoldChange < 0,]
y = y[order(y$padj),]
dec = rownames(y)[1:40]
write(dec,file='top40_decreased_genes.txt',sep='\t')

#heat maps

unt=4
trt=12
sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))        
 
deseq.counts.table = DESeqDataSetFromMatrix(merged.counts, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)

a = assay(rld_HH)[c(inc[1:20],dec[1:20]),]
a = a - rowMeans(a)
a = a[,c(13:16,5:8)]

pdf(file='Ab.v.noAb.top20.heatmap.pdf')
print(pheatmap(a, cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE))
dev.off()

#complement gene set
x = read.table('complement_hallmark.txt',sep='\t')
x = as.vector(x[-c(1,2),])

genes = c()
for (i in x) {
    genes = append(genes,grep(paste0('^',i,'$'),toupper(rownames(noAb.tbi.v.ab.tbi.lattice))))
}

lattice = noAb.tbi.v.ab.tbi.lattice[genes,]
lattice = lattice[order(lattice$padj),]

a = assay(rld_HH)[rownames(lattice)[1:20],]
a = a - rowMeans(a)
a = a[,c(13:16,5:8)]

pdf(file='complement.top20.heatmap.pdf')
print(pheatmap(a, cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE))
dev.off()

#inflammation gene set
x = read.table('inflammation_hallmark.txt',sep='\t')
x = as.vector(x[-c(1,2),])

genes = c()
for (i in x) {
    genes = append(genes,grep(paste0('^',i,'$'),toupper(rownames(noAb.tbi.v.ab.tbi.lattice))))
}

lattice = noAb.tbi.v.ab.tbi.lattice[genes,]
lattice = lattice[order(lattice$padj),]

a = assay(rld_HH)[rownames(lattice)[1:40],]
a = a - rowMeans(a)
a = a[,c(13:16,5:8)]

pdf(file='inflammation.top40.heatmap.pdf')
print(pheatmap(a, cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE))
dev.off()

#volcano plots (p-value < 0.05)
lattice = noAb.sham.v.noAb.tbi.lattice

pdf(file = 'noAb.sham.v.no.Ab.tbi.volcano.plot.pdf')
plot(lattice$log2FoldChange, -log10(lattice$pvalue), cex = 0.8, pch = 20, main = 'TBI with Ablation v. TBI without Ablation Volcano Plot', xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
with(subset(lattice, pvalue < 0.05 & log2FoldChange < 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='blue'))
with(subset(lattice, pvalue < 0.05 & log2FoldChange > 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='red'))
abline(h=-log(0.05,base=10), lty = 5, lwd = 4, col = '#C0C0C0')
dev.off()

lattice = ab.sham.v.ab.tbi.lattice

pdf(file = 'ab.sham.v.ab.tbi.volcano.plot.pdf')
plot(lattice$log2FoldChange, -log10(lattice$pvalue), cex = 0.8, pch = 20, main = 'TBI with Ablation v. TBI without Ablation Volcano Plot', xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
with(subset(lattice, pvalue < 0.05 & log2FoldChange < 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='blue'))
with(subset(lattice, pvalue < 0.05 & log2FoldChange > 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='red'))
abline(h=-log(0.05,base=10), lty = 5, lwd = 4, col = '#C0C0C0')
dev.off()

lattice = noAb.tbi.v.ab.tbi.lattice

pdf(file = 'noAb.tbi.v.ab.tbi.volcano.plot')
plot(lattice$log2FoldChange, -log10(lattice$pvalue), cex = 0.8, pch = 20, main = 'TBI with Ablation v. TBI without Ablation Volcano Plot', xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
with(subset(lattice, pvalue < 0.05 & log2FoldChange < 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='blue'))
with(subset(lattice, pvalue < 0.05 & log2FoldChange > 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='red'))
abline(h=-log(0.05,base=10), lty = 5, lwd = 4, col = '#C0C0C0')
dev.off()

#volcano plot w/ highlights
genes = c('^Arc$', '^C1qc$', '^C1qa$', '^Itgam$', '^C4b$', '^Fn1$')

row.nums = c()
for(i in genes) { 
row.nums = append(row.nums, grep(i,rownames(lattice)))
}

lattice$color = 'black'
lattice[lattice$arfauxin == 'Ablation Activated',]$color = 'red'
lattice[lattice$arfauxin == 'Ablation Repressed',]$color = 'blue'
lattice[row.nums,]$color = 'green'

pdf(file = 'noAb.tbi.v.ab.tbi.highlights.volcano.plot.pdf')

plot(lattice$log2FoldChange, -log10(lattice$pvalue), col =
lattice$color, cex = 0.8, pch = 20, main = 'TBI with Ablation v. TBI without Ablation Volcano Plot', xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
abline(h=-log(0.05,base=10), lty = 5, lwd = 4, col = '#C0C0C0')

dev.off()

##### BAO ######

    #Heatmap for Correlation btw Replicates
    library(gplots)
    pdf("RelationRep.pdf", height = 8, width = 8.5)
    distance <- dist(t(assay(rld_HH)))
    bluwhite <- colorRampPalette(c("Blue3","White"))
    heatmap.2(as.matrix(distance), col = bluwhite, trace = "none", breaks = 200,
              main="Distance between Replicates", lwid = c(5,20), margins = c(11,12), density.info = "none",
              keysize = 1.2, lhei=c(2,10), key.xlab = "Distance")
    dev.off()
    

    process.deseq.df(four_100nM_vs_DMSO,"4hr_100nM_Processed",adj=.001,pval = NA, fc = 0)
    process.deseq.df(four_10nM_vs_DMSO,"4hr_10nM_Processed",adj=.001,pval = NA, fc = 0)
    process.deseq.df(twentyfour_100nM_vs_DMSO,"24hr_100nM_Processed",adj=.001,pval = NA, fc = 0)
    process.deseq.df(twentyfour_10nM_vs_DMSO,"24hr_10nM_Processed",adj=.001,pval = NA, fc = 0)
    
    
#Converting dds to a dataframe and adding differential status label
    #Note that we are using the dplyr package to convert the rownames(gene names) into an 
    #actual column in the dataframe, which makes manipulation in R easier
library(dplyr)
    
df_4_100 <- as.data.frame(four_100nM_vs_DMSO)
df_4_100 <- df_4_100 %>% rownames_to_column("gene")
df_4_100 = df_4_100[order(df_4_100$gene),]
df_4_100 <- mutate(df_4_100, status = factor (case_when(df_4_100$padj < .001 & df_4_100$log2FoldChange > .25 ~ "increased",
                                                    df_4_100$padj < .001 & df_4_100$log2FoldChange < -.25 ~ "decreased",
                                                    TRUE ~ "unchanged")))

df_4_10 <- as.data.frame(four_10nM_vs_DMSO)
df_4_10 <- df_4_10 %>% rownames_to_column("gene")
df_4_10 = df_4_10[order(df_4_10$gene),]
df_4_10 <- mutate(df_4_10, status = factor (case_when(df_4_10$padj < .001 & df_4_10$log2FoldChange > .25 ~ "increased",
                                                        df_4_10$padj < .001 & df_4_10$log2FoldChange < -.25 ~ "decreased",
                                                        TRUE ~ "unchanged")))

df_24_100 <- as.data.frame(twentyfour_100nM_vs_DMSO)
df_24_100 <- df_24_100 %>% rownames_to_column("gene")
df_24_100 = df_24_100[order(df_24_100$gene),]
df_24_100 <- mutate(df_24_100, status = factor (case_when(df_24_100$padj < .001 & df_24_100$log2FoldChange > .25 ~ "increased",
                                                      df_24_100$padj < .001 & df_24_100$log2FoldChange < -.25 ~ "decreased",
                                                      TRUE ~ "unchanged")))

df_24_10 <- as.data.frame(twentyfour_10nM_vs_DMSO)
df_24_10 <- df_24_10 %>% rownames_to_column("gene")
df_24_10 = df_24_10[order(df_24_10$gene),]
df_24_10 <- mutate(df_24_10, status = factor (case_when(df_24_10$padj < .001 & df_24_10$log2FoldChange > .25 ~ "increased",
                                                          df_24_10$padj < .001 & df_24_10$log2FoldChange < -.25 ~ "decreased",
                                                          TRUE ~ "unchanged")))

save(df_24_10, df_24_100, df_4_10, df_4_100, file = "analysis_df.Rdata")
load("analysis_df.Rdata")


#Generating Subset for Comparing Activation/Repression between different conditions
par(mar=c(10,5,5,4))
dev.off()
#Significant peaks from 4_100 graphed by log2FoldChange from 24_100
df_4_100_vs_24_100 <- cbind(df_4_100, df_24_100$log2FoldChange, df_24_100$padj)
colnames(df_4_100_vs_24_100)[8:9] <- c("foldchange_24_100", "p_24_100")

  #ACTIVATED
  df_4_100_vs_24_100_act <- subset (df_4_100_vs_24_100, padj < .001 & log2FoldChange > 0) 
    df_4_100_vs_24_100_act_act <- subset(df_4_100_vs_24_100_act, p_24_100 < .001 & foldchange_24_100 > 0)
    df_4_100_vs_24_100_act_unchanged <- subset(df_4_100_vs_24_100_act, p_24_100 > .2 & abs(foldchange_24_100) < .2)
    df_4_100_vs_24_100_act_other <- subset(df_4_100_vs_24_100_act, p_24_100 < .2 & p_24_100 > .001 | p_24_100 > .2 & abs(foldchange_24_100) > .2)
    df_4_100_vs_24_100_act_dec <- subset(df_4_100_vs_24_100_act, p_24_100 < .001 & foldchange_24_100 < 0)
    
    boxplot(df_4_100_vs_24_100_act_other$foldchange_24_100, 
            df_4_100_vs_24_100_act_unchanged$foldchange_24_100, 
            df_4_100_vs_24_100_act_act$foldchange_24_100,
            df_4_100_vs_24_100_act_dec$foldchange_24_100,
            main = "All 4hr 100nM Romi Activated Genes",
            col=c("lightgrey","darkgrey","firebrick2","dodgerblue"),
            ylab="24hr 100nM Log2 Fold Change",ylim=c(-5,10.5), las=2, xaxt="n")
    lablist.x <- c("24hr 100nM Other","24hr 100nM Unchanged","24hr 100nM Activated","24hr 100nM Repressed")
    axis(4, at=seq(-5, 10, by=5), labels = FALSE)
    text(x = seq(1.1, 4.4, by=1.05), par("usr")[3] - .5, labels = lablist.x, srt = 60, pos = 2, xpd = TRUE, col = c("lightgrey","darkgrey","firebrick2","dodgerblue"))
    abline(h=0,lty=3)
    
  #REPRESSED   
  df_4_100_vs_24_100_dec <- subset(df_4_100_vs_24_100, padj < .001 & log2FoldChange < 0)
    df_4_100_vs_24_100_dec_dec <- subset(df_4_100_vs_24_100_dec, p_24_100 < .001 & foldchange_24_100 < 0)
    df_4_100_vs_24_100_dec_unchanged <- subset(df_4_100_vs_24_100_dec, p_24_100 > .2 & abs(foldchange_24_100) < .2)
    df_4_100_vs_24_100_dec_other <- subset(df_4_100_vs_24_100_dec, p_24_100 < .2 & p_24_100 > .001 | p_24_100 > .2 & abs(foldchange_24_100) > .2)
    df_4_100_vs_24_100_dec_act <- subset(df_4_100_vs_24_100_dec, p_24_100 < .001 & foldchange_24_100 > 0)
    
    boxplot(df_4_100_vs_24_100_dec_other$foldchange_24_100, 
            df_4_100_vs_24_100_dec_unchanged$foldchange_24_100, 
            df_4_100_vs_24_100_dec_act$foldchange_24_100,
            df_4_100_vs_24_100_dec_dec$foldchange_24_100,
            main = "All 4hr 100nM Romi Repressed Genes",
            col=c("lightgrey","darkgrey","firebrick2","dodgerblue"),
            ylab="24hr 100nM Log2 Fold Change",ylim=c(-8,3), las=2, xaxt="n")
    lablist.x <- c("24hr 100nM Other","24hr 100nM Unchanged","24hr 100nM Activated","24hr 100nM Repressed")
    axis(4, at=seq(-8, 2, by=2), labels = FALSE)
    text(x = seq(1.1, 4.4, by=1.05), par("usr")[3] - .4, labels = lablist.x, srt = 60, pos = 2, xpd = TRUE, col = c("lightgrey","darkgrey","firebrick2","dodgerblue"))
    abline(h=0,lty=3)

#Significant peaks from 4_10 graphed by log2FoldChange from 4_100
df_4_10_vs_4_100 <- cbind(df_4_10, df_4_100$log2FoldChange, df_4_100$padj)
colnames(df_4_10_vs_4_100)[8:9] <- c("foldchange_4_100", "p_4_100")
    
    #ACTIVATED
    df_4_10_vs_4_100_act <- subset (df_4_10_vs_4_100, padj < .001 & log2FoldChange > 0) 
    df_4_10_vs_4_100_act_act <- subset(df_4_10_vs_4_100_act, p_4_100 < .001 & foldchange_4_100 > 0)
    df_4_10_vs_4_100_act_unchanged <- subset(df_4_10_vs_4_100_act, p_4_100 > .2 & abs(foldchange_4_100) < .2)
    df_4_10_vs_4_100_act_other <- subset(df_4_10_vs_4_100_act, p_4_100 < .2 & p_4_100 > .001 | p_4_100 > .2 & abs(foldchange_4_100) > .2)
    df_4_10_vs_4_100_act_dec <- subset(df_4_10_vs_4_100_act, p_4_100 < .001 & foldchange_4_100 < 0)
    
    boxplot(df_4_10_vs_4_100_act_other$foldchange_4_100, 
            df_4_10_vs_4_100_act_unchanged$foldchange_4_100, 
            df_4_10_vs_4_100_act_act$foldchange_4_100,
            df_4_10_vs_4_100_act_dec$foldchange_4_100,
            main = "All 4hr 10nM Romi Activated Genes",
            col=c("lightgrey","darkgrey","firebrick2","dodgerblue"),
            ylab="4hr 100nM Log2 Fold Change",ylim=c(-1.5,4), las=2, xaxt="n")
    lablist.x <- c("4hr 100nM Other","4hr 100nM Unchanged","4hr 100nM Activated","4hr 100nM Repressed")
    axis(4, at=seq(-1, 4, by=1), labels = FALSE)
    text(x = seq(1.1, 4.4, by=1.05), par("usr")[3] - .5, labels = lablist.x, srt = 60, pos = 2, xpd = TRUE, col = c("lightgrey","darkgrey","firebrick2","dodgerblue"))
    abline(h=0,lty=3)
    
    #REPRESSED   
    df_4_10_vs_4_100_dec <- subset(df_4_10_vs_4_100, padj < .001 & log2FoldChange < 0)
    df_4_10_vs_4_100_dec_dec <- subset(df_4_10_vs_4_100_dec, p_4_100 < .001 & foldchange_4_100 < 0)
    df_4_10_vs_4_100_dec_unchanged <- subset(df_4_10_vs_4_100_dec, p_4_100 > .2 & abs(foldchange_4_100) < .2)
    df_4_10_vs_4_100_dec_other <- subset(df_4_10_vs_4_100_dec, p_4_100 < .2 & p_4_100 > .001 | p_4_100 > .2 & abs(foldchange_4_100) > .2)
    df_4_10_vs_4_100_dec_act <- subset(df_4_10_vs_4_100_dec, p_4_100 < .001 & foldchange_4_100 > 0)
    
    boxplot(df_4_10_vs_4_100_dec_other$foldchange_4_100, 
            df_4_10_vs_4_100_dec_unchanged$foldchange_4_100, 
            df_4_10_vs_4_100_dec_act$foldchange_4_100,
            df_4_10_vs_4_100_dec_dec$foldchange_4_100,
            main = "All 4hr 10nM Romi Repressed Genes",
            col=c("lightgrey","darkgrey","firebrick2","dodgerblue"),
            ylab="4hr 100nM Log2 Fold Change",ylim=c(-3,1), las=2, xaxt="n")
    lablist.x <- c("4hr 100nM Other","4hr 100nM Unchanged","4hr 100nM Activated","4hr 100nM Repressed")
    axis(4, at=seq(-3, 1, by=1), labels = FALSE)
    text(x = seq(1.1, 4.4, by=1.05), par("usr")[3] - .4, labels = lablist.x, srt = 60, pos = 2, xpd = TRUE, col = c("lightgrey","darkgrey","firebrick2","dodgerblue"))
    abline(h=0,lty=3)    
    
        
#Significant peaks from 24_100 graphed by log2FoldChange from 4_100    
df_24_100_vs_4_100 <- cbind(df_24_100, df_4_100$log2FoldChange, df_4_100$padj)
colnames(df_24_100_vs_4_100)[8:9] <- c("foldchange_4_100", "p_4_100")

    #ACTIVATED
    df_24_100_vs_4_100_act <- subset (df_24_100_vs_4_100, padj < .001 & log2FoldChange > 0) 
    df_24_100_vs_4_100_act_act <- subset(df_24_100_vs_4_100_act, p_4_100 < .001 & foldchange_4_100 > 0)
    df_24_100_vs_4_100_act_unchanged <- subset(df_24_100_vs_4_100_act, p_4_100 > .2 & abs(foldchange_4_100) < .2)
    df_24_100_vs_4_100_act_other <- subset(df_24_100_vs_4_100_act, p_4_100 < .2 & p_4_100 > .001 | p_4_100 > .2 & abs(foldchange_4_100) > .2)
    df_24_100_vs_4_100_act_dec <- subset(df_24_100_vs_4_100_act, p_4_100 < .001 & foldchange_4_100 < 0)
    
    boxplot(df_24_100_vs_4_100_act_other$foldchange_4_100, 
            df_24_100_vs_4_100_act_unchanged$foldchange_4_100, 
            df_24_100_vs_4_100_act_act$foldchange_4_100,
            df_24_100_vs_4_100_act_dec$foldchange_4_100,
            main = "All 24hr 100nM Romi Activated Genes", 
            col=c("lightgrey","darkgrey","firebrick2","dodgerblue"),
            ylab="4hr 100nM Log2 Fold Change",ylim=c(-2.5,8), las=2, xaxt="n")
    lablist.x <- c("4hr 100nM Other","4hr 100nM Unchanged","4hr 100nM Activated","4hr 100nM Repressed")
    axis(4, at=seq(-2, 8, by=2), labels = FALSE)
    text(x = seq(1.1, 4.4, by=1.05), par("usr")[3] - .4, labels = lablist.x, srt = 60, pos = 2, xpd = TRUE, col = c("lightgrey","darkgrey","firebrick2","dodgerblue"))
    abline(h=0,lty=3)
    
    #REPRESSED   
    df_24_100_vs_4_100_dec <- subset (df_24_100_vs_4_100, padj < .001 & log2FoldChange < 0)
    df_24_100_vs_4_100_dec_dec <- subset(df_24_100_vs_4_100_dec, p_4_100 < .001 & foldchange_4_100 < 0)
    df_24_100_vs_4_100_dec_unchanged <- subset(df_24_100_vs_4_100_dec, p_4_100 > .2 & abs(foldchange_4_100) < .2)
    df_24_100_vs_4_100_dec_other <- subset(df_24_100_vs_4_100_dec, p_4_100 < .2 & p_4_100 > .001 | p_4_100 > .2 & abs(foldchange_4_100) > .2)
    df_24_100_vs_4_100_dec_act <- subset(df_24_100_vs_4_100_dec,  p_4_100 < .001 & foldchange_4_100 > 0)
    
    boxplot(df_24_100_vs_4_100_dec_other$foldchange_4_100, 
            df_24_100_vs_4_100_dec_unchanged$foldchange_4_100, 
            df_24_100_vs_4_100_dec_act$foldchange_4_100,
            df_24_100_vs_4_100_dec_dec$foldchange_4_100,
            main = "All 24hr 100nM Romi Repressed Genes", 
            col=c("lightgrey","darkgrey","firebrick2","dodgerblue"),
            ylab="4hr 100nM Log2 Fold Change",ylim=c(-4,2), las=2, xaxt="n")
    lablist.x <- c("4hr 100nM Other","4hr 100nM Unchanged","4hr 100nM Activated","4hr 100nM Repressed")
    axis(4, at=seq(-4, 2, by=1), labels = FALSE)
    text(x = seq(1.1, 4.4, by=1.05), par("usr")[3]-.2, labels = lablist.x, srt = 60, pos = 2, xpd = TRUE, col = c("lightgrey","darkgrey","firebrick2","dodgerblue"))
    abline(h=0,lty=3)

###GENE ONTOLOGY--upload these files to PANTHER to test for enrichment
write.table(df_4_10_vs_4_100_act$gene, "increase_GO_4_10.txt",row.names=F, col.names = F,quote=FALSE)
write.table(df_4_10_vs_4_100_dec$gene,"decrease_GO_4_10.txt", row.names=F, col.names = F, quote=FALSE)
    
write.table(df_4_100_vs_24_100_act$gene, "increase_GO_4_100.txt",row.names=F, col.names = F,quote=FALSE)
write.table(df_4_100_vs_24_100_dec$gene,"decrease_GO_4_100.txt", row.names=F, col.names = F, quote=FALSE)

write.table(df_24_10_vs_4_100_act$gene, "increase_GO_24_10.txt",row.names=F, col.names = F,quote=FALSE)
write.table(df_24_10_vs_4_100_dec$gene,"decrease_GO_24_10.txt", row.names=F, col.names = F, quote=FALSE)

write.table(df_24_100_vs_4_100_act$gene, "increase_GO_24_100.txt",row.names=F, col.names = F,quote=FALSE)
write.table(df_24_100_vs_4_100_dec$gene,"decrease_GO_24_100.txt", row.names=F, col.names = F, quote=FALSE)

par(mar=c(5,32,5,2))
dev.off()

#The panther.txt file comes with several rows of nonsense; delete these rows before running the file in R.
#BIOLOGICAL PROCESSES
GO_inc <- read.table("GO_4_100_bio_inc.txt", header=T, sep="\t")
GO_inc <- GO_inc[,c(-2)]
colnames(GO_inc) <- c("Process","Count","Expected","plusminus","foldchange","pvalue","fdr")
GO_inc$foldchange <- as.numeric(sub("<", "", GO_inc$foldchange))
GO_inc <- GO_inc[GO_inc$fdr<.001,]
GO_inc <- GO_inc[order(-GO_inc$foldchange),]
barplot(GO_inc$foldchange[15:1], names.arg = GO_inc$Process[15:1], las=2, xlim=c(0,3), 
                              col = "firebrick2", main = "Biological Processes Enriched in 4hr 100nM Romi Activated Genes                                                                         ", 
                               horiz=T, las=1, xlab = "Enrichment"); box();

GO_dec <- read.table("GO_4_100_bio_dec.txt", header=T, sep="\t")
GO_dec <- GO_dec[,c(-2)]
colnames(GO_dec) <- c("Process","Count","Expected","plusminus","foldchange","pvalue","fdr")
GO_dec$foldchange <- as.numeric(sub("<", "", GO_dec$foldchange))
GO_dec <- GO_dec[GO_dec$fdr<.001,]
GO_dec <- GO_dec[order(-GO_dec$foldchange),]
barplot(GO_dec$foldchange[15:1], names.arg = GO_dec$Process[15:1], las=2, xlim=c(0,6), 
                              col = "dodgerblue", main = "Biological Processes Enriched in 4hr 100nM Romi Repressed Genes                                                                           ", 
                              horiz=T, las=1, xlab = "Enrichment"); box();

#CELLULAR COMPONENTS
GO_inc <- read.table("GO_4_100_cell_inc.txt", header=T, sep="\t")
GO_inc <- GO_inc[,c(-2)]
colnames(GO_inc) <- c("Process","Count","Expected","plusminus","foldchange","pvalue","fdr")
GO_inc$foldchange <- as.numeric(sub("<", "", GO_inc$foldchange))
GO_inc <- GO_inc[GO_inc$fdr<.001,]
GO_inc <- GO_inc[order(-GO_inc$foldchange),]
barplot(GO_inc$foldchange[15:1], names.arg = GO_inc$Process[15:1], las=2, xlim=c(0,3), 
                              col = "firebrick2", main = "Cellular Components Enriched in 4hr 100nM Romi Activated Genes                                                                       ", 
                              horiz=T, las=1, xlab = "Enrichment"); box()

GO_dec <- read.table("GO_4_100_cell_dec.txt", header=T, sep="\t")
GO_dec <- GO_dec[,c(-2)]
colnames(GO_dec) <- c("Process","Count","Expected","plusminus","foldchange","pvalue","fdr")
GO_dec$foldchange <- as.numeric(sub("<", "", GO_dec$foldchange))
GO_dec <- GO_dec[GO_dec$fdr<.001,]
GO_dec <- GO_dec[order(-GO_dec$foldchange),]
barplot(GO_dec$foldchange[15:1], names.arg = GO_dec$Process[15:1], las=2, xlim=c(0,5), 
                              col = "dodgerblue", main = "Cellular Components Enriched in 4hr 100nM Romi Repressed Genes                                                                              ", 
                              horiz=T, las=1, xlab = "Enrichment"); box();

#MOLECULAR FUNCTIONS
GO_inc <- read.table("GO_4_100_func_inc.txt", header=T, sep="\t")
GO_inc <- GO_inc[,c(-2)]
colnames(GO_inc) <- c("Process","Count","Expected","plusminus","foldchange","pvalue","fdr")
GO_inc$foldchange <- as.numeric(sub("<", "", GO_inc$foldchange))
GO_inc <- GO_inc[GO_inc$fdr<.001,]
GO_inc <- GO_inc[order(-GO_inc$foldchange),]
barplot(GO_inc$foldchange[15:1], names.arg = GO_inc$Process[15:1], las=2, xlim=c(0,2.5), 
                              col = "firebrick2", main = "Molecular Functions Enriched in 4hr 100nM Romi Activated Genes                                                                     ", 
                              horiz=T, las=1, xlab = "Enrichment"); box()

GO_dec <- read.table("GO_4_100_func_dec.txt", header=T, sep="\t")
GO_dec <- GO_dec[,c(-2)]
colnames(GO_dec) <- c("Process","Count","Expected","plusminus","foldchange","pvalue","fdr")
GO_dec$foldchange <- as.numeric(sub("<", "", GO_dec$foldchange))
GO_dec <- GO_dec[GO_dec$fdr<.001,]
GO_dec <- GO_dec[order(-GO_dec$foldchange),]
barplot(GO_dec$foldchange[15:1], names.arg = GO_dec$Process[15:1], las=2, xlim=c(0,6), 
                              col = "dodgerblue", main = "Molecular Functions Enriched in 4hr 100nM Romi Repressed Genes                                                                         ", 
                              horiz=T, las=1, xlab = "Enrichment"); box();

#PATHWAYS
GO_inc <- read.table("GO_4_100_path_inc.txt", header=T, sep="\t")
GO_inc <- GO_inc[,c(-2)]
colnames(GO_inc) <- c("Process","Count","Expected","plusminus","foldchange","pvalue","fdr")
GO_inc$foldchange <- as.numeric(sub("<", "", GO_inc$foldchange))
GO_inc <- GO_inc[order(GO_inc$pvalue),]
GO_inc <- GO_inc[c(1:10),]
GO_inc <- GO_inc[order(GO_inc$foldchange),]
barplot(GO_inc$foldchange, names.arg = GO_inc$Process, las=2, xlim=c(0,4), 
        col = "firebrick2", main = "Pathways Enriched in 4hr 100nM Romi Activated Genes                                                     ", 
        horiz=T, las=1, xlab = "Enrichment"); box();

GO_dec <- read.table("GO_4_100_path_dec.txt", header=T, sep="\t")
GO_dec <- GO_dec[,c(-2)]
colnames(GO_dec) <- c("Process","Count","Expected","plusminus","foldchange","pvalue","fdr")
GO_dec$foldchange <- as.numeric(sub("<", "", GO_dec$foldchange))
GO_dec <- GO_dec[order(GO_dec$pvalue),]
GO_dec <- GO_dec[c(1:10),]
GO_dec <- GO_dec[order(GO_dec$foldchange),]
barplot(GO_dec$foldchange, names.arg = GO_dec$Process, las=2, xlim=c(0,4), 
                              col = "dodgerblue", main = "Pathways Enriched in 4hr 100nM Romi Repressed Genes                                                    ", 
                              horiz=T, las=1, xlab = "Enrichment"); box()

#Various Custom colors for heatmap, the recommended is blocol
blocol = colorRampPalette(c("blue","blue3","black","red3","red"))
redcol = colorRampPalette(c("black","red"))
whocol = colorRampPalette(c("white","red"))
grecol = colorRampPalette(c("green","green4","black","orchid4","orchid"))
orocol = colorRampPalette(c("blue","black","orange"))
hotcol = colorRampPalette(c("white","yellow","orange","red"))

#BART Analysis to identify TFs driving differential expression
bart_4_100_inc <- read.table("C:/School/UVA/Research_Main/RNAseq/BART_Results/4_100_inc_results.txt", header=TRUE)
bart_4_100_dec <- read.table("C:/School/UVA/Research_Main/RNAseq/BART_Results/4_100_dec_results.txt", header=TRUE)

bart_24_100_inc <- read.table("C:/School/UVA/Research_Main/RNAseq/BART_Results/24_100_inc_results.txt", header=TRUE)
bart_24_100_dec <- read.table("C:/School/UVA/Research_Main/RNAseq/BART_Results/24_100_dec_results.txt", header=TRUE)

bart_24_10_inc <- read.table("C:/School/UVA/Research_Main/RNAseq/BART_Results/24_10_inc_results.txt", header=TRUE)
bart_24_10_dec <- read.table("C:/School/UVA/Research_Main/RNAseq/BART_Results/24_10_dec_results.txt", header=TRUE)

### INCREASING 4hr 100nM
bart_4_100_inc <- subset(bart_4_100_inc, irwin_hall_pvalue < .01)
bart_4_100_df <- df_4_100[df_4_100$gene %in% bart_4_100_inc$TF,]
bart_4_100_inc <- bart_4_100_inc[bart_4_100_inc$TF %in% bart_4_100_df$gene,]

bart_4_100_df <- bart_4_100_df[order(bart_4_100_df$gene),]

bart_4_100_inc <- bart_4_100_inc[order(bart_4_100_inc$TF),]

bart <- cbind(bart_4_100_inc,bart_4_100_df)
bart <- bart[,c(1,7,8,9,13)]
bart <- bart[order(-bart$baseMean),]

bart_expression_4_100 <- bart[bart$TF %in% rownames(z),]
barplot(bart_expression_4_100$baseMean, names = bart_expression_4_100$TF, las = 2, ylab = "Relative Expression",
        ylim = c(0,10000),
        col=ifelse(bart_expression_4_100$padj < .001,
                   ifelse(bart_expression_4_100$log2FoldChange > .25, "darkorchid2",
                          ifelse(bart_expression_4_100$log2FoldChange < -.25, "green3", "grey80" )),"grey80"),
        main="Relative Expression of 75 Transcription Factors Predicted to Drive Upregulation in 4hr 100nM Romi");box(); abline(h=540, lty =3, col ="red")
legend("top", 
       legend = c("4hr 100nM Romi Repressed", 
                  "4hr 100nM Romi Activated",
                  "All Other Genes"), 
       fill = c("green3", "darkorchid2","grey80"));  


bart_ZNF <- bart[bart$TF %in% c("ZNF644","ZNF639","ZNF24","ZNF318","ZNF740","ZNF143","ZNF197","ZNF589","ZNF584",
                                "ZKSCAN8","ZC3H8","ZSCAN29","ZBTB40", "ZEB2","CHAMP1"),]
barplot(bart_ZNF$baseMean, names = bart_ZNF$TF, ylab = "Relative Expression", ylim = c(0,1000), las =2,
        col=ifelse(bart_ZNF$padj < .001,
                   ifelse(bart_ZNF$log2FoldChange > .2, "darkorchid2",
                          ifelse(bart_ZNF$log2FoldChange < -.2, "green3", "grey80" )) ,"grey80"),
        main="Relative Expression of Zinc Finger Family Members"); box(); abline(h=540, lty =3, col ="red")

bart_HDAC <- bart[bart$TF %in% c("TBL1XR1","HDAC8","HDAC1","HINFP","KDM5A","KAT8","SETDB1","EPC1"),]
barplot(bart_HDAC$baseMean, names = bart_HDAC$TF, ylab = "Relative Expression", ylim = c(0,6000),
        col=ifelse(bart_HDAC$padj < .001,
                   ifelse(bart_HDAC$log2FoldChange > .2, "darkorchid2",
                          ifelse(bart_HDAC$log2FoldChange < -.2, "green3", "grey80" )) ,"grey80"),
        main="Relative Expression of Histone-related Genes"); box(); abline(h=540, lty =3, col ="red")

bart_ETS <- bart[bart$TF %in% c("ELF1","ETS2","ELK1","GABPA"),]
barplot(bart_ETS$baseMean, names = bart_ETS$TF, ylab = "Relative Expression", ylim = c(0,1400),
        col=ifelse(bart_ETS$padj < .001,
                   ifelse(bart_ETS$log2FoldChange > .2, "darkorchid2",
                          ifelse(bart_ETS$log2FoldChange < -.2, "green3", "grey80" )) ,"grey80"),
        main="Relative Expression of Candidate ETS Family Members"); box();  abline(h=540, lty =3, col ="red")


### HEATMAP FOR BART PREDICTED TFs UPREGULATING 4_100
library(gplots)
df_bart_4_100_inc <- cbind(df_4_10[df_4_10$gene %in% bart_4_100_inc$TF,],
                       df_4_100[df_4_100$gene %in% bart_4_100_inc$TF,],
                       df_24_10[df_24_10$gene %in% bart_4_100_inc$TF,],
                       df_24_100[df_24_100$gene %in% bart_4_100_inc$TF,])

rownames(df_bart_4_100_inc) = df_bart_4_100_inc$gene
df_bart_4_100_inc <- df_bart_4_100_inc[,c(3,11,19,27,7,15,23,31)]
colnames(df_bart_4_100_inc) = c("4_10","4_100","24_10","24_100","4_10_p","4_100_p","24_10_p","24_100_p")
z <- df_bart_4_100_inc
z$stat_4_100 <- ifelse(z$`4_100_p` < .001, ifelse(z$`4_100` > .25, 
                                                  "darkorchid2",ifelse(z$`4_100` < -.25, "green3","gray")),"gray")

z <- z[!rownames(z) %in% c("GATA1","NFE2","ZNF175","ZNF83","TAL1"),]
z[45,9] <- "gray"
z[57,9] <- "gray"

raw_4_100_inc <- rld_HH[row.names(rld_HH) %in% row.names(z),]
raw_4_100_inc <- raw_4_100_inc[order(row.names(raw_4_100_inc))]

x <- scale(as.matrix(assay(rld_HH)))
colcondition <- c("lightgray","lightgray","lightgray","yellow","yellow","yellow","gold","gold","gold",
                  "gray50","gray50","gray50","orange","orange","orange","orange3","orange3","orange3")

pdf("4_100_bart_inc.pdf", width = 14, height = 11.5)
  par(mar=c(27,3,4.1,2.1)) 
  heatmap.2(assay(raw_4_100_inc), density.info = "none",trace = "none", breaks=300, col = blocol,
          lwid = c(5,10), margins = c(2,10),lhei=c(2,11),keysize = 1.2, scale = "row", cexRow = 1.1,
          key.xlab = "Relative Expression", symbreaks = T, 
          hclustfun = function(x) hclust(x, method = "ward.D2"),
          distfun = function(x) as.dist(1-cor(t(x))),
          ColSideColors = colcondition, labCol = F, main = "Expression of 75 Transcription Factors
Predicted to Drive Upregulation in 4hr 100nM Romi",
          RowSideColors = z$stat_4_100); legend("right",
          bty="n", inset = .795,
          legend = c("4hr 100nM Repressed","4hr 100nM Activated","",
                     "4hr DMSO","4hr 10nM","4hr 100nM","24hr DMSO","24hr 10nM","24hr 100nM"), 
          pch = 15, pt.cex=3.2, y.intersp = 1.2, x.intersp = 1.2,
          col=c("green3","darkorchid2","white","lightgray","yellow","gold","grey50","orange","orange3"))
dev.off()

### DECREASING 4hr 100nM
bart_4_100_dec <- subset(bart_4_100_dec, irwin_hall_pvalue < .01)
bart_4_100_dec_df <- df_4_100[df_4_100$gene %in% bart_4_100_dec$TF,]
bart_4_100_dec <- bart_4_100_dec[bart_4_100_dec$TF %in% bart_4_100_dec_df$gene,]

bart_4_100_dec_df <- bart_4_100_dec_df[order(bart_4_100_dec_df$gene),]

bart_4_100_dec <- bart_4_100_dec[order(bart_4_100_dec$TF),]

bart_dec <- cbind(bart_4_100_dec,bart_4_100_dec_df)
bart_dec <- bart_dec[,c(1,7,8,9,13)]
bart_dec <- bart_dec[order(-bart_dec$baseMean),]
bart_dec <- na.omit(bart_dec)

barplot((bart_dec$baseMean), names = bart_dec$TF, las = 2, ylab = "Relative Expression", cex.names=.9,
        ylim = c(0,27000),
        col=ifelse(bart_dec$padj < .001,
                   ifelse(bart_dec$log2FoldChange > .25, "darkorchid2",
                          ifelse(bart_dec$log2FoldChange < -.2, "green3", "grey80" )) ,"grey80"),
        main="Relative Expression of BART Candidate Transcription Factors"); box(); abline(h = median(bart_dec$baseMean), lty=3, col="red")
legend("top", 
                                                                            legend = c("4hr 100nM Romi Repressed", 
                                                                                       "4hr 100nM Romi Activated",
                                                                                       "All Other Genes"), 
                                                                            fill = c("green3", "darkorchid2","grey80"))

ETS_dec <- bart_dec[bart_dec$TF %in% c("ETV6","ETS1","GABPA","ELK1"),]
barplot(ETS_dec$baseMean, names = ETS_dec$TF, ylab = "Relative Expression", ylim = c(0,3000),
        col=ifelse(ETS_dec$padj < .001,
                   ifelse(ETS_dec$log2FoldChange > .25, "darkorchid2",
                          ifelse(ETS_dec$log2FoldChange < -.25, "green3", "grey80" )) ,"grey80"),
        main="Relative Expression of ETS Family Members"); box(); abline(h = median(bart_dec$baseMean), lty=3, col="red")

STAT_dec <- bart_dec[bart_dec$TF %in% c("STAT3","STAT2","STAT5","STAT5A","STAT5B"),]
barplot(STAT_dec$baseMean, names = STAT_dec$TF, ylab = "Relative Expression", ylim = c(0,1200),
        col=ifelse(STAT_dec$padj < .001,
                   ifelse(STAT_dec$log2FoldChange > .25, "darkorchid2",
                          ifelse(STAT_dec$log2FoldChange < -.25, "green3", "grey80" )) ,"grey80"),
        main="Relative Expression of STAT Family Members"); box(); abline(h = median(bart_dec$baseMean), lty=3, col="red")

NFKB_dec <- bart_dec[bart_dec$TF %in% c("NFKB1","NFKB2","RELA","RELB","REL"),]
barplot(NFKB_dec$baseMean, names = NFKB_dec$TF, ylab = "Relative Expression", ylim = c(0,17000),
        col=ifelse(NFKB_dec$padj < .001,
                   ifelse(NFKB_dec$log2FoldChange > .25, "darkorchid2",
                          ifelse(NFKB_dec$log2FoldChange < -.25, "green3", "grey80" )) ,"grey80"),
        main="Relative Expression of NFKB Family Members"); box(); abline(h = median(bart_dec$baseMean), lty=3, col="red")

### HEATMAP FOR BART PREDICTED TFs DOWNREGULATING 4_100                                                                             
df_bart_4_100_dec <- cbind(df_4_10[df_4_10$gene %in% bart_4_100_dec$TF,],
                           df_4_100[df_4_100$gene %in% bart_4_100_dec$TF,],
                           df_24_10[df_24_10$gene %in% bart_4_100_dec$TF,],
                           df_24_100[df_24_100$gene %in% bart_4_100_dec$TF,])



        rownames(df_bart_4_100_dec) = df_bart_4_100_dec$gene
        df_bart_4_100_dec <- df_bart_4_100_dec[,c(3,11,19,27,7,15,23,31)]
        colnames(df_bart_4_100_dec) = c("4_10","4_100","24_10","24_100","4_10_p","4_100_p","24_10_p","24_100_p")
        m <- df_bart_4_100_dec
        m$stat_4_100 <- ifelse(m$`4_100_p` < .001, ifelse(m$`4_100` > .25, 
                                                          "darkorchid2",ifelse(m$`4_100` < -.25, "green3","gray")),"gray")

        m <- m[!rownames(m) %in% c("CDK7","EBF1","MAF","NFATC1","NFE2","SPI1","TBX21","ZBED1"),]

raw_4_100_dec <- rld_HH[row.names(rld_HH) %in% row.names(m),]
raw_4_100_dec <- raw_4_100_dec[order(row.names(raw_4_100_dec))]
pdf("4_100_bart_dec.pdf", width = 14, height = 11.5)
  par(mar=c(27,3,4.1,2.1)) 
  heatmap.2(assay(raw_4_100_dec), density.info = "none",trace = "none", breaks=300, col = blocol,
            lwid = c(5,10), margins = c(2,10),lhei=c(2,11),keysize = 1.2, scale = "row", cexRow = 1.1,
            key.xlab = "Relative Expression", symbreaks = T, 
            hclustfun = function(x) hclust(x, method = "ward.D2"),
            distfun = function(x) as.dist(1-cor(t(x))),
            ColSideColors = colcondition, labCol = F, main = "Expression of 66 Transcription Factors
            Predicted to Drive Downregulation in 4hr 100nM Romi",
            RowSideColors = m$stat_4_100); legend("right",
          bty="n", inset = .795,
          legend = c("4hr 100nM Repressed","4hr 100nM Activated","",
                     "4hr DMSO","4hr 10nM","4hr 100nM","24hr DMSO","24hr 10nM","24hr 100nM"), 
          pch = 15, pt.cex=3.2, y.intersp = 1.2, x.intersp = 1.2,
          col=c("green3","darkorchid2","white","lightgray","yellow","gold","grey50","orange","orange3"))
  dev.off()


### HEATMAP FOR BART PREDICTED TFS SHARED BETWEEN UP AND DOWNREGULATION
gamma <- z[row.names(z) %in% row.names(m),]
gamma <- rbind(m[row.names(m) %in% row.names(z),])

delta <- rld_HH[row.names(rld_HH) %in% row.names(gamma),]
delta <- delta[order(row.names(delta))]
pdf("4_100_bart_updown.pdf",width = 12, height = 8)
par(mar=c(20,2,4.1,2.1)) 
heatmap.2(assay(delta), density.info = "none",trace = "none", breaks=300, col = blocol,
          lwid = c(5,10), margins = c(2,10),lhei=c(2,11),keysize = 1.2, scale = "row", cexRow = 1.1,
          key.xlab = "Relative Expression", symbreaks = T, 
          hclustfun = function(x) hclust(x, method = "ward.D"),
          distfun = function(x) as.dist(1-cor(t(x))),
          ColSideColors = colcondition, labCol = F, main = "Expression of 11 Transcription Factors Predicted to Drive
Both Upregulation and Downregulation in 4hr 100nM Romi",
          RowSideColors = gamma$stat_4_100); legend("right",
          bty="n", inset = .795,
          legend = c("4hr 100nM Repressed","4hr 100nM Activated","",
                     "4hr DMSO","4hr 10nM","4hr 100nM","24hr DMSO","24hr 10nM","24hr 100nM"), 
          pch = 15, pt.cex=2.7, y.intersp = 1.1, x.intersp = 1.2, cex=.9,
          col=c("green3","darkorchid2","white","lightgray","yellow","gold","grey50","orange","orange3"))
dev.off()

theta <- df_4_100[df_4_100$gene %in% row.names(gamma),]
theta <- theta[order(-theta$baseMean),]
barplot(theta$baseMean, names = theta$gene, ylab = "Relative Expression", ylim = c(0,6000),
        col=ifelse(theta$padj < .001,
                   ifelse(theta$log2FoldChange > .25, "darkorchid2",
                          ifelse(theta$log2FoldChange < -.25, "green3", "grey80" )) ,"grey80"),
        main="Relative Expression of 11 TFs Predicted to Drive both Up and Downregulation at 4hr 100nM"); box();

### CDF AND DENSITY FUNCTIONS
library(dgof)

cdf <- df_4_100
cdf_inc <- cdf[cdf$status == "increased",]
cdf_dec <- cdf[cdf$status == "decreased",]
cdf_unc <- na.omit(cdf[cdf$status == "unchanged",])
cdf_total <- rbind(cdf_inc, cdf_dec, cdf_unc)

pdf("CDF_4_100.pdf", width = 7, height = 7)
plot(ecdf(log10(cdf_inc$baseMean)), col ="red", xlab = "Log10 Normalized Count", ylab = "Cumulative Distribution Function",
     main = "Expression CDF of All Genes at 4hr 100nM Romidepsin",lwd = 2); lines(ecdf(log10(cdf_dec$baseMean)), 
        col ="blue",lwd = 2); lines(ecdf(log10(cdf_unc$baseMean)), col = "gray",lwd = 2); legend("right", bty = "n",
          legend = c("Activated","Repressed","Unchanged"), col = c("red","blue","gray"), lwd = 2, cex = .9, inset = .05)
dev.off()

pdf("CDF_4_100_change.pdf", width = 7, height = 7)
plot(ecdf(abs(cdf_inc$log2FoldChange)), col ="red", xlab = "Absolute Log2 Fold Change", ylab = "Cumulative Distribution Function", xlim = c(0,5),
     main = "Fold Change CDF of All Genes at 4hr 100nM Romidepsin",lwd = 2); lines(ecdf(abs(cdf_dec$log2FoldChange)), 
     col ="blue",lwd = 2); lines(ecdf(abs(cdf_unc$log2FoldChange)), col = "gray",lwd = 2); legend("right", bty = "n",
     legend = c("Activated","Repressed","Unchanged"), col = c("red","blue","gray"), lwd = 2, cex = .9, inset = .05)
dev.off()

ks.test(cdf_inc$baseMean, cdf_dec$baseMean, alternative = "greater")

pdf("Density_4_100.pdf", width = 7, height = 7)
plot(density(log10(cdf_inc$baseMean)), col ="red", xlab = "Log10 Normalized Count", ylab = "Density",
     main = "Expression Density of All Genes at 4hr 100nM Romidepsin", lwd = 2); lines(density(log10(cdf_dec$baseMean)), 
        col ="blue",lwd = 2); lines(density(log10(cdf_unc$baseMean)), col = "gray",lwd = 2); legend("topright", bty = "n",
          legend = c("Activated","Repressed","Unchanged"), col = c("red","blue","gray"), lwd = 2, cex = .9, inset = .05)
dev.off()

pdf("Change_4_100.pdf", width = 7, height = 7)
plot(density(abs(cdf_inc$log2FoldChange)), col ="red", xlab = "Absolute Log2 Fold Change", ylab = "Density", xlim = c(0,5), ylim = c(0, 2),
     main = "Fold Change Density of All Genes at 4hr 100nM Romidepsin",lwd = 2); lines(density(abs(cdf_dec$log2FoldChange)), 
  col ="blue",lwd = 2); lines(density(abs(cdf_unc$log2FoldChange)), col = "gray",lwd = 2); legend("topright", bty = "n",
legend = c("Activated","Repressed","Unchanged"), col = c("red","blue","gray"), lwd = 2, cex = .9, inset = .05)
dev.off()

ks.test(abs(cdf_inc$log2FoldChange), abs(cdf_dec$log2FoldChange), alternative = "less")

#BART CDF INCREASED
cdf_z <- cdf[cdf$gene %in% row.names(z),]
cdf_z_inc <- cdf_z[cdf_z$status == "increased",]
cdf_z_dec <- cdf_z[cdf_z$status == "decreased",]
cdf_z_unc <- cdf_z[cdf_z$status == "unchanged",]

    pdf("Density_BART_inc.pdf", width = 7, height = 7)
    plot(density(log10(cdf_z_inc$baseMean)), col ="darkorchid2", xlab = "Log10 Normalized Count", ylab = "Density", xlim = c(0,5),
         main = "Expression Density of 75 TFs 
    Predicted to Drive Upregulation in 4hr 100nM Romi",lwd = 2); lines(density(log10(cdf_z_dec$baseMean)), 
            col ="green3",lwd = 2); lines(density(log10(cdf_z_unc$baseMean)), col = "gray",lwd = 2); legend("topleft", bty = "n",
             legend = c("Activated","Repressed","Unchanged"), col = c("purple","green","gray"), lwd = 2, cex = .9, inset = .05)
    dev.off()

ks.test(cdf_z_inc$baseMean, cdf_z_dec$baseMean, alternative = "greater") 

    
    pdf("Change_BART_inc.pdf", width = 7, height = 7)
    plot(density(abs(cdf_z_inc$log2FoldChange)), col ="darkorchid2", xlab = "Absolute Log2 Fold Change", ylab = "Density", xlim = c(0,4), ylim = c(0, 3),
         main = "Fold Change Density of 75 TFs 
    Predicted to Drive Upregulation in 4hr 100nM Romi",lwd = 2); lines(density(abs(cdf_z_dec$log2FoldChange)), 
          col ="green3",lwd = 2); lines(density(abs(cdf_z_unc$log2FoldChange)), col = "gray",lwd = 2); legend("topright", bty = "n",
          legend = c("Activated","Repressed","Unchanged"), col = c("purple","green","gray"), lwd = 2, cex = .9, inset = .05)
    dev.off()

ks.test(abs(cdf_z_inc$log2FoldChange), abs(cdf_z_dec$log2FoldChange), alternative = "less")  

#BART CDF DECREASED
cdf_m <- cdf[cdf$gene %in% row.names(m),]
cdf_m_inc <- cdf_m[cdf_m$status == "increased",]
cdf_m_dec <- cdf_m[cdf_m$status == "decreased",]
cdf_m_unc <- cdf_m[cdf_m$status == "unchanged",]

    pdf("Density_BART_dec.pdf", width = 7, height = 7)
    plot(density(log10(cdf_m_inc$baseMean)), col ="darkorchid2", xlab = "Log10 Normalized Count", ylab = "Density", xlim = c(0,5),
         main = "Expression Density of 66 TFs 
    Predicted to Drive Downregulation in 4hr 100nM Romi",lwd = 2); lines(density(log10(cdf_m_dec$baseMean)), 
          col ="green3",lwd = 2); lines(density(log10(cdf_m_unc$baseMean)), col = "gray",lwd = 2); legend("topleft", bty = "n",
           legend = c("Activated","Repressed","Unchanged"), col = c("purple","green","gray"), lwd = 2, cex = .9, inset = .05)
    dev.off()

ks.test(cdf_m_inc$baseMean, cdf_m_dec$baseMean, alternative = "greater")    
    
    pdf("Change_BART_dec.pdf", width = 7, height = 7)
    plot(density(abs(cdf_m_inc$log2FoldChange)), col ="darkorchid2", xlab = "Absolute Log2 Fold Change", ylab = "Density", xlim = c(0,3), ylim = c(0, 2.7),
         main = "Fold Change Density of 66 TFs 
    Predicted to Drive Downregulation in 4hr 100nM Romi",lwd = 2); lines(density(abs(cdf_m_dec$log2FoldChange)), 
         col ="green3",lwd = 2); lines(density(abs(cdf_m_unc$log2FoldChange)), col = "gray",lwd = 2); legend("topright", bty = "n",
          legend = c("Activated","Repressed","Unchanged"), col = c("purple","green","gray"), lwd = 2, cex = .9, inset = .05)
    dev.off()

ks.test(abs(cdf_m_inc$log2FoldChange), abs(cdf_m_dec$log2FoldChange), alternative = "less")


###ROTATING ANIMATION OF 4hr 100nM Data (p-value, basemean, foldchange)
library(lattice)
library(animation)
cloud(-log(padj) ~ log10(baseMean)*log2FoldChange, data = cdf_total, group = status, col = c("blue", "red", "gray"),
      xlab = "", ylab = "", zlab = "", screen=list(z=50,x=-60), scales=list(draw=FALSE))

angles <- seq(from=1,to=360,by=.5)
draw.plot <- function(angles)
{
  for(i in 1:length(angles))
  {
    print(cloud(-log(padj) ~ log10(baseMean)*log2FoldChange, data = cdf_total, group = status, col = c("blue", "red", "gray"),
                xlab = "", ylab = "", zlab = "", screen=list(z=angles[i],x=-60),scales=list(draw=FALSE)))
    setTxtProgressBar(txtProgressBar(style=3),i/length(angles))
  }
}

saveGIF(draw.plot(angles), interval = 1/30, movie.name="whatisthis.GIF",
        ani.height=640,ani.width=640)

