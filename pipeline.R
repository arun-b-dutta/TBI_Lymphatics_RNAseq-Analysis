library(lattice)
library(DESeq2)
library(pheatmap)
library(GSA)
library(seq2pathway)

source('https://raw.githubusercontent.com/guertinlab/seqOutBias/master/docs/R/seqOutBias_functions.R')
source('https://raw.githubusercontent.com/mjg54/znf143_pro_seq_analysis/master/docs/ZNF143_functions.R')

setwd('RNAseq_analysis')

categorize.deseq.df <- function(df, fdr = 0.05, log2fold = 0.0, treat
= 'Auxin') {

     df.activated = data.frame(matrix(nrow = 0, ncol = 0))
     df.repressed = data.frame(matrix(nrow = 0, ncol = 0))

     if (nrow(df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold,]) != 0) {
     	df.activated = df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold,]
	df.activated$response = paste(treat, 'Activated')
	}

     if (nrow(df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold,]) != 0) {
     	df.repressed = df[df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold,]
	df.repressed$response = paste(treat, 'Repressed')
	}
    
    df.unchanged = df[df$padj > 0.5 & !is.na(df$padj) & abs(df$log2FoldChange) < 0.25,]
    df.unchanged$response = paste(treat, 'Unchanged')

    df.dregs = df[!(df$padj < fdr & !is.na(df$padj) & df$log2FoldChange > log2fold) &
                  !(df$padj < fdr & !is.na(df$padj) & df$log2FoldChange < -log2fold) &
                  !(df$padj > 0.5 & !is.na(df$padj) &
    		  abs(df$log2FoldChange) < 0.25), ]
    df.dregs$response = paste(treat, 'All Other Genes')
    
    df.effects.lattice = 
    rbind(df.activated, 
          df.unchanged, 
          df.repressed, 
          df.dregs)
    df.effects.lattice$response = factor(df.effects.lattice$response)
    df.effects.lattice$response = relevel(df.effects.lattice$response, ref = paste(treat, 'Unchanged'))
    df.effects.lattice$response = relevel(df.effects.lattice$response, ref = paste(treat, 'All Other Genes'))
    return(df.effects.lattice)
}

plotPCAlattice <- function(df, file = 'PCA_lattice.pdf') {  
  perVar = round(100 * attr(df, "percentVar"))
  df = data.frame(cbind(df, sapply(strsplit(as.character(df$name), '_rep'), '[', 1)))
  colnames(df) = c(colnames(df)[1:(ncol(df)-1)], 'unique_condition')
  print(df)
  #get colors and take away the hex transparency
  color.x = substring(rainbow(length(unique(df$unique_condition))), 1,7) 
  
  df$color = NA
  df$alpha.x = NA
  df$alpha.y = NA
  df$colpal = NA
  
  for (i in 1:length(unique(df$unique_condition))) {
    
    df[df$unique_condition == unique(df$unique_condition)[[i]],]$color = color.x[i]   
    #gives replicates for unique condition
    reps_col<- df[df$unique_condition == unique(df$unique_condition)[[i]],]
    #gives number of replicates in unique condition
    replicates.x = nrow(reps_col)
    alx <- rev(seq(0.2, 1, length.out = replicates.x))
    
    #count transparency(alx), convert alx to hex(aly), combain color and transparency(cp)
    for(rep in 1:replicates.x) {
    
      na <- reps_col[rep, ]$name
      df[df$name == na, ]$alpha.x = alx[rep]
      aly = as.hexmode(round(alx * 255))
      df[df$name == na, ]$alpha.y = aly[rep]
      cp = paste0(color.x[i], aly)
      df[df$name == na, ]$colpal = cp[rep]
      #print(df)
    }
  }
  colpal = df$colpal
  df$name = gsub('_', ' ', df$name)
  df$name <- factor(df$name, levels=df$name, order=TRUE)
  pdf(file, width=6, height=6, useDingbats=FALSE)
  print(xyplot(PC2 ~ PC1, groups = name, data=df,
               xlab = paste('PC1: ', perVar[1], '% variance', sep = ''),
               ylab = paste('PC2: ', perVar[2], '% variance', sep = ''),
               par.settings = list(superpose.symbol = list(pch = c(20), col=colpal)),
               pch = 20, cex = 1.7,
               auto.key = TRUE,
               col = colpal))
  dev.off()
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
unt=4
trt=12
sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))

deseq.counts.table = DESeqDataSetFromMatrix(merged.counts, DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)
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
             groups=ab.sham.v.ab.tbi.lattice$response,
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
             groups=noAb.tbi.v.ab.tbi.lattice$response,
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
             groups=noAb.sham.v.noAb.tbi.lattice$response,
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
             groups=noAb.sham.v.Ab.sham.lattice$response,
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

#GSEA

#command line
#wget http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/7.0/msigdb.v7.0.symbols.gmt

msigdb = GSA.read.gmt('msigdb.v7.0.symbols.gmt')

lattice = na.omit(noAb.tbi.v.ab.tbi.lattice)
lattice = lattice[lattice$pvalue < 0.05 & lattice$log2FoldChange > 0,]

inc.pval.0.05.gsea = FisherTest_MsigDB(msigdb,toupper(rownames(lattice)),'mm10')

x = data.frame(GeneSet = inc.pval.0.05.gsea$GeneSet, padj = inc.pval.0.05.gsea$FDR)
write.table(x,file='inc.gene.sets.pval.0.05.txt',sep='\t')

lattice = na.omit(noAb.tbi.v.ab.tbi.lattice)
lattice = lattice[lattice$pvalue < 0.05 & lattice$log2FoldChange < 0,]

dec.pval.0.05.gsea = FisherTest_MsigDB(msigdb,toupper(rownames(lattice)),'mm10')

x = data.frame(GeneSet = dec.pval.0.05.gsea$GeneSet, padj = dec.pval.0.05.gsea$FDR)
write.table(x,file='dec.gene.sets.pval.0.05.txt',sep='\t')

gene.sets = c("GO_REGULATION_OF_IMMUNE_SYSTEM_PROCESS","GO_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
              "GO_REGULATION_OF_IMMUNE_RESPONSE","GO_INNATE_IMMUNE_RESPONSE",
              "GO_REGULATION_OF_IMMUNE_EFFECTOR_PROCESS","GO_IMMUNE_EFFECTOR_PROCESS",
              "GO_REGULATION_OF_INNATE_IMMUNE_RESPONSE","REACTOME_INNATE_IMMUNE_SYSTEM",
              "HALLMARK_COMPLEMENT","GO_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")

row.nums = c()
for (i in gene.sets) {
    row.nums = append(row.nums,grep(i,inc.pval.0.05.gsea$GeneSet))
}

x = inc.pval.0.05.gsea[row.nums,]
x = x[order(x$FDR),]
x$order = 10:1

pdf(file = 'inc.gsea.0.05.pval.pdf',width=9)
print(barchart((-log2(FDR)) ~ reorder(GeneSet,order), data = x,
               main = "Enriched Gene Sets in Increased Genes (0.05 p-val)",
               xlab = "Gene Set",
               ylab = "-log2(padj)",
               col = 'red',
               #ylim = c(0,1),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()

#heat maps

unt=4
trt=4
sample.conditions = factor(c(rep("untreated",unt), rep("treated",trt)), levels=c("untreated","treated"))        
 
deseq.counts.table = DESeqDataSetFromMatrix(merged.counts[,c(13:16,5:8)], DataFrame(sample.conditions), ~ sample.conditions);
colData(deseq.counts.table)$condition<-factor(colData(deseq.counts.table)$sample.conditions, levels=c('untreated','treated'));
dds = DESeq(deseq.counts.table)
rld_HH = rlogTransformation(dds)

a = assay(rld_HH)[c(inc[1:20],dec[1:20]),]
a = a - rowMeans(a)

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

pdf(file='inflammation.top40.heatmap.pdf')
print(pheatmap(a, cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE))
dev.off()

#DAM genes
x = read.table('DAM.genes.txt',sep='\t')
x = as.vector(x[,1])
x[18] = 'Lilrb4a'
x = x[-c(7,19)]

#pull the row numbers corresponding to complement genes and put them in a list object called 'genes'
genes = c()
for (i in x) {
    genes = append(genes,grep(paste0('^',toupper(i),'$'),toupper(rownames(noAb.tbi.v.ab.tbi.lattice))))
}

#pull the appropriate rows from the lattice object and reorder them by padj
lattice = noAb.tbi.v.ab.tbi.lattice[genes,]
lattice = lattice[order(lattice$padj),]

#define object 'a' as the normalized, log2-transformed read counts of just the genes and samples of interest
a = assay(rld_HH)[rownames(lattice),]
#center data by subtracting mean of each row
a = a - rowMeans(a)

#plot and save heat map
pdf(file='DAM.heatmap.pdf')
print(pheatmap(a, cluster_rows=FALSE, show_rownames=TRUE,cluster_cols=FALSE))
dev.off()

#volcano plots (p-value < 0.05)
volcano.plot <- function(lattice, thresh = 0.05, treat='Ablation', title="", highlights = c()) {
    name = strsplit(deparse(substitute(lattice)),'.lattice')[[1]][1]

    #for each plot we set the variable 'lattice' to be the comparison of interest and add a 'color' column 
    #that will dictate data point color
    lattice$color = 'black'
    lattice[lattice$treat == paste0(treat,' Activated'),]$color = 'red'
    lattice[lattice$treat == paste0(treat,' Repressed'),]$color = 'blue'

    #if a vector of 'highlight' genes is inputted, then ...
    if (length(highlights) > 0) {
        name = paste0(name,'.highlights')
        highlights = paste0('^',highlights,'$')
        row.nums = c()
        for(i in highlights) { 
            row.nums = append(row.nums, grep(i,rownames(lattice)))
        }
        lattice[row.nums,]$color = 'green'
    }

    pdf(file = paste0(name,'.volcano.plot.pdf'))
    
    #plot the log2FoldChange on the x-axis and the -log10(pvalue) on the y-axis
    plot(lattice$log2FoldChange, -log10(lattice$pvalue), xlim = c(-5,6), col = lattice$color, cex = 0.8, pch = 20, 
     main = title, xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
    #add a line to signify the significance threshold
    abline(h=-log(thresh,base=10), lty = 5, lwd = 4, col = '#C0C0C0')

    dev.off()
}

volcano.plot(noAb.tbi.v.ab.tbi.lattice, treat='Ablation', title='TBI No Ablation v. Ablation')
volcano.plot(ab.sham.v.ab.tbi.lattice, treat='TBI', title='Ablation Sham v. TBI')
volcano.plot(noAb.sham.v.noAb.tbi.lattice, treat ='TBI', title='No ablation Sham v. TBI')

x = c('Arc', 'C1qc', 'C1qa', 'Itgam', 'C4b', 'Fn1')
volcano.plot(noAb.tbi.v.ab.tbi.lattice, treat='Ablation', title='TBI No Ablation v. Ablation', highlights = x)

#old way - DON'T USE
lattice = noAb.sham.v.noAb.tbi.lattice

pdf(file = 'noAb.sham.v.no.Ab.tbi.volcano.plot.pdf')
plot(lattice$log2FoldChange, -log10(lattice$pvalue), cex = 0.8, pch = 20, main = 'No Ablation Sham v. TBI', xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
with(subset(lattice, pvalue < 0.05 & log2FoldChange < 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='blue'))
with(subset(lattice, pvalue < 0.05 & log2FoldChange > 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='red'))
abline(h=-log(0.05,base=10), lty = 5, lwd = 4, col = '#C0C0C0')
dev.off()

lattice = ab.sham.v.ab.tbi.lattice

pdf(file = 'ab.sham.v.ab.tbi.volcano.plot.pdf')
plot(lattice$log2FoldChange, -log10(lattice$pvalue), cex = 0.8, pch = 20, main = 'Ablation Sham v. TBI', xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
with(subset(lattice, pvalue < 0.05 & log2FoldChange < 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='blue'))
with(subset(lattice, pvalue < 0.05 & log2FoldChange > 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='red'))
abline(h=-log(0.05,base=10), lty = 5, lwd = 4, col = '#C0C0C0')
dev.off()

lattice = noAb.tbi.v.ab.tbi.lattice

pdf(file = 'noAb.tbi.v.ab.tbi.volcano.plot')
plot(lattice$log2FoldChange, -log10(lattice$pvalue), cex = 0.8, pch = 20, main = 'TBI No Ablation v. Ablation', xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
with(subset(lattice, pvalue < 0.05 & log2FoldChange < 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='blue'))
with(subset(lattice, pvalue < 0.05 & log2FoldChange > 0),points(log2FoldChange,-log10(pvalue), cex = 0.8, pch = 20, col='red'))
abline(h=-log(0.05,base=10), lty = 5, lwd = 4, col = '#C0C0C0')
dev.off()

genes = c('^Arc$', '^C1qc$', '^C1qa$', '^Itgam$', '^C4b$', '^Fn1$')

row.nums = c()
for(i in genes) { 
row.nums = append(row.nums, grep(i,rownames(lattice)))
}

lattice$color = 'black'
lattice[lattice$treat == 'Ablation Activated',]$color = 'red'
lattice[lattice$treat == 'Ablation Repressed',]$color = 'blue'
lattice[row.nums,]$color = 'green'

pdf(file = 'noAb.tbi.v.ab.tbi.highlights.volcano.plot.pdf')

plot(lattice$log2FoldChange, -log10(lattice$pvalue), col =
lattice$color, cex = 0.8, pch = 20, main = 'TBI No Ablation v. Ablation', xlab = 'log2FoldChange', ylab = '-log10(pvalue)')
abline(h=-log(0.05,base=10), lty = 5, lwd = 4, col = '#C0C0C0')

dev.off()
