#NoAb_Sham - light blue
#Ab_Sham - dark grey
#NoAB_TBI - blue
#Ab_TBI - orange

library(lattice)

setwd("~/Desktop/Ashley_RNAseq")

qc <- read.table('qc_metrics.txt',sep='\t',header=TRUE) 

rownames(qc) = c('Ab Sham 1','Ab Sham 2','Ab Sham 3','Ab Sham 4',
                 'Ab TBI 1','Ab TBI 2','Ab TBI 3','Ab TBI 4',
                 'noAb Sham 1','noAb Sham 2','noAb Sham 3','noAb Sham 4',
                 'noAb TBI 1','noAb TBI 2','noAb TBI 3','noAb TBI 4')
qc = qc[,-1]

qc$prop.junk = 1 - (qc$without.junk / qc$all.counts)
qc$prop.dups = 1 - (qc$unique.counts / qc$without.junk)

load('all.counts.Rdata')
qc$prop.in.features = colSums(all.counts[1:55421,]) / colSums(all.counts)

#qc$order = 1:25

qc$colors = c(rep('dark grey',4),rep('orange',4),rep('light blue',4),rep('blue',4))

pdf('total.reads.pdf', width=10, height=4) 
print(barchart(all.counts ~ rownames(qc),  data = qc,
               main = "Total Reads",
               xlab = "Sample",
               ylab = "Read Count",
               col = qc$colors,
               ylim = c(0,100000000),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()

pdf('proportion.junk.pdf', width=10, height=4) 
print(barchart(prop.junk ~ rownames(qc),  data = qc,
               main = "Proportion Junk",
               xlab = "Sample",
               ylab = "Proportion Junk",
               col = qc$colors,
               ylim = c(0,1),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()

pdf('proportion.duplicates.pdf', width=10, height=4) 
print(barchart(prop.dups ~ rownames(qc),  data = qc,
               main = "Proportion PCR Duplicates",
               xlab = "Sample",
               ylab = "Proportion PCR Duplicates",
               col = qc$colors,
               ylim = c(0,1),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()

pdf('proportion.in.features.pdf', width=10, height=4) 
print(barchart(prop.in.features ~ rownames(qc),  data = qc,
               main = "Proportion In Features",
               xlab = "Sample",
               ylab = "Proportion In Features",
               col = qc$colors,
               ylim = c(0,1),
               scales = list(x = list(rot = 60)),
               horiz = FALSE))
dev.off()
