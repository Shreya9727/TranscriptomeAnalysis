setwd('path to working directory')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library("DESeq2")
library("ggplot2")


file = read.table("emt_brca.tsv", sep = '\t',header = TRUE, row.names=1, check.names =TRUE)
dim(file)
head(file)

file2 = round(2**file)-1
dim(file2)
countdata = as.matrix(file2)
head(colnames(countdata))

condition = factor(c(rep("normal",292), rep("tumor",1097)))
coldata = data.frame(row.names = colnames(countdata), condition)
coldata

ddsFull = DESeqDataSetFromMatrix(countData = countdata, colData = coldata, 
                                 design =~condition)
ddsFull

dds = DESeq(ddsFull)
res <- results( dds )
summary(res)

res_ordered = res[order(res$padj),]
res_d = as.data.frame(res_ordered)
res_d
res_d$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_d$diffexpressed[res_d$log2FoldChange > 0.6 & res_d$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_d$diffexpressed[res_d$log2FoldChange < -0.6 & res_d$pvalue < 0.05] <- "DOWN"

res_d <- cbind(rownames(res_d), data.frame(res_d, row.names=NULL))
res_d
colnames(res_d)[1] <- "genes"
res_d
plot1 = ggplot(res_d, aes(x = log2FoldChange, y = -log10(padj)))+ 
  geom_point(aes(colour=diffexpressed),size=2, alpha=1) + 
  scale_colour_manual(values=c("red","blue","green"))+
  labs(title = "Breast cancer - Differential expression")+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8)
plot1
#plot1+geom_text(data=(res_d), aes(label=genes)) 
a = res_d[which(res_d$diffexpressed=='UP'),]
b = res_d[which(res_d$diffexpressed=='DOWN'),]
write.table(a, file = 'emt_up.csv', sep = '\t', row.names=FALSE, col.names=TRUE)
write.table(b, file = 'emt_down.csv', sep = '\t', row.names=FALSE, col.names=TRUE)



