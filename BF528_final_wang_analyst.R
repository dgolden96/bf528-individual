library('org.Mm.eg.db')
library('biomaRt')
library('dplyr')
library('tidyverse')
library('limma')
library('edgeR')
library('ggplot2')
library('ggrepel')

expr_mat <- read.table("C:/Users/dango/Documents/liver-normalization-rma.txt", sep = "\t", row.names = 1, header = TRUE, as.is = TRUE)
rma <- as.data.frame(expr_mat)
samples <- read.csv("C:/Users/dango/Documents/group_6_mic_info.csv", as.is = TRUE)

rma.subset1 <- rma[paste0('X',samples$array_id[samples$chemical=='FLUCONAZOLE' | (samples$chemical=='Control' & samples$vehicle=='CORN_OIL_100_%')])]
sample_filter1 <- dplyr::filter(samples, (samples$chemical %in% c("FLUCONAZOLE","Control") & samples$vehicle == "CORN_OIL_100_%"))
design1 <- model.matrix(
  ~factor(
    sample_filter1$chemical,
    levels=c('Control','FLUCONAZOLE')
  )
)
colnames(design1) <- c('Intercept','FLUCONAZOLE')
fit1 <- lmFit(rma.subset1, design1)
fit1 <- eBayes(fit1)
t1 <- topTable(fit1, coef=2, n=nrow(rma.subset1), adjust='BH')
t1_filt <- dplyr::filter(t1, adj.P.Val < 0.05)
t1_dblfilt <- dplyr::filter(t1_filt, abs(logFC) > 1.5)
top10_1_dblfilt <- dplyr::arrange(t1_dblfilt, P.Value)[1:10,]
top10_1 <- dplyr::arrange(t1_filt, P.Value)[1:10,]
#write.csv(top10_1, "fluconazole_top10.csv")
write.csv(top10_1_dblfilt, "fluco_dblfilt.csv")

rma.subset2 <- rma[paste0('X',samples$array_id[samples$chemical=='3-METHYLCHOLANTHRENE' | (samples$chemical=='Control' & samples$vehicle=='CMC_.5_%')])]
sample_filter2 <- dplyr::filter(samples, (samples$chemical %in% c("3-METHYLCHOLANTHRENE","Control") & samples$vehicle == "CMC_.5_%"))
design2 <- model.matrix(
  ~factor(
    sample_filter2$chemical,
    levels=c('Control','3-METHYLCHOLANTHRENE')
  )
)
colnames(design2) <- c('Intercept','3-METHYLCHOLANTHRENE')
fit2 <- lmFit(rma.subset2, design2)
fit2 <- eBayes(fit2)
t2 <- topTable(fit2, coef=2, n=nrow(rma.subset2), adjust='BH')
t2_filt <- dplyr::filter(t2, adj.P.Val < 0.05)
t2_dblfilt <- dplyr::filter(t2_filt, abs(logFC) > 1.5)
top10_2_dblfilt <- dplyr::arrange(t2_dblfilt, P.Value)[1:10,]
top10_2 <- dplyr::arrange(t2_filt, P.Value)[1:10,]
#write.csv(top10_2, "3_methylcholanthrene_top10.csv")
write.csv(top10_2_dblfilt, "meth_dblfilt.csv")

rma.subset3 <- rma[paste0('X',samples$array_id[samples$chemical=='PIRINIXIC_ACID' | (samples$chemical=='Control' & samples$vehicle=='CMC_.5_%')])]
sample_filter3 <- dplyr::filter(samples, (samples$chemical %in% c("PIRINIXIC_ACID","Control") & samples$vehicle == "CMC_.5_%"))
design3 <- model.matrix(
  ~factor(
    sample_filter3$chemical,
    levels=c('Control','PIRINIXIC_ACID')
  )
)
colnames(design3) <- c('Intercept','PIRINIXIC_ACID')
fit3 <- lmFit(rma.subset3, design3)
fit3 <- eBayes(fit3)
t3 <- topTable(fit3, coef=2, n=nrow(rma.subset3), adjust='BH')
t3_filt <- dplyr::filter(t3, adj.P.Val < 0.05)
t3_dblfilt <- dplyr::filter(t3_filt, abs(logFC) > 1.5)
top10_3_dblfilt <- dplyr::arrange(t3_dblfilt, P.Value)[1:10,]
top10_3 <- dplyr::arrange(t3_filt, P.Value)[1:10,]
#write.csv(top10_3, "pirinixic_acid.csv")
write.csv(top10_3_dblfilt, "pirin_dblfilt.csv")

hist1 <- ggplot(t1_filt, aes(x=logFC)) +
  geom_histogram(binwidth=0.25, col='white') +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Count') +
  ggtitle('Differential Gene Expression for Fluconazole Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(hist1)
#ggsave(filename='fluconazole_hist.png', plot=hist1)

hist1dbl <- ggplot(t1_dblfilt, aes(x=logFC)) +
  geom_histogram(binwidth=0.25, col='white') +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Count') +
  ggtitle('Differential Gene Expression for Fluconazole Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(hist1dbl)
#ggsave(filename='fluconazole_histdbl.png', plot=hist1dbl)

hist2 <- ggplot(t2_filt, aes(x=logFC)) +
  geom_histogram(binwidth=0.25, col='white') +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Count') +
  ggtitle('Differential Gene Expression for 3-Methylcholanthrene Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(hist2)
#ggsave(filename='3_methylcholanthrene_hist.png', plot=hist2)

hist2dbl <- ggplot(t2_dblfilt, aes(x=logFC)) +
  geom_histogram(binwidth=0.25, col='white') +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Count') +
  ggtitle('Differential Gene Expression for 3-Methylcholanthrene Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(hist2dbl)
#ggsave(filename='3_methylcholanthrene_histdbl.png', plot=hist2dbl)

hist3 <- ggplot(t3_filt, aes(x=logFC)) +
  geom_histogram(binwidth=0.25, col='white') +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Count') +
  ggtitle('Differential Gene Expression for Pirinixic Acid Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(hist3)
#ggsave(filename='pirinixic_acid_hist.png', plot=hist3)

hist3dbl <- ggplot(t3_dblfilt, aes(x=logFC)) +
  geom_histogram(binwidth=0.25, col='white') +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Count') +
  ggtitle('Differential Gene Expression for Pirinixic Acid Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(hist3dbl)
#ggsave(filename='pirinixic_acid_histdbl.png', plot=hist3dbl)

scatter1 <- ggplot(t1_filt, aes(x=logFC,y=P.Value)) +
  geom_point() +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Nominal p-value') +
  ggtitle('Log2FC vs. p-value from DE Analysis of Fluconazole Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(scatter1)
#ggsave(filename='fluconazole_scatter.png', plot=scatter1)

scatter1dbl <- ggplot(t1_dblfilt, aes(x=logFC,y=P.Value)) +
  geom_point() +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Nominal p-value') +
  ggtitle('Log2FC vs. p-value from DE Analysis of Fluconazole Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(scatter1dbl)
#ggsave(filename='fluconazole_scatterdbl.png', plot=scatter1dbl)

scatter2 <- ggplot(t2_filt, aes(x=logFC,y=P.Value)) +
  geom_point() +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Nominal p-value') +
  ggtitle('Log2FC vs. p-value from DE Analysis of 3-Methylcholanthrene Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(scatter2)
#ggsave(filename='3_methylcholanthrene_scatter.png', plot=scatter2)

scatter2dbl <- ggplot(t2_dblfilt, aes(x=logFC,y=P.Value)) +
  geom_point() +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Nominal p-value') +
  ggtitle('Log2FC vs. p-value from DE Analysis of 3-Methylcholanthrene Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(scatter2dbl)
#ggsave(filename='3_methylcholanthrene_scatterdbl.png', plot=scatter2dbl)

scatter3 <- ggplot(t3_filt, aes(x=logFC,y=P.Value)) +
  geom_point() +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Nominal p-value') +
  ggtitle('Log2FC vs. p-value from DE Analysis of Pirinixic Acid Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(scatter3)
#ggsave(filename='pirinixic_acid_scatter.png', plot=scatter3)

scatter3dbl <- ggplot(t3_dblfilt, aes(x=logFC,y=P.Value)) +
  geom_point() +
  theme_bw() +
  labs(x='Log2 Fold Change',y='Nominal p-value') +
  ggtitle('Log2FC vs. p-value from DE Analysis of Pirinixic Acid Treatment') +
  theme(plot.title = element_text(hjust = 0.5))
#print(scatter3dbl)
#ggsave(filename='pirinixic_acid_scatterdbl.png', plot=scatter3dbl)

# Implementing the log2 fold change cutoff simply reduces the usable number of results for the 3-methylcholanthrene treatment too much.

fluco_deseq <- read.csv("C:/Users/dango/Documents/CAR_PXR.csv", header = TRUE)
fluco_deseq_top10 <- head(dplyr::arrange(fluco_deseq, pvalue), 10)
#write.csv(fluco_deseq_top10, "fluco_deseq.csv")
meth_deseq <- read.csv("C:/Users/dango/Documents/AhR.csv", header = TRUE)
meth_deseq_top10 <- head(dplyr::arrange(meth_deseq, pvalue), 10)
#write.csv(meth_deseq_top10, "meth_deseq.csv")
pirin_deseq <- read.csv("C:/Users/dango/Documents/PPARA.csv", header = TRUE)
pirin_deseq_top10 <- head(dplyr::arrange(pirin_deseq, pvalue), 10)
#write.csv(pirin_deseq_top10, "pirin_deseq.csv")
fluco_limma <- t1_filt
meth_limma <- t2_filt
pirin_limma <- t3_filt

fluco_deseq_refseq <- fluco_deseq_top10$X
meth_deseq_refseq <- meth_deseq_top10$X
pirin_deseq_refseq <- pirin_deseq_top10$X

#fluco_deseq <- dplyr::filter(fluco_deseq, abs(log2FoldChange) > 1.5)
#meth_deseq <- dplyr::filter(meth_deseq, abs(log2FoldChange) > 1.5)
#pirin_deseq <- dplyr::filter(pirin_deseq, abs(log2FoldChange) > 1.5)
#fluco_limma <- dplyr::filter(fluco_limma, abs(logFC) > 1.5)
#meth_limma <- dplyr::filter(meth_limma, abs(logFC) > 1.5)
#pirin_limma <- dplyr::filter(pirin_limma, abs(logFC) > 1.5)

fluco_limma <- tibble::rownames_to_column(fluco_limma, "PROBEID")
meth_limma <- tibble::rownames_to_column(meth_limma, "PROBEID")
pirin_limma <- tibble::rownames_to_column(pirin_limma, "PROBEID")

fluco_deseq_up <- dplyr::filter(fluco_deseq, baseMean > median(baseMean))
meth_deseq_up <- dplyr::filter(meth_deseq, baseMean > median(baseMean))
pirin_deseq_up <- dplyr::filter(pirin_deseq, baseMean > median(baseMean))
fluco_limma_up <- dplyr::filter(fluco_limma, AveExpr > median(AveExpr))
meth_limma_up <- dplyr::filter(meth_limma, AveExpr > median(AveExpr))
pirin_limma_up <- dplyr::filter(pirin_limma, AveExpr > median(AveExpr))

fluco_deseq_down <- dplyr::filter(fluco_deseq, baseMean < median(baseMean))
meth_deseq_down <- dplyr::filter(meth_deseq, baseMean < median(baseMean))
pirin_deseq_down <- dplyr::filter(pirin_deseq, baseMean < median(baseMean))
fluco_limma_down <- dplyr::filter(fluco_limma, AveExpr < median(AveExpr))
meth_limma_down <- dplyr::filter(meth_limma, AveExpr < median(AveExpr))
pirin_limma_down <- dplyr::filter(pirin_limma, AveExpr < median(AveExpr))

refseq_affy_map <- read.csv("C:/Users/dango/Documents/refseq_affy_map.csv")
fluco_deseq_merge <- drop_na(merge(refseq_affy_map, fluco_deseq, by.x="REFSEQ", by.y="X"))
meth_deseq_merge <- drop_na(merge(refseq_affy_map, meth_deseq, by.x="REFSEQ", by.y="X"))
pirin_deseq_merge <- drop_na(merge(refseq_affy_map, pirin_deseq, by.x="REFSEQ", by.y="X"))
fluco_limma_merge <- drop_na(merge(refseq_affy_map, fluco_limma, by="PROBEID"))
meth_limma_merge <- drop_na(merge(refseq_affy_map, meth_limma, by="PROBEID"))
pirin_limma_merge <- drop_na(merge(refseq_affy_map, pirin_limma, by="PROBEID"))

fluco_deseq_down_merge <- drop_na(merge(refseq_affy_map, fluco_deseq_down, by.x="REFSEQ", by.y="X"))
meth_deseq_down_merge <- drop_na(merge(refseq_affy_map, meth_deseq_down, by.x="REFSEQ", by.y="X"))
pirin_deseq_down_merge <- drop_na(merge(refseq_affy_map, pirin_deseq_down, by.x="REFSEQ", by.y="X"))
fluco_deseq_up_merge <- drop_na(merge(refseq_affy_map, fluco_deseq_up, by.x="REFSEQ", by.y="X"))
meth_deseq_up_merge <- drop_na(merge(refseq_affy_map, meth_deseq_up, by.x="REFSEQ", by.y="X"))
pirin_deseq_up_merge <- drop_na(merge(refseq_affy_map, pirin_deseq_up, by.x="REFSEQ", by.y="X"))

fluco_limma_down_merge <- drop_na(merge(refseq_affy_map, fluco_limma_down, by="PROBEID"))
meth_limma_down_merge <- drop_na(merge(refseq_affy_map, meth_limma_down, by="PROBEID"))
pirin_limma_down_merge <- drop_na(merge(refseq_affy_map, pirin_limma_down, by="PROBEID"))
fluco_limma_up_merge <- drop_na(merge(refseq_affy_map, fluco_limma_up, by="PROBEID"))
meth_limma_up_merge <- drop_na(merge(refseq_affy_map, meth_limma_up, by="PROBEID"))
pirin_limma_up_merge <- drop_na(merge(refseq_affy_map, pirin_limma_up, by="PROBEID"))

names(fluco_deseq_merge)[names(fluco_deseq_merge)=="REFSEQ"] <- "REFSEQ.d"
names(fluco_deseq_merge)[names(fluco_deseq_merge)=="PROBEID"] <- "PROBEID.d"
names(meth_deseq_merge)[names(meth_deseq_merge)=="REFSEQ"] <- "REFSEQ.d"
names(meth_deseq_merge)[names(meth_deseq_merge)=="PROBEID"] <- "PROBEID.d"
names(pirin_deseq_merge)[names(pirin_deseq_merge)=="REFSEQ"] <- "REFSEQ.d"
names(pirin_deseq_merge)[names(pirin_deseq_merge)=="PROBEID"] <- "PROBEID.d"

names(fluco_limma_merge)[names(fluco_limma_merge)=="REFSEQ"] <- "REFSEQ.l"
names(fluco_limma_merge)[names(fluco_limma_merge)=="PROBEID"] <- "PROBEID.l"
names(meth_limma_merge)[names(meth_limma_merge)=="REFSEQ"] <- "REFSEQ.l"
names(meth_limma_merge)[names(meth_limma_merge)=="PROBEID"] <- "PROBEID.l"
names(pirin_limma_merge)[names(pirin_limma_merge)=="REFSEQ"] <- "REFSEQ.l"
names(pirin_limma_merge)[names(pirin_limma_merge)=="PROBEID"] <- "PROBEID.l"

names(fluco_deseq_up_merge)[names(fluco_deseq_up_merge)=="REFSEQ"] <- "REFSEQ.d"
names(fluco_deseq_up_merge)[names(fluco_deseq_up_merge)=="PROBEID"] <- "PROBEID.d"
names(meth_deseq_up_merge)[names(meth_deseq_up_merge)=="REFSEQ"] <- "REFSEQ.d"
names(meth_deseq_up_merge)[names(meth_deseq_up_merge)=="PROBEID"] <- "PROBEID.d"
names(pirin_deseq_up_merge)[names(pirin_deseq_up_merge)=="REFSEQ"] <- "REFSEQ.d"
names(pirin_deseq_up_merge)[names(pirin_deseq_up_merge)=="PROBEID"] <- "PROBEID.d"

names(fluco_limma_up_merge)[names(fluco_limma_up_merge)=="REFSEQ"] <- "REFSEQ.l"
names(fluco_limma_up_merge)[names(fluco_limma_up_merge)=="PROBEID"] <- "PROBEID.l"
names(meth_limma_up_merge)[names(meth_limma_up_merge)=="REFSEQ"] <- "REFSEQ.l"
names(meth_limma_up_merge)[names(meth_limma_up_merge)=="PROBEID"] <- "PROBEID.l"
names(pirin_limma_up_merge)[names(pirin_limma_up_merge)=="REFSEQ"] <- "REFSEQ.l"
names(pirin_limma_up_merge)[names(pirin_limma_up_merge)=="PROBEID"] <- "PROBEID.l"

names(fluco_deseq_down_merge)[names(fluco_deseq_down_merge)=="REFSEQ"] <- "REFSEQ.d"
names(fluco_deseq_down_merge)[names(fluco_deseq_down_merge)=="PROBEID"] <- "PROBEID.d"
names(meth_deseq_down_merge)[names(meth_deseq_down_merge)=="REFSEQ"] <- "REFSEQ.d"
names(meth_deseq_down_merge)[names(meth_deseq_down_merge)=="PROBEID"] <- "PROBEID.d"
names(pirin_deseq_down_merge)[names(pirin_deseq_down_merge)=="REFSEQ"] <- "REFSEQ.d"
names(pirin_deseq_down_merge)[names(pirin_deseq_down_merge)=="PROBEID"] <- "PROBEID.d"

names(fluco_limma_down_merge)[names(fluco_limma_down_merge)=="REFSEQ"] <- "REFSEQ.l"
names(fluco_limma_down_merge)[names(fluco_limma_down_merge)=="PROBEID"] <- "PROBEID.l"
names(meth_limma_down_merge)[names(meth_limma_down_merge)=="REFSEQ"] <- "REFSEQ.l"
names(meth_limma_down_merge)[names(meth_limma_down_merge)=="PROBEID"] <- "PROBEID.l"
names(pirin_limma_down_merge)[names(pirin_limma_down_merge)=="REFSEQ"] <- "REFSEQ.l"
names(pirin_limma_down_merge)[names(pirin_limma_down_merge)=="PROBEID"] <- "PROBEID.l"

fluco_merge <- drop_na(merge(fluco_deseq_merge, fluco_limma_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
fluco_plusminus <- dplyr::filter(fluco_merge, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))
meth_merge <- drop_na(merge(meth_deseq_merge, meth_limma_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
meth_plusminus <- dplyr::filter(meth_merge, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))
pirin_merge <- drop_na(merge(pirin_deseq_merge, pirin_limma_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
pirin_plusminus <- dplyr::filter(pirin_merge, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))

fluco_merge_top10 <- dplyr::arrange(fluco_merge, )

fluco_merge_up <- drop_na(merge(fluco_deseq_up_merge, fluco_limma_up_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
fluco_plusminus_up <- dplyr::filter(fluco_merge_up, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))
meth_merge_up <- drop_na(merge(meth_deseq_up_merge, meth_limma_up_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
meth_plusminus_up <- dplyr::filter(meth_merge_up, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))
pirin_merge_up <- drop_na(merge(pirin_deseq_up_merge, pirin_limma_up_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
pirin_plusminus_up <- dplyr::filter(pirin_merge_up, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))

fluco_merge_down <- drop_na(merge(fluco_deseq_down_merge, fluco_limma_up_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
fluco_plusminus_down <- dplyr::filter(fluco_merge_down, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))
meth_merge_down <- drop_na(merge(meth_deseq_down_merge, meth_limma_up_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
meth_plusminus_down <- dplyr::filter(meth_merge_down, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))
pirin_merge_down <- drop_na(merge(pirin_deseq_down_merge, pirin_limma_up_merge, by.x="REFSEQ.d", by.y="REFSEQ.l"))
pirin_plusminus_down <- dplyr::filter(pirin_merge_down, ((logFC > 0 & log2FoldChange > 0) | (logFC < 0 & log2FoldChange < 0)))

N <- nrow(refseq_affy_map)

n0_fluco <- nrow(fluco_plusminus)
n1_fluco <- nrow(fluco_deseq)
n2_fluco <- nrow(fluco_limma)
# Actually had to do some algebra to get this. Wild!
nx_fluco <- abs((n0_fluco*N-n1_fluco*n2_fluco)/(n0_fluco+N-n1_fluco-n2_fluco))
C_fluco <- 2 * nx_fluco / (n1_fluco + n2_fluco)
#corrected_n0_fluco <- 2024.9
#fluco_con <- (2*(corrected_n0_fluco))/(n1_fluco + n2_fluco)

n0_fluco_up <- nrow(fluco_plusminus_up)
n1_fluco_up <- nrow(fluco_deseq_up)
n2_fluco_up <- nrow(fluco_limma_up)
# Actually had to do some algebra to get this. Wild!
nx_fluco_up <- abs((n0_fluco_up*N-n1_fluco_up*n2_fluco_up)/(n0_fluco_up+N-n1_fluco_up-n2_fluco_up))
C_fluco_up <- 2 * nx_fluco_up / (n1_fluco_up + n2_fluco_up)
#corrected_n0_fluco_up <- 2019.8
#fluco_con_up <- (2*(corrected_n0_fluco_up))/(n1_fluco_up + n2_fluco_up)

n0_fluco_down <- nrow(fluco_plusminus_down)
n1_fluco_down <- nrow(fluco_deseq_down)
n2_fluco_down <- nrow(fluco_limma_down)
# Actually had to do some algebra to get this. Wild!
nx_fluco_down <- abs((n0_fluco_down*N-n1_fluco_down*n2_fluco_down)/(n0_fluco_down+N-n1_fluco_down-n2_fluco_down))
C_fluco_down <- 2 * nx_fluco_down / (n1_fluco_down + n2_fluco_down)
#corrected_n0_fluco_down <- 231.8
#fluco_con_down <- (2*(corrected_n0_fluco_down))/(n1_fluco_down + n2_fluco_down)

n0_meth <- nrow(meth_plusminus)
n1_meth <- nrow(meth_deseq)
n2_meth <- nrow(meth_limma)
# Actually had to do some algebra to get this. Wild!
nx_meth <- abs((n0_meth*N-n1_meth*n2_meth)/(n0_meth+N-n1_meth-n2_meth))
C_meth <- 2 * nx_meth / (n1_meth + n2_meth)
#corrected_n0_fluco <- 2024.9
#fluco_con <- (2*(corrected_n0_fluco))/(n1_fluco + n2_fluco)

n0_meth_up <- nrow(meth_plusminus_up)
n1_meth_up <- nrow(meth_deseq_up)
n2_meth_up <- nrow(meth_limma_up)
# Actually had to do some algebra to get this. Wild!
nx_meth_up <- abs((n0_meth_up*N-n1_meth_up*n2_meth_up)/(n0_meth_up+N-n1_meth_up-n2_meth_up))
C_meth_up <- 2 * nx_meth_up / (n1_meth_up + n2_meth_up)
#corrected_n0_meth_up <- 27.4
#meth_con_up <- (2*(corrected_n0_meth_up))/(n1_meth_up + n2_meth_up)

n0_meth_down <- nrow(meth_plusminus_down)
n1_meth_down <- nrow(meth_deseq_down)
n2_meth_down <- nrow(meth_limma_down)
# Actually had to do some algebra to get this. Wild!
nx_meth_down <- abs((n0_meth_down*N-n1_meth_down*n2_meth_down)/(n0_meth_down+N-n1_meth_down-n2_meth_down))
C_meth_down <- 2 * nx_meth_down / (n1_meth_down + n2_meth_down)
#corrected_n0_meth_down <- 27.4
#meth_con_down <- (2*(corrected_n0_meth_down))/(n1_meth_down + n2_meth_down)

n0_pirin <- nrow(pirin_plusminus)
n1_pirin <- nrow(pirin_deseq)
n2_pirin <- nrow(pirin_limma)
# Actually had to do some algebra to get this. Wild!
nx_pirin <- abs((n0_pirin*N-n1_pirin*n2_pirin)/(n0_pirin+N-n1_pirin-n2_pirin))
C_pirin <- 2 * nx_pirin / (n1_pirin + n2_pirin)

n0_pirin_up <- nrow(pirin_plusminus_up)
n1_pirin_up <- nrow(pirin_deseq_up)
n2_pirin_up <- nrow(pirin_limma_up)
# Actually had to do some algebra to get this. Wild!
nx_pirin_up <- abs((n0_pirin_up*N-n1_pirin_up*n2_pirin_up)/(n0_pirin_up+N-n1_pirin_up-n2_pirin_up))
C_pirin_up <- 2 * nx_pirin_up / (n1_pirin_up + n2_pirin_up)
#corrected_n0_pirin_up <- 4072.7
#pirin_con_up <- (2*(corrected_n0_pirin_up))/(n1_pirin_up + n2_pirin_up)

n0_pirin_down <- nrow(pirin_plusminus_down)
n1_pirin_down <- nrow(pirin_deseq_down)
n2_pirin_down <- nrow(pirin_limma_down)
# Actually had to do some algebra to get this. Wild!
nx_pirin_down <- abs((n0_pirin_down*N-n1_pirin_down*n2_pirin_down)/(n0_pirin_down+N-n1_pirin_down-n2_pirin_down))
C_pirin_down <- 2 * nx_pirin_down / (n1_pirin_down + n2_pirin_down)
#corrected_n0_pirin_down <- 4072.7
#pirin_con_down <- (2*(corrected_n0_pirin_down))/(n1_pirin_down + n2_pirin_down)

treatments <- c("Fluconazole","3-Methylcholanthrene", "Pirinixic Acid")
treatments_triple <- c("Fluconazole","Fluconazole","Fluconazole","3-Methylcholanthrene","3-Methylcholanthrene","3-Methylcholanthrene",
                       "Pirinixic Acid","Pirinixic Acid","Pirinixic Acid")
concordances <- c(C_fluco, C_meth, C_pirin)
status <- c("ABOVE MEDIAN", "BELOW MEDIAN", "OVERALL", "ABOVE MEDIAN", "BELOW MEDIAN", "OVERALL", "ABOVE MEDIAN", "BELOW MEDIAN", "OVERALL")
concordances_triple <- c(C_fluco_up, C_fluco_down, C_fluco, C_meth_up, C_meth_down, C_meth, C_pirin_up, C_pirin_down, C_pirin)
scatter_df <- data.frame(treatments, concordances)
barplot_df <- data.frame(treatments_triple, status, concordances_triple)


DEG_counts_deseq <- c(n1_fluco, n1_meth, n1_pirin)
scatter_deseq <- ggplot(scatter_df, aes(x=DEG_counts_deseq, y=concordances)) +
  geom_point() +
  theme_bw() +
  geom_text_repel(aes(label = treatments)) + 
  xlab("Number of Differentially Expressed Genes") +
  ylab("Concordance") +
  ggtitle("Concordance of Chemical Treatments in DESeq2 Analysis") +
  theme(plot.title = element_text(hjust = 0.5))

DEG_counts_limma <- c(n2_fluco, n2_meth, n2_pirin)
scatter_limma <- ggplot(scatter_df, aes(x=DEG_counts_limma, y=concordances)) +
  geom_point() +
  theme_bw() +
  geom_text_repel(aes(label = treatments)) + 
  xlab("Number of Differentially Expressed Genes") +
  ylab("Concordance") +
  ggtitle("Concordance of Chemical Treatments in Limma Analysis") +
  theme(plot.title = element_text(hjust = 0.5))

conc_barplot <- ggplot(barplot_df, aes(x=treatments_triple, y=concordances_triple, fill=status)) +
  geom_bar(stat="identity", position="dodge") +
  theme_bw() +
  xlab("Treatments") +
  ylab("Concordance") + 
  ggtitle("Effect of Expression Level on Concordance") +
  theme(plot.title = element_text(hjust = 0.5))

#print(conc_barplot)
