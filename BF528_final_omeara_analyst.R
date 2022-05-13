p4p7 <- read.table("C:/Users/dango/Documents/gene_exp.diff", header = TRUE)
p4p7_sort <- dplyr::arrange(p4p7, q_value)
p4p7_hist <- ggplot(p4p7_sort, aes(x=log2.fold_change.)) +
  geom_histogram(binwidth = 0.5) + 
  theme_bw() +
  xlab("Log2 Fold Change") +
  ylab("Count") + 
  ggtitle("Histogram of Log2 Fold Change Values for All Rows") +
  theme(plot.title = element_text(hjust = 0.5))
print(p4p7_hist)
p4p7_sig <- dplyr::filter(p4p7_sort, significant == "yes")
p4p7_hist_sig <- ggplot(p4p7_sig, aes(x=log2.fold_change.)) +
  geom_histogram(binwidth = 0.5) + 
  theme_bw() +
  xlab("Log2 Fold Change") +
  ylab("Count") + 
  ggtitle("Histogram of Log2 Fold Change Values for Significant Rows") +
  theme(plot.title = element_text(hjust = 0.5))
print(p4p7_hist_sig)
p4p7_sig_plus <- dplyr::filter(p4p7_sig, log2.fold_change. > 0)
p4p7_sig_minus <- dplyr::filter(p4p7_sig, log2.fold_change. < 0)
plus_genes <- p4p7_sig_plus$gene
minus_genes <- p4p7_sig_minus$gene
write(plus_genes, file = "plus_genes")
write(minus_genes, file = "minus_genes")
p4p7_fcsort <- dplyr::arrange(p4p7_sig, abs(log2.fold_change.))
p4p7_fcsort <- head(p4p7_fcsort, 457)
  # need to exclude two rows with missing values
p4p7_top10abs <- tail(p4p7_fcsort, 10)
  # names, FPKM values, log fold change, p-value, and q-value
top10_final <- dplyr::select(p4p7_top10abs, gene, value_1, value_2, log2.fold_change.,p_value, q_value)
write.csv(top10_final, "project2_top10.csv")
