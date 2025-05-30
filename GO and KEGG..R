# =======================================
# GO and KEGG pathway enrichment analysis
# =======================================

# Install and load library

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Convert the gene list from excel file to Entrez ID
gene_list <- Core_targets_of_compounds$`Gene symbol`

gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# Biological Process (BP) analysis============================================================
go_bp_results <- enrichGO(gene = gene_ids$ENTREZID, 
                          OrgDb = org.Hs.eg.db, 
                          keyType = "ENTREZID", 
                          ont = "BP",  # Biological Process
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05,
                          readable = TRUE)


# Cellular Component (CC) analysis
go_cc_results <- enrichGO(gene = gene_ids$ENTREZID, 
                          OrgDb = org.Hs.eg.db, 
                          keyType = "ENTREZID", 
                          ont = "CC",  # Cellular Component
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05,
                          readable = TRUE)


# Molecular Function (MF) analysis
go_mf_results <- enrichGO(gene = gene_ids$ENTREZID, 
                          OrgDb = org.Hs.eg.db, 
                          keyType = "ENTREZID", 
                          ont = "MF",  # Molecular Function
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05,
                          readable = TRUE)
# KEGG pathway analysis
kegg_results <- enrichKEGG(gene = gene_ids$ENTREZID,
                          organism     = 'hsa',
                          pvalueCutoff = 0.05)
kegg_results <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Visualization of Dot_plot
pdf("Go_KEGG_Dotplots.pdf", width = 8, height = 6)

print(dotplot(go_bp_results, showCategory = 10))
print(dotplot(go_cc_results, showCategory = 10))
print(dotplot(go_mf_results, showCategory = 10))
print(barplot(kegg_results, showCategory = 20))

dev.off()

# Save the results
dir.create("GO_Plots", showWarnings = FALSE)

go_bp_results <- pairwise_termsim(go_bp_results)
go_cc_results <- pairwise_termsim(go_cc_results)
go_mf_results <- pairwise_termsim(go_mf_results)
kegg_results <- pairwise_termsim(kegg_results)

write.csv(go_bp_results, "GO_enrichment_results_BP.csv")
write.csv(go_cc_results, "GO_enrichment_results_CC.csv")
write.csv(go_mf_results, "GO_enrichment_results_MF.csv")
write.csv(kegg_results, "KEGG_enrichment_results.csv")


# Emapplot and Cnetplot========================================================
## ---- BP ----
pdf("GO_Plots/emapplot_BP.pdf", width = 10, height = 8)
emapplot(go_bp_results, showCategory = 20, layout = "kk")
dev.off()

pdf("GO_Plots/cnetplot_BP.pdf", width = 12, height = 10)
cnetplot(go_bp_results, showCategory = 10, circular = FALSE, colorEdge = TRUE)
dev.off()

## ---- CC ----
pdf("GO_Plots/emapplot_CC.pdf", width = 10, height = 8)
emapplot(go_cc_results, showCategory = 20, layout = "kk")
dev.off()

pdf("GO_Plots/cnetplot_CC.pdf", width = 12, height = 10)
cnetplot(go_cc_results, showCategory = 10, circular = FALSE, colorEdge = TRUE)
dev.off()

## ---- MF  ----
pdf("GO_Plots/emapplot_MF.pdf", width = 10, height = 8)
emapplot(go_mf_results, showCategory = 20, layout = "kk")
dev.off()

pdf("GO_Plots/cnetplot_MF.pdf", width = 12, height = 10)
cnetplot(go_mf_results, showCategory = 10, circular = FALSE, colorEdge = TRUE)
dev.off()

## ---- KEGG  ----
pdf("GO_Plots/emapplot_KEGG.pdf", width = 10, height = 8)
emapplot(kegg_results, showCategory = 20, layout = "kk")
dev.off()

pdf("GO_Plots/cnetplot_KEGG.pdf", width = 12, height = 10)
cnetplot(kegg_results, showCategory = 10, circular = FALSE, colorEdge = TRUE)
dev.off()
