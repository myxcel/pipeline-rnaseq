library(clusterProfiler)
library(org.Mm.eg.db)   # base de données d'annotations pour la souris
library(ggplot2)
library(dplyr)

# ── 1. Charger les résultats DESeq2 ────────────────
res_df <- read.csv("results/deseq2/DE_results.csv", row.names = 1)

# Listes de gènes up et down-régulés
genes_up   <- rownames(filter(res_df, padj < 0.05, log2FoldChange >  1))
genes_down <- rownames(filter(res_df, padj < 0.05, log2FoldChange < -1))
all_genes  <- rownames(res_df)   # background = tous les gènes testés

# ── 2. Convertir les symboles en Entrez IDs ────────
# clusterProfiler a besoin des Entrez Gene IDs (format NCBI)
convert_ids <- function(gene_symbols) {
  bitr(gene_symbols,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Mm.eg.db)$ENTREZID
}

entrez_up      <- convert_ids(genes_up)
entrez_down    <- convert_ids(genes_down)
entrez_universe <- convert_ids(all_genes)

# ── 3. Enrichissement GO — Processus Biologiques ───
ego_up <- enrichGO(
  gene          = entrez_up,
  universe      = entrez_universe,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",             # Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE              # afficher les symboles plutôt que Entrez IDs
)

ego_down <- enrichGO(
  gene          = entrez_down,
  universe      = entrez_universe,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1,
  readable      = TRUE
)

# ── 4. Figure : Dot plot UP ────────────────────────
dot_up <- dotplot(ego_up,
  showCategory = 20,
  title        = "GO Biological Process — Gènes up-régulés (adulte vs embryonnaire)",
  font.size    = 9
) + theme_bw()

ggsave("results/deseq2/GO_dotplot_UP.png", dot_up, width = 10, height = 8, dpi = 300)

dot_down <- dotplot(ego_down,
  showCategory = 20,
  title        = "GO Biological Process — Gènes down-régulés",
  font.size    = 9
) + theme_bw()

ggsave("results/deseq2/GO_dotplot_DOWN.png", dot_down, width = 10, height = 8, dpi = 300)

# ── 5. Exporter les tables ─────────────────────────
write.csv(as.data.frame(ego_up),   "results/deseq2/GO_up.csv",   quote = FALSE)
write.csv(as.data.frame(ego_down), "results/deseq2/GO_down.csv", quote = FALSE)

cat("✅ Enrichissement GO terminé.\n")
