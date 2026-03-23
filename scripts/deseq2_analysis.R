# ─────────────────────────────────────────────────
# Script DESeq2 — Analyse différentielle RNAseq
# Contexte : cortex murin en développement
# ─────────────────────────────────────────────────

library(DESeq2)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# ── 1. Charger la matrice de comptage ──────────────
counts_raw <- read.table(
  "results/counts/counts_matrix.txt",
  header = TRUE, skip = 1, row.names = 1
)

# Garder uniquement les colonnes de comptage (enlever Chr, Start, End, Strand, Length)
counts <- counts_raw[, 6:ncol(counts_raw)]

# Nettoyer les noms de colonnes (enlever le chemin du BAM)
colnames(counts) <- gsub("results.aligned.|_Aligned.*", "", colnames(counts))

# ── 2. Table des métadonnées (coldata) ─────────────
# ADAPTE les noms et conditions à ton dataset réel !
coldata <- data.frame(
  sample    = colnames(counts),
  condition = factor(c("embryonic", "embryonic", "embryonic",
                        "adult",     "adult",     "adult")),
  row.names = colnames(counts)
)

# Vérification : les colonnes de counts correspondent aux lignes de coldata
stopifnot(all(rownames(coldata) == colnames(counts)))

# ── 3. Créer l'objet DESeq2 ────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),   # DESeq2 requiert des entiers
  colData   = coldata,
  design    = ~ condition       # comparer par condition
)

# Filtrer les gènes à très faible expression (< 10 reads au total)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# ── 4. Normalisation et test statistique ───────────
dds <- DESeq(dds)
# Cette ligne fait TOUT : normalisation par size factors,
# estimation de la dispersion, test du ratio de vraisemblance (Wald test)

# ── 5. Extraire les résultats ──────────────────────
res <- results(dds,
  contrast     = c("condition", "adult", "embryonic"),
  alpha        = 0.05,
  lfcThreshold = 0
)

# Rétrécissement du log2FC (plus stable pour les gènes à faible comptage)
res_shrunk <- lfcShrink(dds,
  coef = "condition_adult_vs_embryonic",
  type = "apeglm"
)

# Résumé
summary(res)

# Exporter les résultats triés par p-value ajustée
res_df <- as.data.frame(res_shrunk) %>%
  arrange(padj) %>%
  filter(!is.na(padj))

write.csv(res_df, "results/deseq2/DE_results.csv", quote = FALSE)

# ── 6. Figure 1 : PCA ──────────────────────────────
# Transformation VST (variance-stabilizing) pour la visualisation
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pca_var  <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(vjust = -0.8, size = 3) +
  xlab(paste0("PC1: ", pca_var[1], "% variance")) +
  ylab(paste0("PC2: ", pca_var[2], "% variance")) +
  ggtitle("PCA — Cortex murin : embryonnaire vs adulte") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("embryonic" = "#2E6DA4", "adult" = "#E84A1F"))

ggsave("results/deseq2/PCA_plot.png", pca_plot, width = 7, height = 6, dpi = 300)

# ── 7. Figure 2 : Volcano plot ─────────────────────
# Gènes d'intérêt neurodéveloppement à annoter (alignés sur les travaux Godin)
genes_of_interest <- c("Satb2", "Sox2", "Eomes", "Tbr1", "Ctip2", "Pax6")

volcano <- EnhancedVolcano(res_df,
  lab              = rownames(res_df),
  x                = "log2FoldChange",
  y                = "padj",
  selectLab        = genes_of_interest,
  pCutoff          = 0.05,
  FCcutoff         = 1,
  pointSize        = 2,
  labSize          = 4,
  title            = "Cortex murin : adulte vs embryonnaire",
  subtitle         = "DESeq2 — log2FC > 1, FDR < 5%",
  caption          = paste("Gènes significatifs :",
                            sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                                na.rm = TRUE)),
  col              = c("grey80", "grey50", "#2E6DA4", "#E84A1F"),
  colAlpha         = 0.8,
  legendPosition   = "right"
)

ggsave("results/deseq2/Volcano_plot.png", volcano, width = 10, height = 8, dpi = 300)

# ── 8. Figure 3 : Heatmap top 50 gènes ────────────
top50 <- res_df %>%
  filter(padj < 0.05) %>%
  slice_min(padj, n = 50) %>%
  rownames()

mat <- assay(vsd)[top50, ]
mat <- mat - rowMeans(mat)   # centrer par gène

annotation_col <- data.frame(
  Condition = coldata$condition,
  row.names = rownames(coldata)
)

pheatmap(mat,
  annotation_col   = annotation_col,
  show_rownames    = TRUE,
  show_colnames    = TRUE,
  cluster_cols     = TRUE,
  cluster_rows     = TRUE,
  color            = colorRampPalette(c("#2E6DA4", "white", "#E84A1F"))(100),
  fontsize_row     = 7,
  main             = "Top 50 gènes différentiels — Cortex murin",
  filename         = "results/deseq2/Heatmap_top50.png",
  width            = 10,
  height           = 12
)

cat("✅ Analyse DESeq2 terminée. Fichiers dans results/deseq2/\n")
