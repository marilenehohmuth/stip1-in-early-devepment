#-- Load
library(readr); library(Matrix); library(Seurat)

#-- Data
matrix_file <- "/dados/marcelo/c_felix/matrix.mtx"
features_file <- "/dados/marcelo/c_felix/features.tsv"
counts <- readMM(matrix_file)
features <- read_delim(features_file, delim = "\t", col_names = FALSE, show_col_types = FALSE)
rownames(counts) <- features$V1
meta.data <- read_delim("/dados/marcelo/c_felix/meta.tab", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
rownames(meta.data) <- meta.data$cell
colnames(counts) <- meta.data$cell
seurat <- CreateSeuratObject(counts = counts, 
                             min.cells = 0, 
                             min.features = 0, 
                             meta.data = meta.data)
seurat@meta.data$celltype <- meta.data$celltype
seurat@meta.data$umapX <- meta.data$umapX
seurat@meta.data$umapY <- meta.data$umapY
ggplot(seurat@meta.data, aes(x = umapX, y = umapY, color = celltype)) +
  geom_point()
seurat <- NormalizeData(seurat)
Idents(seurat) <- "celltype"
markers <- FindAllMarkers(seurat, only.pos = TRUE)
markers <- markers[which(markers$p_val_adj < 0.05), ]
# save(seurat, markers, file = "/dados/marcelo/c_felix/my_image.RData")

#-- Analysis
load(file = "/dados/marcelo/c_felix/my_image.RData")
markers <- markers[which(markers$pct.1 > 0.25), ]
markers <- markers[which(markers$avg_log2FC > 1), ]
all_genes <- sort(rownames(seurat))
files <- c("/dados/marcelo/c_felix/HTKO_vs_WT_DESeq2_filtered_DEGs_pvalue0.05.csv", 
           "/dados/marcelo/c_felix/TGA_vs_WT_DESeq2_filtered_DEGs_pvalue0.05.csv", 
           "/dados/marcelo/c_felix/TPR1_vs_WT_DESeq2_filtered_DEGs_pvalue0.05.csv")
my_list_01 <- list()
for (file in 1:length(files)) {
  results <- read_csv(file = files[file])
  results <- results[which(results$log2FoldChange > 1), ]
  results <- results[order(results$log2FoldChange, decreasing = TRUE), ]
  results <- results[order(results$padj, decreasing = FALSE), ]
  deg_genes <- results$...1
  ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  results <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                   filters = "mgi_symbol",
                   values = deg_genes,
                   mart = ensembl)
  deg_genes <- results$ensembl_gene_id
  my_list_02 <- list()
  for (i in 1:length(names(table(markers$cluster)))) {
    genes_in_term <- sort(markers[which(markers$cluster == names(table(markers$cluster))[i]), ]$gene)
    A <- length(intersect(genes_in_term, deg_genes))
    B <- length(setdiff(deg_genes, genes_in_term))
    C <- length(setdiff(genes_in_term, deg_genes))
    D <- length(setdiff(all_genes, union(genes_in_term, deg_genes)))
    fisher_test <- fisher.test(matrix(c(A, B, C, D), nrow = 2,
                                      dimnames = list(c("DEG", "Not_DEG"),
                                                      c("In_Term", "Not_In_Term"))))
    my_list_02[[i]] <- fisher_test$p.value
  }
  tmp <- data.frame(x = names(table(markers$cluster)), 
                    y = unlist(my_list_02))
  tmp$p.adjust <- p.adjust(tmp$y)
  tmp <- tmp[order(tmp$x, decreasing = FALSE), ]
  tmp <- tmp[order(tmp$p.adjust, decreasing = FALSE), ]
  my_list_01[[file]] <- tmp[which(tmp$p.adjust < 0.05), ]
}
tmp_02 <- my_list_01[[2]]
tmp_03 <- my_list_01[[3]]
all_terms <- unique(c(tmp_02$x, tmp_03$x))
p_adjust_matrix <- data.frame(
  Term = all_terms,
  p.adjust_tmp_02 = sapply(all_terms, function(term) {
    if (term %in% tmp_02$x) {
      tmp_02$p.adjust[tmp_02$x == term]
    } else {
      NA
    }
  }),
  p.adjust_tmp_03 = sapply(all_terms, function(term) {
    if (term %in% tmp_03$x) {
      tmp_03$p.adjust[tmp_03$x == term]
    } else {
      NA
    }
  })
)
p_adjust_matrix <- p_adjust_matrix[order(p_adjust_matrix$Term, decreasing = FALSE), ]
tmp <- p_adjust_matrix[, 2:3]
tmp <- -log10(tmp)
colnames(tmp) <- c("TGA_vs_WT", "TPR1_vs_WT")
tmp$Cluster <- rownames(tmp)
tmp <- tmp %>%
  pivot_longer(cols = starts_with("T"), 
               names_to = "Condition", 
               values_to = "Value")
png("/dados/marcelo/Camila data and mouse gastrulation and early organogenesis.png", 
    width = 6, height = 6, units = "in", res = 300)
ggplot(tmp, aes(x = Condition, y = Cluster, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "gray", limits = c(0, max(tmp$Value, na.rm = TRUE))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "", x = "Contrast", y = "Clusters", fill = expression(-log[10](adj. ~ italic(p)-value)))
dev.off()
