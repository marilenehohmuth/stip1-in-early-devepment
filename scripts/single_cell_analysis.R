# ANALYSIS OF A PUBLICLY AVAILABLE SINGLE-CELL RNA-SEQ DATASET ##############
# @ This script generates the plots present in the first main figure of the # 
# @ manuscript + supplemental material.                                     #
#############################################################################

##########################
#### Loading packages ####
##########################

library(GEOquery) # v2.68.0
library(fusca) # v1.3.1
library(ggplot2) # v3.4.4
library(dplyr) # v1.1.4
library(RColorBrewer) # v1.1_3
library(stringr) # v1.5.1
library(ggrepel) # v0.9.4
library(ggpubr) # v0.6.0
library(clusterProfiler) # v4.8.1
library(org.Mm.eg.db) # v3.17.0
library(corrplot) # v0.92
library(cowplot) # v1.1.1

###################################
#### Downloading data from GEO ####
###################################

# Download raw data.
getGEOSuppFiles("GSE45719", baseDir = paste0(getwd(), "/data/"))

# Extract files.
dir.create(paste0(getwd(), "/data/GSE45719/files/"))
untar(
  paste0(getwd(), "/data/GSE45719/GSE45719_RAW.tar"), 
  exdir = paste0(getwd(), "/data/GSE45719/files/")
)

# Iterate over each downloaded file to create a count matrix and metadata table with all cells...
count_matrix <- data.frame(matrix(nrow = 22958, ncol = 0))
metadata <- data.frame(matrix(nrow = 0, ncol = 0))
for(file in list.files(paste0(getwd(), "/data/GSE45719/files/"))) {
  
  # Get elements from file name.
  elements <- unlist(str_split(file, pattern = "_"))
  stage <- elements[2]
  cell_id <- paste(elements[!grepl("expression", elements)], collapse = "_")

  # Get developmental stage of a given cell from file name.
  stage_plus_cell_id <- paste0(stage, "_", cell_id)

  # Read file of a given cell.
  cell <- read.csv(
    gzfile(paste0(getwd(), "/data/GSE45719/files/", file)),
    sep = "\t",
    header = TRUE,
    row.names = NULL
  )

  # Extract reads of a given cell.
  cell_reads <- as.data.frame(cell$reads)
  colnames(cell_reads) <- cell_id
  rownames(cell_reads) <- make.unique(cell$X.Gene_symbol, sep = ".")

  # Create dataframe with cell metadata.
  cell_info <- data.frame(stage)
  rownames(cell_info) <- cell_id

  #identical(rownames(cell_reads), rownames(count_matrix))

  # Add cell to count matrix.
  count_matrix <- cbind(count_matrix, cell_reads)

  # Add cell to metadata table.
  metadata <- rbind(metadata, cell_info)
}

#########################
#### Processing data ####
#########################

metadata <- metadata %>% mutate(stage_refined = case_when(
  stage == "16cell" ~ "16-cell",
  stage == "4cell" ~ "4-cell",
  stage == "8cell" ~ "8-cell",
  stage == "BXC" ~ "Liver",
  stage == "C57twocell" ~ "2-cell (C57)",
  stage == "early2cell" ~ "Early 2-cell",
  stage == "earlyblast" ~ "Early blastocyst",
  stage == "late2cell" ~ "Late 2-cell",
  stage == "lateblast" ~ "Late blastocyst",
  stage == "mid2cell" ~ "Intermediate 2-cell",
  stage == "midblast" ~ "Intermediate blastocyst",
  stage == "zy1" | stage == "zy2" | stage == "zy3" | stage == "zy4" ~ "Zygote",
  stage == "fibroblast" ~ "Fibroblast"
))

# Filter out fibroblasts and liver cells.
metadata <- metadata[!(metadata$stage_refined %in% c("Fibroblast", "Liver")),]
count_matrix <- count_matrix %>% dplyr::select(rownames(metadata))

# Create CellRouter object.
cellrouter <- CreateCellRouter(
  as.matrix(count_matrix), 
  min.cells = 0, 
  min.genes = 0, 
  is.expr = 0
)

cellrouter@assays$RNA@rawdata <- cellrouter@assays$RNA@rawdata[rownames(cellrouter@assays$RNA@ndata), colnames(cellrouter@assays$RNA@ndata)]

# Add stage information to the CellRouter object.
cellrouter@assays$RNA@sampTab$stage <- metadata$stage[match(rownames(cellrouter@assays$RNA@sampTab), rownames(metadata))]
cellrouter@assays$RNA@sampTab$stage_refined <- metadata$stage_refined[match(rownames(cellrouter@assays$RNA@sampTab), rownames(metadata))]
cellrouter@assays$RNA@sampTab <- cellrouter@assays$RNA@sampTab %>% mutate(stage_refined_color = case_when(
  stage_refined == "16-cell" ~ "#A6CEE3",
  stage_refined == "4-cell" ~ "#1F78B4",
  stage_refined == "8-cell" ~ "#B2DF8A",
  stage_refined == "2-cell (C57)" ~ "#33A02C",
  stage_refined == "Early 2-cell" ~ "#FB9A99",
  stage_refined == "Early blastocyst" ~ "#E31A1C",
  stage_refined == "Late 2-cell" ~ "#FDBF6F",
  stage_refined == "Late blastocyst" ~ "#FF7F00",
  stage_refined == "Intermediate 2-cell" ~ "#CAB2D6", 
  stage_refined == "Intermediate blastocyst" ~ "#6A3D9A",
  stage_refined == "Zygote" ~ "#FFFF99"
))

# Normalize data.
cellrouter <- Normalize(cellrouter)

# Scale data according to all genes.
cellrouter <- scaleData(cellrouter, blocksize = nrow(cellrouter@assays$RNA@ndata))

# Perform Principal Component Analysis (PCA).
cellrouter <- computePCA(cellrouter, num.pcs = 50, seed = 42)

# Plot standard deviation of the principal components.
pdf(
  paste0(getwd(), '/results/single_cell/PCA_sdev_PCs.pdf'),
  width = 6,
  height = 4
)
plot(cellrouter@pca$sdev, xlab='PC', ylab='Standard deviation of PC')
dev.off()

# Determine order of stages.
cellrouter@assays$RNA@sampTab$stage_refined <- factor(
  cellrouter@assays$RNA@sampTab$stage_refined, 
  levels = c('Zygote', '2-cell (C57)', 'Early 2-cell', 'Intermediate 2-cell', 'Late 2-cell', '4-cell', '8-cell', '16-cell', 'Early blastocyst', 'Intermediate blastocyst', 'Late blastocyst') 
)

# Extract PCA coordinates & add stage information.
dim_red_coords <- cellrouter@rdimension %>% dplyr::select(PC1, PC2)
dim_red_coords$stage_refined <- cellrouter@assays$RNA@sampTab$stage_refined[match(rownames(dim_red_coords), rownames(cellrouter@assays$RNA@sampTab))]

# Get median of PC1 and PC2 values to have positions for labels in PCA plot.
dim_red_coords$median_x <- NA
dim_red_coords$median_y <- NA
for(stage in dim_red_coords$stage_refined %>% unique) {
  median_x <- dim_red_coords$PC1[dim_red_coords$stage_refined == stage] %>% median
  median_y <- dim_red_coords$PC2[dim_red_coords$stage_refined == stage] %>% median
  dim_red_coords$median_x[dim_red_coords$stage_refined == stage] <- as.numeric(median_x)
  dim_red_coords$median_y[dim_red_coords$stage_refined == stage] <- as.numeric(median_y)
}
dim_red_coords_unique <- dim_red_coords[!duplicated(dim_red_coords$median_x),]

# Plot PCA.
pdf(
  paste0(getwd(), '/results/single_cell/PCA.pdf'),
  width = 5.5,
  height = 4
)
ggplot(
  dim_red_coords,
  aes(x = PC1, y = PC2, color = stage_refined)
) +
  geom_point(alpha = 0.5) +
  theme_void() +
  scale_color_manual(
    values = c(
      "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
      "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99"
    ),
    name = "Pre-implantation stage"
  ) +
  geom_text_repel(
    data = dim_red_coords_unique,
    aes(x = median_x, y = median_y, label = stage_refined),
    color = "black",
    show.legend = FALSE
  )
dev.off()

# Build KNN.
cellrouter <- buildKNN(
  cellrouter, 
  k = 15, 
  column.ann = 'stage_refined', 
  num.pcs = 20, 
  sim.type = 'jaccard'
)

# Plot KNN.
plotKNN(
  cellrouter, 
  reduction.type = 'pca', 
  column.ann = 'stage_refined', 
  column.color = 'stage_refined_color',
  width = 5, 
  height = 3.5, 
  filename = paste0(getwd(), '/results/single_cell/KNN.pdf')
)

write.table(
  cellrouter@graph$edges, 
  file = paste0(getwd(), "/results/single_cell/cell_edge_weighted_network.txt"), 
  sep = '\t', 
  row.names = FALSE, 
  col.names = FALSE,  
  quote = FALSE
)

# Add gene expression to dataframe containing PCA coordinates.
genes <- c("Stip1", "Sox2", "Pou5f1", "Nanog", "Stat3", "Klf4", "Smad1", "Zfp42", "Nr0b1", "Esrrb")
transposed_scaled <- as.matrix(t(cellrouter@assays$RNA@scale.data)) %>% as.data.frame

dim_red_coords_genes <- data.frame()
for(cell in transposed_norm_counts %>% rownames) {
  for(gene in genes) {
    PC1 <- dim_red_coords$PC1[rownames(dim_red_coords) == cell]
    PC2 <- dim_red_coords$PC2[rownames(dim_red_coords) == cell]
    gene_exp <- transposed_scaled[,gene][rownames(transposed_scaled) == cell]
    dim_red_coords_genes <- rbind(dim_red_coords_genes, c(cell, PC1, PC2, gene, gene_exp))
  }
}
colnames(dim_red_coords_genes) <- c("cell", "PC1", "PC2", "gene", "gene_exp")

pdf(
  paste0(getwd(), '/results/single_cell/PCA_Stip1_and_Pluripotency_genes.pdf'),
  width = 10,
  height = 5
)
ggplot(
  dim_red_coords_genes %>% dplyr::arrange(gene_exp),
  aes(x = as.numeric(PC1), y = as.numeric(PC2), color = as.numeric(gene_exp))
) +
  geom_point() +
  theme_classic() +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "italic")
  ) +
  scale_color_gradient2(
    midpoint = 0,
    mid = "white",
    low = "blue",
    high = "orange",
    name = "Scaled gene\nexpression\n"
  ) +
  facet_wrap(~gene, nrow = 2) +
  xlab("PC1") +
  ylab("PC2")
dev.off()

# Plot Stip1 expression across developmental stages.
transposed_norm_counts <- t(cellrouter@assays$RNA@ndata) %>% as.data.frame
cellrouter@assays$RNA@sampTab$Stip1 <- transposed_norm_counts$Stip1[match(rownames(cellrouter@assays$RNA@sampTab), rownames(transposed_norm_counts))]
pdf(
  paste0(getwd(), '/results/single_cell/violin_Stip1_across_stages.pdf'),
  width = 6, 
  height = 4
)
ggplot(
  cellrouter@assays$RNA@sampTab,
  aes(x = stage_refined, y = Stip1)
) +
  geom_violin(aes(fill = stage_refined), scale = "width") +
  geom_boxplot(outlier.shape = NA, width = 0.25, fill = "white") +
  scale_fill_manual(
    values = cellrouter@assays$RNA@sampTab$stage_refined_color %>% unique,
    name = "Stage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
  ) +
  xlab("Pre-implantation stage") +
  ylab(expression("Normalized"~italic("Stip1")~"expression")) +
  stat_compare_means( # Compare each stage against all other stages!
    label = "p.signif", 
    method = "wilcox.test", 
    ref.group = ".all.",
    size = 2,
  )
dev.off()

#############################################
#### Performing cell trajectory analysis ####
#############################################

find_trajectories <- function(
  cellrouter_,
  source,
  targets
) {

  cellrouter_ <- findPathsRJava(
    cellrouter_, 
    sources = source, 
    targets = targets, 
    column = 'stage_refined', 
    maindir = paste0(getwd(), '/results/single_cell/'), 
    method = 'graph'
  )
  
  # Process trajectories.
  cellrouter_ <- processTrajectories(
    cellrouter_, 
    genes = rownames(cellrouter_@assays$RNA@ndata), 
    path.rank = 'rank', 
    num.cells = 0, 
    neighs = 5,
    column.ann = 'stage_refined', 
    column.color = 'stage_refined_color'
  )

  # Find Spearman correlations with pseudotime.
  cellrouter_ <- correlationPseudotime(cellrouter_, type = 'spearman')

  # Get top genes.
  cellrouter_ <- topGenes(cellrouter_, 0.8, 0.1)

  # Smooth dynamics.
  cellrouter_ <- smoothDynamics(cellrouter_, names = unique(names(cellrouter_@pathsinfo$distr)))

  # Custer genes according to pseudotime.
  cellrouter_ <- clusterGenesPseudotime(cellrouter_, 5)

  # Return processed CellRouter object.
  return(cellrouter_)
}

# Find trajectories between the zygote and the late blastocyst stages & plot
# the expression of Stip1 + core pluripotency genes.
cellrouter <- find_trajectories(
  cellrouter_ = cellrouter,
  source = 'Zygote',
  targets = setdiff(as.vector(cellrouter@assays$RNA@sampTab$stage_refined), 'Zygote')
)
plottrajectories(
  object = cellrouter,
  trajectories = paste0('Zygote.Late blastocyst'),
  geneList = c('Stip1', 'Nanog', 'Pou5f1', 'Sox2'),
  rescale = TRUE,
  width = 5,
  height = 3,
  filename = paste0(getwd(), '/results/single_cell/trajectories_Zygote_to_Late_blastocyst.pdf')
)

# Find trajectories between the early blastocyst and the late blastocyst stages & plot
# the expression of Stip1 + core + extended core pluripotency genes.
cellrouter <- find_trajectories(
  cellrouter_ = cellrouter,
  source = 'Early blastocyst',
  target = 'Late blastocyst'
)
plottrajectories(
  object = cellrouter,
  trajectories = paste0('Early blastocyst.Late blastocyst'),
  geneList = c("Stip1", "Sox2", "Pou5f1", "Nanog", "Stat3", "Klf4", "Smad1", "Zfp42", "Nr0b1", "Esrrb"),
  rescale = TRUE,
  width = 5,
  height = 3,
  filename = paste0(getwd(), '/results/single_cell/trajectories_Early_blastocyst_to_Late_blastocyst.pdf')
)

# Find trajectories between the early blastocyst and the intermediate blastocyst stages & plot
# the expression of Stip1 + core + extended core pluripotency genes.
cellrouter <- find_trajectories(
  cellrouter_ = cellrouter,
  source = 'Early blastocyst',
  target = 'Intermediate blastocyst'
)
plottrajectories(
  object = cellrouter,
  trajectories = paste0('Early blastocyst.Intermediate blastocyst'),
  geneList = c("Sox2"),
  rescale = TRUE,
  width = 5,
  height = 3,
  filename = paste0(getwd(), '/results/single_cell/trajectories_Early_blastocyst_to_Intermediate_blastocyst.pdf')
)

# Find trajectories between the intermediate blastocyst and the late blastocyst stages & plot
# the expression of Stip1 + core + extended core pluripotency genes.
find_trajectories(
  cellrouter_ = cellrouter,
  source = 'Intermediate blastocyst',
  target = setdiff(as.vector(cellrouter@assays$RNA@sampTab$stage_refined), 'Intermediate blastocyst')
)
plottrajectories(
  object = cellrouter,
  trajectories = paste0('Intermediate blastocyst.Late blastocyst'),
  geneList = c("Stip1", "Sox2", "Pou5f1", "Nanog", "Stat3", "Klf4", "Smad1", "Zfp42", "Nr0b1", "Esrrb"),
  rescale = TRUE,
  width = 5,
  height = 3,
  filename = paste0(getwd(), '/results/single_cell/trajectories_Intermediate_blastocyst_to_Late_blastocyst.pdf')
)

save(cellrouter, file = paste0(getwd(), "/results/single_cell/cellrouter_object.R"))

####################################################################
#### Correlation between Stip1 and all genes (blastocyst stage) ####
####################################################################

# Subset blastocyst cells.
blastocyst_metadata <- metadata[grepl("blastocyst", metadata$stage_refined) == TRUE,]
blastocyst_counts <- count_matrix %>% dplyr::select(rownames(blastocyst_metadata))

# Create CellRouter object with blastocyst cells.
blastocyst_cellrouter <- CreateCellRouter(
  as.matrix(blastocyst_counts),
  min.genes = 0,
  min.cells = 0,
  is.expr = 0
)

# Normalize & scale blastocyst data.
blastocyst_cellrouter <- Normalize(blastocyst_cellrouter)
blastocyst_cellrouter <- scaleData(blastocyst_cellrouter, blocksize = nrow(blastocyst_cellrouter@assays$RNA@ndata))

# Calculate the Pearson correlation between Stip1 and all genes.
corr <- data.frame(matrix(0, ncol = 4))
colnames(corr) <- c('gene1', 'gene2', 'correlation', 'pvalue')
i <- 1
for(gene in rownames(blastocyst_cellrouter@assays$RNA@ndata)){
  c <- cor.test(
    as.numeric(blastocyst_cellrouter@assays$RNA@ndata[gene,]),
    as.numeric(blastocyst_cellrouter@assays$RNA@ndata['Stip1',]),
    method = "pearson"
  )
  corr[i,'gene1'] <- 'Stip1'
  corr[i,'gene2'] <- gene
  corr[i,'correlation'] <- c$estimate # Correlation coefficient
  corr[i,'pvalue'] <- c$p.value # Statistical significance
  i <- i + 1
}

# Remove genes with NAs and order dataframe by correlation coefficient.
corr <- corr[!is.na(corr$correlation),]
corr <- corr[order(corr$correlation, decreasing = TRUE),]

# Save significant results to output file.
write.csv(corr[corr$pvalue <= 0.05,], paste0(getwd(), "/results/single_cell/blastocyst_corr_Stip1_vs_all_genes_pvalue0.05.csv"))

# Remove Stip1 as obviously Stip1's correlation with itself will be 1.
corr <- corr[corr$gene2 != "Stip1",]
corr$position <- 1:nrow(corr)
corr <- corr %>% dplyr::mutate(status = case_when(
  correlation > 0 & pvalue <= 0.05 ~ "Positive",
  correlation < 0 & pvalue <= 0.05 ~ "Negative",
  pvalue > 0.05 ~ "Non-significant",
))

# Plot correlations.
pdf(paste0(getwd(), "/results/single_cell/corrplot.pdf"), width = 5, height = 3)
ggplot(
  corr,
  aes(x = position, y = correlation, color = status)
) +
  geom_point(size = 0.2) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  ) +
  xlab("Ranked genes") +
  ylab("Pearson correlation coefficient") +
  scale_color_manual(
    values = c("royalblue1", "snow3", "indianred1"),
    name = "Correlation"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = corr[corr$gene2 %in% c("Nanog"),],
    aes(label = gene2),
    size = 2,
    color = "black",
    nudge_x = 6,
    min.segment.length = 0,
    segment.size = 0.2,
    seed = 42
  ) +
  geom_text_repel(
    data = corr[corr$gene2 %in% c("Stat3", "Zfp42"),],
    aes(label = gene2),
    size = 2,
    color = "black",
    nudge_x = -2,
    min.segment.length = 0,
    segment.size = 0.2,
    seed = 42
  ) +
  geom_text_repel(
    data = corr[corr$gene2 %in% c("Klf4"),],
    aes(label = gene2),
    size = 2,
    color = "black",
    nudge_x = 4,
    min.segment.length = 0,
    segment.size = 0.2,
    seed = 42
  ) +
  geom_text_repel(
    data = corr[corr$gene2 %in% c("Esrrb"),],
    aes(label = gene2),
    size = 2,
    color = "black",
    nudge_x = -8,
    min.segment.length = 0,
    segment.size = 0.2,
    seed = 42
  ) +
  geom_text_repel(
    data = corr[corr$gene2 %in% c("Pou5f1"),],
    aes(label = gene2),
    size = 2,
    color = "black",
    nudge_x = 5,
    min.segment.length = 0,
    segment.size = 0.2,
    seed = 42
  ) +
  geom_text_repel(
    data = corr[corr$gene2 %in% c("Sox2"),],
    aes(label = gene2),
    size = 2,
    color = "black",
    nudge_x = 5,
    min.segment.length = 0,
    segment.size = 0.2,
    seed = 42
  ) +
  geom_text_repel(
    data = corr[corr$gene2 %in% c("Smad1", "Nr0b1"),],
    aes(label = gene2),
    size = 2,
    color = "black",
    nudge_x = -3,
    min.segment.length = 0,
    segment.size = 0.2,
    seed = 42
  )
dev.off()

# Create ranked gene list.
ranked_genelist <- corr$correlation
names(ranked_genelist) <- corr$gene2
ranked_genelist <- sort(ranked_genelist, decreasing = TRUE)
     
# Create gene lists for Over-Representation Analysis (ORA).
genelist_pos <- names(ranked_genelist)[ranked_genelist > 0]
genelist_neg <- names(ranked_genelist)[ranked_genelist < 0]
genelist <- list(
  "Positive" = genelist_pos,
  "Negative" = genelist_neg
)

# Perform ORA using the Gene Ontology (GO) database.
go_res <- compareCluster(
  genelist,
  fun = "enrichGO",
  ont = "ALL",
  keyType = "SYMBOL",
  OrgDb = "org.Mm.eg.db",
  pvalueCutoff = 0.05
)
go_res <- as.data.frame(go_res)

# Save ORA results to output file.
write.csv(
  go_res, 
  paste0(getwd(), "/results/single_cell/blastocyst_ORA_enrichGO_compareCluster_p0.05_results_table.csv")
)

# Terms to highlight on the plot.
terms <- c(
  "DNA replication",
  "G1/S transition of mitotic cell cycle",
  "G2/M transition of mitotic cell cycle",
  "proteasome core complex",
  "proteasomal protein catabolic process",
  "signal transduction by p53 class mediator",
  "p53 binding",
  "intrinsic apoptotic signaling pathway by p53 class mediator",
  "Wnt signaling pathway",
  "histone modification",
  "mRNA processing",
  "DNA repair",
  "stem cell population maintenance",
  "stem cell proliferation",
  "stem cell division",
  "cell fate commitment"
)

# Plot enrichment results.
pdf(paste0(getwd(), "/results/single_cell/blastocyst_ORA_enrichGO_compareCluster_p0.05_selected_terms.pdf"), width = 12, height = 5)
ggplot(
  go_res[go_res$Description %in% terms,],
  aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)), fill = Cluster, size = Count) 
) +
  geom_point(shape = 21, color = "black") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  xlab(expression(-log[10]~"(Adjusted p-value)")) +
  ylab("Gene Ontology (GO) term") +
  ggtitle("") +
  scale_fill_manual(
    values = c("indianred1", "royalblue1"),
    labels = c("Positive correlations", "Negative correlations"),
    name = "Cluster"
  ) +
  scale_size_continuous(range = c(1,5), name = "# Genes") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed")
dev.off()

##################################################################################
#### Correlation between Stip1 and core pluripotency genes (blastocyst stage) ####
##################################################################################

# Look specifically at the correlation between Stip1 and the core pluripotency genes.
stip1_vs_core <- blastocyst_cellrouter@assays$RNA@ndata[c("Stip1", "Nanog", "Pou5f1", "Sox2"),]
stip1_vs_core <- as.data.frame(t(stip1_vs_core))

# Plot Stip1 vs Sox2.
pdf(paste0(getwd(), "/results/single_cell/blastocyst_Stip1_vs_Sox2.pdf"), width = 4, height = 3)
ggscatter(
  stip1_vs_core, 
  x = "Stip1", 
  y = "Sox2",
  add = "reg.line",
  add.params = list(color = "blue", fill = "lightgray"),
  conf.int = TRUE
) + 
  stat_cor(method = "pearson", label.x = 1, label.y = 1.35) + 
  font("xlab", face = "italic") +
  font("ylab", face = "italic")
dev.off()

# Plot Stip1 vs Nanog.
pdf(paste0(getwd(), "/results/single_cell/blastocyst_Stip1_vs_Nanog.pdf"), width = 4, height = 3)
ggscatter(
  stip1_vs_core, 
  x = "Stip1", 
  y = "Nanog",
  add = "reg.line",
  add.params = list(color = "blue", fill = "lightgray"),
  conf.int = TRUE
) + 
  stat_cor(method = "pearson", label.x = 1, label.y = 1) + 
  font("xlab", face = "italic") +
  font("ylab", face = "italic")
dev.off()

# Plot Stip1 vs Pou5f1.
pdf(paste0(getwd(), "/results/single_cell/blastocyst_Stip1_vs_Pou5f1.pdf"), width = 4, height = 3)
ggscatter(
  stip1_vs_core, 
  x = "Stip1", 
  y = "Pou5f1",
  add = "reg.line",
  add.params = list(color = "blue", fill = "lightgray"),
  conf.int = TRUE
) + 
  stat_cor(method = "pearson", label.x = 1, label.y = 2.6) + 
  font("xlab", face = "italic") +
  font("ylab", face = "italic")
dev.off()

###################################################################################
#### Correlation between Stip1 and core pluripotency genes (across all stages) ####
###################################################################################

stip1_vs_pluripotency_per_stage <- data.frame()
for(stage in cellrouter@assays$RNA@sampTab$stage_refined %>% unique) {

  # Get normalized counts of Stip1 and core/extended core pluripotency genes.
  gene_exp <- cellrouter@assays$RNA@ndata[rownames(cellrouter@assays$RNA@ndata) %in% c("Stip1", "Sox2", "Pou5f1", "Nanog", "Stat3", "Klf4", "Smad1", "Zfp42", "Nr0b1", "Esrrb"),]
  
  # Keep only cells that belong to a given developmental stage.
  gene_exp <- gene_exp[,colnames(gene_exp) %in% rownames(cellrouter@assays$RNA@sampTab)[cellrouter@assays$RNA@sampTab$stage_refined == stage]]

  for(gene in rownames(gene_exp)){
    c <- cor.test(
      as.numeric(gene_exp[gene,]),
      as.numeric(gene_exp["Stip1",]),
      method = "pearson"
    )
    stip1_vs_pluripotency_per_stage <- rbind(
      stip1_vs_pluripotency_per_stage, c(stage, "Stip1", gene, c$estimate, c$p.value)
    )
  }
}
colnames(stip1_vs_pluripotency_per_stage) <- c("stage_refined", "stip1", "gene2", "corr", "pvalue")

stip1_vs_pluripotency_per_stage$stage_refined <- factor(
  stip1_vs_pluripotency_per_stage$stage_refined, 
  levels = c("Zygote", "2-cell (C57)", "Early 2-cell", "Intermediate 2-cell", "Late 2-cell", "4-cell", "8-cell",
             "16-cell", "Early blastocyst", "Intermediate blastocyst", "Late blastocyst")
)

stip1_vs_pluripotency_per_stage$significant <- ifelse(stip1_vs_pluripotency_per_stage$pvalue <= 0.05, TRUE, FALSE)
stip1_vs_pluripotency_per_stage$significant[stip1_vs_pluripotency_per_stage$corr == 1] <- TRUE

pdf(
  paste0(getwd(), "/results/single_cell/correlation_Stip1_vs_pluripotency_genes_across_stages.pdf"), 
  width = 12, 
  height = 4
)
ggplot(
  stip1_vs_pluripotency_per_stage[is.na(stip1_vs_pluripotency_per_stage$pvalue) == FALSE & stip1_vs_pluripotency_per_stage$gene2 != "Stip1",],
  aes(x = stage_refined, y = as.numeric(corr))
) +
  geom_point(shape = 21, aes(color = significant)) +
  geom_line(aes(group = gene2)) +
  facet_wrap(~gene2, nrow = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text.x = element_text(size = 12, face = "italic")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("snow3", "firebrick2"),
    name = "Significant"
  ) +
  xlab("Pre-implantation stage") +
  ylab("Pearson correlation")
dev.off()

############################################################################
#### Correlation between Stip1 and heat shock genes (across all stages) ####
############################################################################

stip1_vs_heatShock_per_stage <- data.frame()
for(stage in cellrouter@assays$RNA@sampTab$stage_refined %>% unique) {

  # Get normalized counts of Stip1 and core/extended core pluripotency genes.
  gene_exp <- cellrouter@assays$RNA@ndata[rownames(cellrouter@assays$RNA@ndata) %in% c("Stip1", "Hspa4l", "Hspa13", "Hspa14", "Hspa9", "Hspa2", "Hspa4", "Hsp90aa1", "Hspa1a", "Hspa1b", "Hsp90ab1", "Hspa5", "Hspa12a", "Hsp90b1", "Hspa8"),]
  
  # Keep only cells that belong to a given developmental stage.
  gene_exp <- gene_exp[,colnames(gene_exp) %in% rownames(cellrouter@assays$RNA@sampTab)[cellrouter@assays$RNA@sampTab$stage_refined == stage]]

  for(gene in rownames(gene_exp)){
    c <- cor.test(
      as.numeric(gene_exp[gene,]),
      as.numeric(gene_exp["Stip1",]),
      method = "pearson"
    )
    stip1_vs_heatShock_per_stage <- rbind(
      stip1_vs_heatShock_per_stage, c(stage, "Stip1", gene, c$estimate, c$p.value)
    )
  }
}
colnames(stip1_vs_heatShock_per_stage) <- c("stage_refined", "stip1", "gene2", "corr", "pvalue")

stip1_vs_heatShock_per_stage$stage_refined <- factor(
  stip1_vs_heatShock_per_stage$stage_refined, 
  levels = c("Zygote", "2-cell (C57)", "Early 2-cell", "Intermediate 2-cell", "Late 2-cell", "4-cell", "8-cell",
             "16-cell", "Early blastocyst", "Intermediate blastocyst", "Late blastocyst")
)

stip1_vs_heatShock_per_stage$significant <- ifelse(stip1_vs_heatShock_per_stage$pvalue <= 0.05, TRUE, FALSE)
stip1_vs_heatShock_per_stage$significant[stip1_vs_heatShock_per_stage$corr == 1] <- TRUE

pdf(
  paste0(getwd(), "/results/single_cell/correlation_Stip1_vs_heatShock_genes_across_stages.pdf"), 
  width = 14, 
  height = 4
)
ggplot(
  stip1_vs_heatShock_per_stage[is.na(stip1_vs_heatShock_per_stage$pvalue) == FALSE & stip1_vs_heatShock_per_stage$gene2 != "Stip1",],
  aes(x = stage_refined, y = as.numeric(corr))
) +
  geom_point(shape = 21, aes(color = significant)) +
  geom_line(aes(group = gene2)) +
  facet_wrap(~gene2, nrow = 2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text.x = element_text(size = 12, face = "italic")
  ) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("snow3", "firebrick2"),
    name = "Significant"
  ) +
  xlab("Pre-implantation stage") +
  ylab("Pearson correlation")
dev.off()