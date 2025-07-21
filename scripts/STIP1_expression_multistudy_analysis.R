#-- Load
library("Seurat"); library("GEOquery"); library("reticulate")
scanpy <- import("scanpy")
numpy <- import("numpy")
pandas <- import("pandas")

#-- Analysis
adata <- scanpy$read_h5ad("/dados/marcelo/Marilene_Holmuth_Lopes/adata_PRJEB11202.h5ad")
counts <- t(adata$X)
rownames(counts) <- adata$var_names$to_list()
colnames(counts) <- adata$obs_names$to_list()
meta.data <- adata$obs
counts <- as.matrix(counts)
E_MTAB_3929_sdrf <- read_delim("Marilene_Holmuth_Lopes/E-MTAB-3929.sdrf.txt", 
                               delim = "\t", 
                               escape_double = FALSE, 
                               trim_ws = TRUE)
meta.data <- E_MTAB_3929_sdrf[match(colnames(counts), E_MTAB_3929_sdrf$`Comment[ENA_SAMPLE]`), ]
seurat <- CreateSeuratObject(counts = counts)
seurat$`Characteristics[developmental stage]` <- meta.data$`Characteristics[developmental stage]`
seurat$`Characteristics[treatment]` <- meta.data$`Characteristics[treatment]`
seurat$`Characteristics[inferred lineage]` <- meta.data$`Characteristics[inferred lineage]`
seurat[["RNA"]]$counts[grepl("STIP1", rownames(seurat[["RNA"]]$counts), fixed = TRUE), , drop = FALSE]
seurat <- NormalizeData(seurat)
VlnPlot(seurat, 
        features = "ENSG00000168439-STIP1", 
        group.by = "Characteristics[developmental stage]", 
        split.by = "Characteristics[inferred lineage]")
data_to_plot <- FetchData(seurat, vars = c("ENSG00000168439-STIP1", "Characteristics[developmental stage]", "Characteristics[inferred lineage]"))
colnames(data_to_plot) <- c("expression", "developm_stage", "inferred_lineage")
png("/dados/marcelo/Marilene_Holmuth_Lopes/Petropoulos et al. 2016.png", width = 6, height = 4, units = 'in', res = 300)
ggplot(data_to_plot, aes(x = developm_stage, y = expression, fill = inferred_lineage)) +
  geom_boxplot() +
  geom_point(alpha = 0.10) +
  labs(title = "Expression of STIP1",
       x = "Developmental Stage",
       y = "Expression Level") +
  theme_minimal() +
  facet_grid(. ~ inferred_lineage, scales = "free_x") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

adata <- scanpy$read_h5ad("/dados/marcelo/Marilene_Holmuth_Lopes/adata_PRJNA562548.h5ad")
counts <- t(adata$X)
rownames(counts) <- adata$var_names$to_list()
colnames(counts) <- adata$obs_names$to_list()
meta.data <- adata$obs
counts <- as.matrix(counts)
gse <- getGEO("GSE136447", GSEMatrix = TRUE)
gse <- gse[[1]]
sample_info <- pData(gse)
characteristics_list <- sample_info[, grep("^characteristics_ch", colnames(sample_info))]
counts <- counts[, match(rownames(characteristics_list), colnames(counts))]
meta.data <- characteristics_list
seurat <- CreateSeuratObject(counts = counts)
seurat$characteristics_ch1 <- meta.data$characteristics_ch1
seurat$characteristics_ch1.1 <- meta.data$characteristics_ch1.1
seurat <- NormalizeData(seurat)
data_to_plot <- FetchData(seurat, vars = c("ENSG00000168439-STIP1", "characteristics_ch1", "characteristics_ch1.1"))
colnames(data_to_plot) <- c("expression", "age", "cell_type")
table(data_to_plot$age)
data_to_plot$age <- factor(data_to_plot$age, levels = c("age: embryo invitro day 6", 
                                                        "age: embryo invitro day 7", 
                                                        "age: embryo invitro day 8", 
                                                        "age: embryo invitro day 9", 
                                                        "age: embryo invitro day 10", 
                                                        "age: embryo invitro day 12", 
                                                        "age: embryo invitro day 13.5", 
                                                        "age: embryo invitro day 14"))
png("/dados/marcelo/Marilene_Holmuth_Lopes/Xiang et al. 2020.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(data_to_plot, aes(x = age, y = expression, fill = cell_type)) +
  geom_boxplot() +
  geom_point(alpha = 0.10) +
  labs(title = "Expression of STIP1",
       x = "Developmental Stage",
       y = "Expression Level") +
  theme_minimal() +
  facet_grid(. ~ cell_type) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 70, hjust = 1))
dev.off()

adata <- scanpy$read_h5ad("/dados/marcelo/Marilene_Holmuth_Lopes/adata_PRJEB40781.h5ad")
counts <- t(adata$X)
rownames(counts) <- adata$var_names$to_list()
colnames(counts) <- adata$obs_names$to_list()
meta.data <- adata$obs
counts <- as.matrix(counts)
E_MTAB_9388_sdrf <- read_delim("Marilene_Holmuth_Lopes/E-MTAB-9388.sdrf.txt", 
                               delim = "\t", 
                               escape_double = FALSE, 
                               trim_ws = TRUE)

meta.data <- E_MTAB_9388_sdrf[match(colnames(counts), E_MTAB_9388_sdrf$`Comment[ENA_SAMPLE]`), ]
seurat <- CreateSeuratObject(counts = counts)
seurat$`Characteristics[sampling site]` <- meta.data$`Characteristics[sampling site]`
seurat$`Characteristics[inferred cell type - authors labels]` <- meta.data$`Characteristics[inferred cell type - authors labels]`
seurat[["RNA"]]$counts[grepl("STIP1", rownames(seurat[["RNA"]]$counts), fixed = TRUE), , drop = FALSE]
seurat <- NormalizeData(seurat)
data_to_plot <- FetchData(seurat, vars = c("ENSG00000168439-STIP1", "Characteristics[sampling site]", "Characteristics[inferred cell type - authors labels]"))
colnames(data_to_plot) <- c("expression", "sampling_site", "inferred_cell_type")
png("/dados/marcelo/Marilene_Holmuth_Lopes/Tyser et al. 2021.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(data_to_plot, aes(x = inferred_cell_type, y = expression, fill = sampling_site)) +
  geom_boxplot() +
  geom_point(alpha = 0.10) +
  labs(title = "Expression of STIP1",
       x = "Inferred Cell Type",
       y = "Expression Level") +
  theme_minimal() +
  facet_grid(. ~ sampling_site) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 70, hjust = 1))
dev.off()

GSE134571 <- readRDS("~/Marilene_Holmuth_Lopes/GSE134571_Posterior48h_H9_Amnion_Merged.rds")
GSE134571 <- UpdateSeuratObject(GSE134571)
GSE134571$inferred_cell_type <- ifelse(GSE134571$RNA_snn_res.0.3 == 0, "Transwell-AMLC", 
                                       ifelse(GSE134571$RNA_snn_res.0.3 == 1, "MeLC2", 
                                              ifelse(GSE134571$RNA_snn_res.0.3 == 2, "Human_ES_cell", 
                                                     ifelse(GSE134571$RNA_snn_res.0.3 == 3, "MeLC1", 
                                                            ifelse(GSE134571$RNA_snn_res.0.3 == 4, "hPGCLC", 
                                                                   ifelse(GSE134571$RNA_snn_res.0.3 == 5, "AMLC", NA))))))
GSE134571[["RNA"]]$counts[grepl("STIP1", rownames(GSE134571[["RNA"]]$counts), fixed = TRUE), , drop = FALSE]
data_to_plot <- FetchData(GSE134571, vars = c("STIP1", "old.ident", "inferred_cell_type"))
colnames(data_to_plot) <- c("expression", "identity", "inferred_cell_type")
png("/dados/marcelo/Marilene_Holmuth_Lopes/Zheng et al. 2019.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(data_to_plot, aes(x = inferred_cell_type, y = expression, fill = identity)) +
  geom_boxplot() +
  geom_point(alpha = 0.10) +
  labs(title = "Expression of STIP1",
       x = "Inferred Cell Type",
       y = "Expression Level") +
  theme_minimal() +
  facet_grid(. ~ identity) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 70, hjust = 1))
dev.off()

adata <- scanpy$read_h5ad("/dados/marcelo/Marilene_Holmuth_Lopes/adata_PRJNA720968.h5ad")
counts <- t(adata$X)
rownames(counts) <- adata$var_names$to_list()
colnames(counts) <- adata$obs_names$to_list()
meta.data <- adata$obs
counts <- as.matrix(counts)
gse <- getGEO("GSE171820", GSEMatrix = TRUE)
gse <- gse[[1]]
sample_info <- pData(gse)
characteristics_list <- sample_info[, grep("^characteristics_ch", colnames(sample_info))]
meta.data <- characteristics_list
seurat <- CreateSeuratObject(counts = counts)
seurat$developmental_stage <- meta.data$characteristics_ch1.3
seurat$lineage <- meta.data$characteristics_ch1.1
seurat[["RNA"]]$counts[grepl("STIP1", rownames(seurat[["RNA"]]$counts), fixed = TRUE), , drop = FALSE]
seurat <- NormalizeData(seurat)
data_to_plot <- FetchData(seurat, vars = c("ENSG00000168439-STIP1", "developmental_stage", "lineage"))
colnames(data_to_plot) <- c("expression", "developmental_stage", "lineage")
data_to_plot$developmental_stage <- ifelse(data_to_plot$developmental_stage == "time point: 3 days after stem cells seeded to form cell aggregates", "Day 3", 
                                           ifelse(data_to_plot$developmental_stage == "time point: 4 days after stem cells seeded to form cell aggregates", "Day 4", 
                                                  ifelse(data_to_plot$developmental_stage == "time point: Embryonic day 5", "Day 5", 
                                                         ifelse(data_to_plot$developmental_stage == "time point: Embryonic day 6", "Day 6", 
                                                                ifelse(data_to_plot$developmental_stage == "time point: Embryonic day 7", "Day 7", NA)))))
data_to_plot$developmental_stage <- factor(data_to_plot$developmental_stage, levels = c("Day 3", 
                                                                                        "Day 4", 
                                                                                        "Day 5", 
                                                                                        "Day 6", 
                                                                                        "Day 7"))
png("/dados/marcelo/Marilene_Holmuth_Lopes/Yanagida et al. 2021.png", width = 8, height = 4, units = 'in', res = 300)
ggplot(data_to_plot, aes(x = developmental_stage, y = expression, fill = lineage)) +
  geom_boxplot() +
  geom_point(alpha = 0.10) +
  labs(title = "Expression of STIP1",
       x = "Developmental Stage",
       y = "Expression Level") +
  theme_minimal() +
  facet_wrap(. ~ lineage, nrow = 2) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 70, hjust = 1))
dev.off()
