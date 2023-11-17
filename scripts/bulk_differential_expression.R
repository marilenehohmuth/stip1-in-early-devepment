# ANALYSIS OF IN-HOUSE RNA-SEQ DATA ##########################################
# @ This script generates the plots present in the second main figure of the # 
# @ manuscript + supplemental material.                                      #
##############################################################################


##########################
#### Loading packages ####
##########################

# R v4.0.3
library(DESeq2) # v1.30.1
library(dplyr) # v1.0.10 
library(ggplot2) # v3.3.6
library(RColorBrewer) # v1.1_3
library(ggrepel) # v0.9.1


#################################
#### Defining some functions ####
#################################

find_degs <- function(
     counts, # Count matrix.
     metadata, # Metadata table.
     analysis_design, # Design for DESeq2 analysis (e.g. ~genotype).
     column, # Name of column that contains the conditions to be compared.
     condition_of_interest, # Name of condition of interest (e.g. TGA).
     control_condition # Name of control condition (e.g. WT).
) {
     # Create DESeq2 object.
     dds <- DESeqDataSetFromMatrix(
          countData = counts,
          colData = metadata,
          design = analysis_design
     )

     # Filter out genes that have very low expression across samples.
     dds <- dds[rowSums(counts(dds)) > 1,]

     # Run differential gene expression analysis.
     dds <- estimateSizeFactors(dds)
     dds <- DESeq(dds)

     # Save DESeq2-normalized counts to output file.
     write.csv(
          counts(dds, normalized = TRUE), 
          paste0(getwd(), "/results/bulk/", condition_of_interest, "_vs_", control_condition, "/DESeq2-normalized_counts.csv")
     )

     # Extract results with the contrast of interest.
     res_unshrunken <- results(
          dds, 
          contrast = c(column, condition_of_interest, control_condition), 
          alpha = 0.05
     )

     # Make plot with unshrunken results & save it to output file.
     pdf(
          paste0(getwd(), "/results/bulk/", condition_of_interest, "_vs_", control_condition, "/DESeq2_MAplot_unshrunken_results.pdf"),
          height = 6, 
          width = 7
     )
     plotMA(res_unshrunken)
     dev.off()

     # Perform log fold change shrinkage.
     res_shrunken <- lfcShrink(
          dds = dds, 
          res = res_unshrunken, 
          contrast = c(column, condition_of_interest, control_condition), 
          type = "normal"
     )

     # Make plot with shrunken results & save it to output file.
     pdf(
          paste0(getwd(), "/results/bulk/", condition_of_interest, "_vs_", control_condition, "/DESeq2_MAplot_shrunken_results.pdf"),
          height = 6, 
          width = 7
     )
     plotMA(res_shrunken)
     dev.off()

     res_shrunken <- as.data.frame(res_shrunken)

     # Save unfiltered results to output file.
     write.csv(
          res_shrunken, 
          paste0(getwd(), "/results/bulk/", condition_of_interest, "_vs_", control_condition, "/DESeq2_unfiltered_DEGs.csv")
     )

     # Create dataframe with only significant results.
     deg_sig <- res_shrunken[res_shrunken$pvalue <= 0.05,]
     
     # Save filtered results to output file.
     write.csv(
          deg_sig, 
          paste0(getwd(), "/results/bulk/", condition_of_interest, "_vs_", control_condition, "/DESeq2_filtered_DEGs_pvalue0.05.csv")
     )

     return(res_shrunken)
}

make_volcano <- function(
     degs, # DESeq2 output.
     lfc_threshold, # Log fold change threshold.
     genelist, # List of genes to be highlighted.
     title, # Title for volcano plot.
     nudge_x_ # Horizontal adjustment for point labels.
) {

     # Classify genes.
     degs <- degs %>% mutate(status = case_when(
          pvalue <= 0.05 & log2FoldChange >= lfc_threshold ~ "Upregulated",
          pvalue <= 0.05 & log2FoldChange <= -lfc_threshold ~ "Downregulated",
          pvalue > 0.05 ~ "Non-significant",
          abs(log2FoldChange) < lfc_threshold ~ "Non-significant",
     ))

     # Subset genes to be highlighted on the volcano plot.
     degs_to_highlight <- degs[rownames(degs) %in% genelist,]

     # Create volcano plot.
     volcano <- ggplot(
          degs,
          aes(x = log2FoldChange, y = -log10(pvalue), color = status)
     ) +
          geom_point(alpha = 0.5) +
          theme_bw() +
          theme(
               axis.text = element_text(size = 10),
               axis.title = element_text(size = 12),
               legend.text = element_text(size = 10),
               legend.title = element_text(size = 12),
               plot.title = element_text(size = 12, face = 'bold', hjust = 0.5)
          ) +
          scale_color_manual(
               name = '',
               values = c('royalblue1', 'snow3', 'indianred1'),
               labels = c(
                    paste0('Downregulated (n=', nrow(degs[degs$status == "Downregulated",]), ')'),
                    paste0('Non-significant (n=', nrow(degs[degs$status == "Non-significant",]), ')'),
                    paste0('Upregulated (n=', nrow(degs[degs$status == "Upregulated",]), ')')
               )
          ) +
          geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = 'dashed') +
          geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
          geom_text_repel(
               data = degs_to_highlight %>% filter(status == "Downregulated"),
               aes(label = degs_to_highlight %>% filter(status == "Downregulated") %>% rownames),
               show.legend = FALSE,
               size = 2,
               color = 'black',
               nudge_x = -nudge_x_
          ) +
          geom_text_repel(
               data = degs_to_highlight %>% filter(status == "Upregulated"),
               aes(label = degs_to_highlight %>% filter(status == "Upregulated") %>% rownames),
               show.legend = FALSE,
               size = 2,
               color = 'black',
               nudge_x = nudge_x_
          ) +
          xlab(expression(log[2]~'(Fold change)')) +
          ylab(expression(-log[10]~'(p-value)')) +
          ggtitle(title)

     # Return volcano plot.
     return(volcano)
}


#########################################
#### Reading & processing input data ####
#########################################

#### TGA vs WT ####

# Load counts.
oe_counts <- read.csv(
     paste0(getwd(), "/data/TGA/TGA_vs_WT_raw_counts.csv"),
     header = TRUE,
     row.names = 1,
     sep = ";"
)

# Load metadata.
oe_metadata <- read.csv(
     paste0(getwd(), "/data/TGA/TGA_vs_WT_metadata.csv"),
     header = TRUE,
     row.names = 1,
     sep = ";"
)

# Find differentially expressed genes.
oe_degs <- find_degs(
     counts = oe_counts,
     metadata = oe_metadata,
     analysis_design = as.formula(paste("~", "Genotype")),
     column = "Genotype",
     condition_of_interest = "TGA",
     control_condition = "WT"
)

# Make volcano plot.
pdf(
     paste0(getwd(), "/results/bulk/TGA_vs_WT/DESeq2_volcano_plot.pdf"),
     width = 5, 
     height = 3
)
make_volcano(
     degs = oe_degs,
     lfc_threshold = 0,
     genelist = c("Stip1", "Lefty2", "Mid1", "Gata6", "Klf2", "Bmp4", "Cdkn1b", "Anxa1", "Hspa12a", "Tfcp2l1"),
     title = expression(STIP1^'TgA'~italic('versus')~'WT'),
     nudge_x_ = 0.5
)
dev.off()


#### HTKO vs WT ####

# Read counts.
htko_counts <- read.csv(
     paste0(getwd(), "/data/HTKO/HTKO_vs_WT_raw_counts.csv"),
     header = TRUE,
     row.names = 1,
     sep = ";"
)

# Read metadata.
htko_metadata <- read.csv(
     paste0(getwd(), "/data/HTKO/HTKO_vs_WT_metadata.csv"),
     header = TRUE,
     row.names = 1,
     sep = ";"
)

# Find differentially expressed genes.
htko_degs <- find_degs(
     counts = htko_counts,
     metadata = htko_metadata,
     analysis_design = as.formula(paste("~", "Genotype")),
     column = "Genotype",
     condition_of_interest = "HTKO",
     control_condition = "WT"
)

# Make volcano plot.
pdf(
     paste0(getwd(), "/results/bulk/HTKO_vs_WT/DESeq2_volcano_plot.pdf"),
     width = 5, 
     height = 3
)
make_volcano(
     degs = htko_degs,
     lfc_threshold = 0,
     genelist = c("Stip1", "Usp11", "Foxo3", "Hspa1b"),
     title = expression(STIP1^'+/-'~italic('versus')~'WT'),
     nudge_x_ = 0
)
dev.off()


#### TGA vs HTKO ####

# Read counts.
oe_counts <- read.csv(
     paste0(getwd(), "/data/TGA/TGA_vs_HTKO_raw_counts.csv"),
     header = TRUE,
     row.names = 1,
     sep = ";"
)

# Read metadata.
oe_metadata <- read.csv(
     paste0(getwd(), "/data/TGA/TGA_vs_HTKO_metadata.csv"),
     header = TRUE,
     row.names = 1,
     sep = ";"
)

# Find differentially expressed genes.
oe_degs <- find_degs(
     counts = oe_counts,
     metadata = oe_metadata,
     analysis_design = as.formula(paste("~", "Genotype")),
     column = "Genotype",
     condition_of_interest = "TGA",
     control_condition = "HTKO"
)

# Make volcano plot.
pdf(
     paste0(getwd(), "/results/bulk/TGA_vs_HTKO/DESeq2_volcano_plot.pdf"),
     width = 5, 
     height = 3
)
make_volcano(
     degs = oe_degs,
     lfc_threshold = 0,
     genelist = c("Stip1", "Hspa1b"),
     title = expression(STIP1^'TgA'~italic('versus')~STIP1^'+/-'),
     nudge_x_ = 0.5
)
dev.off()


#### TPR1 vs WT ####

# Read counts.
delta_counts <- read.csv(
     paste0(getwd(), "/data/TPR1/TPR1_vs_WT_raw_counts.csv"),
     header = TRUE,
     row.names = 1,
     sep = ";"
)

# Read metadata.
delta_metadata <- read.csv(
     paste0(getwd(), "/data/TPR1/TPR1_vs_WT_metadata.csv"),
     header = TRUE,
     row.names = 1,
     sep = ";"
)

# Find differentially expressed genes.
delta_degs <- find_degs(
     counts = delta_counts,
     metadata = delta_metadata,
     analysis_design = as.formula(paste("~", "Genotype")),
     column = "Genotype",
     condition_of_interest = "TPR1",
     control_condition = "WT"
)

# Make volcano plot.
pdf(
     paste0(getwd(), "/results/bulk/TPR1_vs_WT/DESeq2_volcano_plot.pdf"),
     width = 5, 
     height = 3
)
make_volcano(
     degs = delta_degs,
     lfc_threshold = 0,
     genelist = c("Stip1", "Lefty1", "Lefty2", "Klf2", "Foxo4", "Otx2", "Gata4", "Gata6", "Fgfr4", "Hspa1a", "Hspa1b"),
     title = expression(STIP1^{Delta~'TPR1'}~italic('versus')~'WT'),
     nudge_x_ = 0.5
)
dev.off()