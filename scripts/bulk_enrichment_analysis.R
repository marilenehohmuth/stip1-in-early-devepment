# ANALYSIS OF IN-HOUSE BULK RNA-SEQ DATA ######################################
# @ This script generates plots present in the second and fourth main figures # 
# @ of the manuscript + supplemental material.                                #
###############################################################################
#                      #
# Functional profiling #
#                      #
########################


##########################
#### Loading packages ####
##########################

# R v4.0.3
library(clusterProfiler) # v3.18.1
library(org.Mm.eg.db) # v3.12.0
library(ggplot2) # v3.3.6
library(dplyr) # v1.0.10 


#################################
#### Defining some functions ####
#################################

process_degs <- function(
     deg_df, # Dataframe containing DEGs as output by DESeq2.
     comparison # Conditions that are being compared.
) {
     # Create ranked gene list.
     ranked_genelist <- deg_df$log2FoldChange
     names(ranked_genelist) <- deg_df$X
     ranked_genelist <- sort(ranked_genelist, decreasing = TRUE)
     
     # Create gene lists for Over-Representation Analysis (ORA).
     genelist_pos <- names(ranked_genelist)[ranked_genelist > 0]
     genelist_neg <- names(ranked_genelist)[ranked_genelist < 0]
     genelist <- list(
          "Upregulated" = genelist_pos,
          "Downregulated" = genelist_neg
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

     # Save ORA results to output file.
     write.csv(
          as.data.frame(go_res), 
          paste0(getwd(), "/results/bulk/", comparison, "/", comparison, "_ORA_enrichGO_compareCluster_p0.05_results_table.csv")
     )

     go_res <- as.data.frame(go_res)

     return(go_res)
}

make_dotplot <- function(
     go_res, # Output of compareCluster.
     title, # Title for bar plot.
     terms # Terms to plot
) {
     plot <- ggplot(
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
     ggtitle(title) +
     scale_fill_manual(
          values = c("indianred1", "royalblue1"),
          labels = c("Upregulated", "Downregulated"),
          name = "Cluster"
     ) +
     scale_size_continuous(range = c(1,5), name = "# Genes") +
     geom_vline(xintercept = -log10(0.05), linetype = "dashed")

     return(plot)

}


#########################################
#### Reading & processing input data ####
#########################################

#### HTKO vs WT ####

# Read file.
htko_df <- read.csv(
     paste0(getwd(), "/results/bulk/HTKO_vs_WT/DESeq2_filtered_DEGs_pvalue0.05.csv"), 
     header = TRUE,
     row.names = NULL
)

# Perform over-representation analysis (ORA).
htko_go <- process_degs(
     deg_df = htko_df,
     comparison = "HTKO_vs_WT"
)

# Select terms to plot.
htko_terms_to_plot <- c(
     "embryonic placenta development",
     "reproductive system development",
     "epithelial cell proliferation",
     "cellular response to oxidative stress",
     "nuclear transport", 
     "regulation of G1/S transition Of mitotic cell cycle",
     "stress-activated MAPK cascade",
     "trophectodermal cell differentiation",
     "lyase activity",
     "secondary active transmembrane transporter activity",
     "microtubule binding",
     "symporter activity",
     "tubulin binding",
     "chaperone-mediated protein transport",
     "negative regulation of DNA-binding transcription factor activity",
     "mRNA processing",
     "regulation of RNA splicing"
)

# Determine order of clusters.
htko_go$Cluster <- factor(
  htko_go$Cluster, 
  levels = c('Upregulated', 'Downregulated') 
)

# Create plot with ORA results & save it to output file.
htko_dotplot <- make_dotplot(
     go_res = htko_go,
     title = expression(STIP1^'+/-'~italic('versus')~'WT'),
     terms = htko_terms_to_plot
)
pdf(paste0(getwd(), "/results/bulk/HTKO_vs_WT/HTKO_vs_WT_ORA_enrichGO_compareCluster_p0.05_selected_terms.pdf"), width = 10, height = 5)
htko_dotplot
dev.off()

#### TGA vs WT ####

# Read file.
tga_df <- read.csv(
     paste0(getwd(), "/results/bulk/TGA_vs_WT/DESeq2_filtered_DEGs_pvalue0.05.csv"), 
     header = TRUE,
     row.names = NULL
)

# Perform over-representation analysis (ORA).
tga_go <- process_degs(
     deg_df = tga_df,
     comparison = "TGA_vs_WT"
)

# Select terms to plot.
tga_terms_to_plot <- c(
     "regulation of cell cycle arrest",
     "regulation of apoptotic signaling pathway",
     "regulation of cell cycle G1/S phase transition",
     "signal transduction by p53 class mediator",
     "cellular response to radiation",
     "cellular response to hypoxia",
     "endothelial cell migration",
     "cellular response to UV",
     "formation of primary germ layer",
     "nuclear transport",
     "cell fate commitment",
     "signal transduction involved in mitotic G1 DNA damage checkpoint",
     "negative regulation of cell development"
)

# Determine order of clusters.
tga_go$Cluster <- factor(
  tga_go$Cluster, 
  levels = c('Upregulated', 'Downregulated') 
)

# Create plot with ORA results & save it to output file.
tga_dotplot <- make_dotplot(
     go_res = tga_go,
     title = expression(STIP1^'TgA'~italic('versus')~'WT'),
     terms = tga_terms_to_plot
)
pdf(paste0(getwd(), "/results/bulk/TGA_vs_WT/TGA_vs_WT_ORA_enrichGO_compareCluster_p0.05_selected_terms.pdf"), width = 10, height = 5)
tga_dotplot
dev.off()


#### TGA vs HTKO ####

# Read file.
tga_df <- read.csv(
     paste0(getwd(), "/results/bulk/TGA_vs_HTKO/DESeq2_filtered_DEGs_pvalue0.05.csv"), 
     header = TRUE,
     row.names = NULL
)

# Perform over-representation analysis (ORA).
tga_go <- process_degs(
     deg_df = tga_df,
     comparison = "TGA_vs_HTKO"
)

# Select terms to plot.
tga_terms_to_plot <- c(
     "intrinsic apoptotic signaling pathway by p53 class mediator",
     "signal transduction by p53 class mediator",
     "positive regulation of apoptotic signaling pathway",
     "regulation of intrinsic apoptotic signaling pathway",
     "Wnt signaling pathway",
     "regulation of cell cycle G1/S phase transition",
     "placenta development",
     "embryonic placenta development",
     "neural tube development",
     "endoderm development",
     "muscle organ development",
     "cell fate commitment",
     "regulation of muscle cell differentiation"
)

# Determine order of clusters.
tga_go$Cluster <- factor(
  tga_go$Cluster, 
  levels = c('Upregulated', 'Downregulated') 
)

# Create plot with ORA results & save it to output file.
tga_dotplot <- make_dotplot(
     go_res = tga_go,
     title = expression(STIP1^'TgA'~italic('versus')~STIP1^'+/-'),
     terms = tga_terms_to_plot
)
pdf(paste0(getwd(), "/results/bulk/TGA_vs_HTKO/TGA_vs_HTKO_ORA_enrichGO_compareCluster_p0.05_selected_terms.pdf"), width = 11.5, height = 5)
tga_dotplot
dev.off()


#### TPR1 vs WT ####

# Read file.
tpr1_df <- read.csv(
     paste0(getwd(), "/results/bulk/TPR1_vs_WT/DESeq2_filtered_DEGs_pvalue0.05.csv"), 
     header = TRUE,
     row.names = NULL
)

# Perform over-representation analysis (ORA).
tpr1_go <- process_degs(
     deg_df = tpr1_df,
     comparison = "TPR1_vs_WT"
)

# Select terms to plot.
tpr1_terms_to_plot <- c(
     "cell cycle arrest",
     "negative regulation of cell cycle",
     "positive regulation of proteolysis",
     "actin filament organization",
     "regulation of apoptotic signaling pathway",
     "negative regulation of nervous system development",
     "DNA repair",
     "regulation of cell morphogenesis involved in differentiation",
     "proteasome-mediated ubiquitin-dependent protein catabolic process",
     "Rho protein signal transduction",
     "autophagy",
     "epithelial tube morphogenesis",
     "morphogenesis of a branching structure",
     "cell fate commitment",
     "DNA-binding transcription repressor activity",
     "epithelial tube morphogenesis"
)

# Determine order of clusters.
tpr1_go$Cluster <- factor(
  tpr1_go$Cluster, 
  levels = c('Upregulated', 'Downregulated') 
)

# Create plot with ORA results & save it to output file.
tpr1_dotplot <- make_dotplot(
     go_res = tpr1_go,
     title = expression(STIP1^{Delta~'TPR1'}~italic('versus')~'WT'),
     terms = tpr1_terms_to_plot
)
pdf(paste0(getwd(), "/results/bulk/TPR1_vs_WT/TPR1_vs_WT_ORA_enrichGO_compareCluster_p0.05_selected_terms.pdf"), width = 10, height = 5)
tpr1_dotplot
dev.off()