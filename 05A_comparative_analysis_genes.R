# --------------------
# Import packages
# --------------------
library(data.table)
library(readr)
library(dplyr)
library(ggplot2)


# --------------------
# Input variables
# --------------------
mitocarta <- "after_mitocarta"

DESEQ_DIR <- "./data/05_deseq"
DATA_OUT_DIR <- "./data/06_comparative_analysis/gene-level"
OUT_DIR <- "./figuras/comparative_analysis/gene-level"

comparisons <- list(
  # 1) Young adult males vs females
  c("SKA111", "G5_vs_G9", "male_vs_female", "M_20_39_vs_F_20_39"),
  c("SKA113", "G7_vs_G11", "male_vs_female", "M_20_39_vs_F_20_39"),
  
  # 2) Old males vs females
  c("SKA111", "G1_vs_G2", "male_vs_female", "M_50_59_vs_F_50_59"),
  c("SKA113", "G3_vs_G4", "male_vs_female", "M_50_59_vs_F_50_59"),
  c("SKA111", "G1_vs_G2", "male_vs_female", "M_60_79_vs_F_60_79"),
  c("SKA113", "G3_vs_G4", "male_vs_female", "M_60_79_vs_F_60_79"),
  
  # 3) Males ORQ vs SHAM
  c("SKA111", "G6_vs_G5", "male", "M_40_49_vs_M_20_39"),
  # c("SKA111", "G6_vs_G5", "male", "M_50_59_vs_M_20_39"),
  # c("SKA111", "G6_vs_G5", "male", "M_60_79_vs_M_20_39"),
  c("SKA113", "G8_vs_G7", "male", "M_40_49_vs_M_20_39"),
  # c("SKA113", "G8_vs_G7", "male", "M_50_59_vs_M_20_39"),
  # c("SKA113", "G8_vs_G7", "male", "M_60_79_vs_M_20_39"),
  
  # 4) Females OVX 1 month vs SHAM
  c("SKA111", "G10_vs_G9", "female", "F_40_49_vs_F_20_39"),
  c("SKA113", "G12_vs_G11", "female", "F_40_49_vs_F_20_39"),
  
  # 5) Females OVX 7 months vs SHAM
  c("SKA111", "G14_vs_G13", "female", "F_50_59_vs_F_20_39"),
  c("SKA111", "G14_vs_G13", "female", "F_60_79_vs_F_20_39")
)

for(comp in comparisons) {
  genotype_mouse <- comp[1]
  group_mouse <- comp[2]
  
  genotype_human <- comp[3]
  group_human <- comp[4]
  

  # --------------------
  # Read DESeq2 tables
  # --------------------
  # Mouse
  res_mouse <- read.csv(file.path(DESEQ_DIR, "mouse", 
                                  "liver", mitocarta,
                                  genotype_mouse,
                                  paste0("deseq_all_", group_mouse, ".tsv")), 
                        header=TRUE, sep='\t')
  
  # Human
  res_human <- read.csv(file.path(DESEQ_DIR, "human", 
                                  "liver", mitocarta,
                                  genotype_human,
                                  paste0("deseq_all_", group_human, ".tsv")), 
                        header=TRUE, sep='\t')
  
  # Remove version from human Ensembl ID
  res_human <- res_human %>%
    mutate(
      Human_gene_stable_ID = gsub("\\..*", "", ensembl)
    )
  
  
  # --------------------
  # Define DEGs
  # --------------------
  # Threshold padj < 0.05 and |log2FC| > 0.379
  
  # Mouse
  deg_mouse <- res_mouse %>%
    filter(!is.na(padj),
           padj < 0.05,
           abs(log2FoldChange) > 0.379)
  
  # Human
  deg_human <- res_human %>%
    filter(!is.na(padj),
           padj < 0.05,
           abs(log2FoldChange) > 0.379)
  
  
  # --------------------
  # MOUSE 2 HUMAN ortholog table
  # --------------------
  mouse2human <- read.delim("./data/mouse2human.txt")
  mouse2human[mouse2human == ""] <- NA
  
  mouse2human <- mouse2human %>%
    dplyr::rename(
      Mouse_gene_stable_ID  = Gene.stable.ID,
      Mouse_gene_name       = Gene.name,
      Human_gene_stable_ID  = Human.gene.stable.ID,
      Human_gene_name       = Human.gene.name
    )
  
  # Join DESeq2 mouse table with ortholog information (using mouse symbol)
  deg_mouse_annot <- deg_mouse %>%
    left_join(mouse2human,
              by = c("symbol" = "Mouse_gene_name"))
  
  # Check how many mouse DEGs do NOT have a human ortholog
  # sum(is.na(deg_mouse_annot$Human_gene_stable_ID))
  
  
  # --------------------
  # Build DEG lists in "human gene space"
  # --------------------
  # Mouse DEGs mapped to human IDs, by direction
  mouse_up <- deg_mouse_annot %>%
    filter(log2FoldChange > 0,
           !is.na(Human_gene_stable_ID)) %>%
    pull(Human_gene_stable_ID) %>%
    unique()
  
  mouse_down <- deg_mouse_annot %>%
    filter(log2FoldChange < 0,
           !is.na(Human_gene_stable_ID)) %>%
    pull(Human_gene_stable_ID) %>%
    unique()
  
  # Human DEGs (already in human ID space), by direction
  human_up <- deg_human %>%
    filter(log2FoldChange > 0) %>%
    pull(Human_gene_stable_ID) %>%
    unique()
  
  human_down <- deg_human %>%
    filter(log2FoldChange < 0) %>%
    pull(Human_gene_stable_ID) %>%
    unique()
  
  
  # --------------------
  # Direction-specific overlaps
  # --------------------
  # Up in mouse, up in human (conserved activation)
  up_up <- intersect(mouse_up, human_up)
  
  # Down in mouse, down in human (conserved repression)
  down_down <- intersect(mouse_down, human_down)
  
  # Up in mouse, down in human (opposite response)
  up_down <- intersect(mouse_up, human_down)
  
  # Down in mouse, up in human (opposite response)
  down_up <- intersect(mouse_down, human_up)
  
  
  # --------------------
  # Build overlap tables with human gene symbols (for interpretation)
  # --------------------
  # Map human IDs to human symbols
  human_id2symbol <- res_human %>%
    dplyr::select(Human_gene_stable_ID, symbol_human = symbol)
  
  # Helper function: get a data frame for a given overlap
  get_overlap_df <- function(ids_vector, label) {
    tibble(Human_gene_stable_ID = ids_vector) %>%
      left_join(human_id2symbol, by = "Human_gene_stable_ID") %>%
      mutate(overlap_type = label)
  }
  
  overlap_up_up_df <- get_overlap_df(up_up, "mouse_up & human_up")
  overlap_down_down_df <- get_overlap_df(down_down, "mouse_down & human_down")
  overlap_up_down_df <- get_overlap_df(up_down, "mouse_up & human_down")
  overlap_down_up_df <- get_overlap_df(down_up, "mouse_down & human_up")
  
  # Bind all overlaps into one table
  all_overlaps_df <- bind_rows(
    overlap_up_up_df,
    overlap_down_down_df,
    overlap_up_down_df,
    overlap_down_up_df
  )
  
  
  # --------------------
  # Species-specific DEGs (optional, simple version)
  # --------------------
  # Mouse DEGs mapped to human IDs (any direction)
  mouse_deg_human_ids <- unique(c(mouse_up, mouse_down))
  
  # Human DEGs (any direction)
  human_deg_ids <- unique(c(human_up, human_down))
  
  # Human-specific DEGs: DEGs in human with no ortholog DEG in mouse
  human_specific_ids <- setdiff(human_deg_ids, mouse_deg_human_ids)
  
  # Mouse-specific DEGs (in human ID space): mapped mouse DEGs that are not DEG in human
  mouse_specific_ids <- setdiff(mouse_deg_human_ids, human_deg_ids)
  
  # Add symbols
  human_specific_df <- tibble(Human_gene_stable_ID = human_specific_ids) %>%
    left_join(human_id2symbol, by = "Human_gene_stable_ID")
  
  mouse_specific_df <- tibble(Human_gene_stable_ID = mouse_specific_ids) %>%
    left_join(human_id2symbol, by = "Human_gene_stable_ID")
  
  
  # --------------------
  # Save data
  # --------------------
  dir.create(file.path(DATA_OUT_DIR, genotype_human), showWarnings = FALSE, recursive = TRUE)
  
  # Save overlaps table
  write.table(
    all_overlaps_df,
    file = file.path(DATA_OUT_DIR, genotype_human,
                     paste0("mouse_human_DEG_overlaps_", 
                            group_mouse, "_AND_", group_human, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Save human-specific DEGs
  write.table(
    human_specific_df,
    file = file.path(DATA_OUT_DIR, genotype_human,
                     paste0("human_specific_DEGs_", 
                            group_mouse, "_AND_", group_human, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Save mouse-specific DEGs (mapped to human IDs)
  write.table(
    mouse_specific_df,
    file = file.path(DATA_OUT_DIR, genotype_human,
                     paste0("mouse_specific_DEGs_", 
                            group_mouse, "_AND_", group_human, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  dir.create(file.path(OUT_DIR, genotype_human), showWarnings = FALSE, recursive = TRUE)
  
  if (nrow(all_overlaps_df) > 0) {
    # Count genes in each overlap category
    overlap_counts <- all_overlaps_df %>%
      group_by(overlap_type) %>%
      summarise(n_genes = n(), .groups = "drop")
    
    p_bar_overlap <- ggplot(overlap_counts,
                            aes(x = overlap_type, y = n_genes, fill = overlap_type)) +
      geom_col(width = 0.5) +
      theme_minimal() +
      labs(
        title = paste0("Overlapping orthologous DEGs:\n", group_mouse, " and ", group_human),
        x = "Overlap category",
        y = "Number of genes"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 15),
            legend.position = "none")
    
    ggsave(
      filename = file.path(OUT_DIR, genotype_human,
                           paste0("barplot_DEG_overlap_",
                                  group_mouse, "_AND_", group_human, ".png")),
      plot = p_bar_overlap,
      width = 7, height = 5, dpi = 300
    )
  }
  
  
  if (nrow(all_overlaps_df) > 0 || nrow(mouse_specific_df) > 0 || nrow(human_specific_df) > 0) {
    n_shared <- length(unique(all_overlaps_df$Human_gene_stable_ID))
    n_mouse_sp <- nrow(mouse_specific_df)
    n_human_sp <- nrow(human_specific_df)
    
    shared_specific_counts <- tibble(
      category = c("Shared DEGs", "Mouse specific DEGs", "Human specific DEGs"),
      n_genes  = c(n_shared, n_mouse_sp, n_human_sp)
    )
    
    p_bar_shared <- ggplot(shared_specific_counts,
                           aes(x = category, y = n_genes, fill = category)) +
      geom_col(width = 0.5) +
      theme_minimal() +
      labs(
        title = paste0("Shared and species specific orthologous DEGs:\n", group_mouse, " and ", group_human),
        x = "",
        y = "Number of genes"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 15),
            legend.position = "none")
    
    ggsave(
      filename = file.path(OUT_DIR, genotype_human,
                           paste0("barplot_shared_vs_specific_DEGs_",
                                  group_mouse, "_AND_", group_human, ".png")),
      plot   = p_bar_shared,
      width  = 7, height = 5, dpi = 300
    )
  }
}
