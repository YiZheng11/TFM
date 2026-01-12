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

GSEA_DIR <- "./extra"
DATA_OUT_DIR <- "./data/06_comparative_analysis/functional-level"
OUT_DIR <- "./figuras/comparative_analysis/functional-level"

comparisons <- list(
  # 1) Young adult males vs females
  c("SKA111", "G5_vs_G9", "male_vs_female", "M_20_39_vs_F_20_39"),
  c("SKA111", "G5_vs_G9", "male_vs_female", "M_40_49_vs_F_40_49"),
  c("SKA113", "G7_vs_G11", "male_vs_female", "M_20_39_vs_F_20_39"),
  c("SKA113", "G7_vs_G11", "male_vs_female", "M_40_49_vs_F_40_49"),
  
  # 2) Old males vs females
  c("SKA111", "G1_vs_G2", "male_vs_female", "M_50_59_vs_F_50_59"),
  # c("SKA113", "G3_vs_G4", "male_vs_female", "M_50_59_vs_F_50_59"),
  c("SKA111", "G1_vs_G2", "male_vs_female", "M_60_79_vs_F_60_79"),
  # c("SKA113", "G3_vs_G4", "male_vs_female", "M_60_79_vs_F_60_79"),
  
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
  gsea_mouse <- read.csv(file.path(GSEA_DIR, "mouse", 
                                  "liver", mitocarta,
                                  genotype_mouse, "gsea",
                                  paste0("GSEA_final_terms_", group_mouse, ".tsv")), 
                        header=TRUE, sep='\t')
  
  # Human
  gsea_human <- read.csv(file.path(GSEA_DIR, "human", 
                                  "liver", mitocarta,
                                  genotype_human, "gsea",
                                  paste0("GSEA_final_terms_", group_human, ".tsv")), 
                        header=TRUE, sep='\t')
  
  
  gsea_both <- gsea_mouse %>%
    dplyr::select(ID,
                  Description_mouse = Description,
                  NES_mouse = NES,
                  padj_mouse = p.adjust) %>%
    full_join(
      gsea_human %>%
        dplyr::select(ID,
                      Description_human = Description,
                      NES_human = NES,
                      padj_human = p.adjust),
      by = "ID"
    )
  
  gsea_both <- gsea_both %>%
    mutate(
      pattern = case_when(
        !is.na(padj_mouse) & !is.na(padj_human) &
          padj_mouse < 0.05 & padj_human < 0.05 & NES_mouse > 0 & NES_human > 0 ~ "both_up",
        !is.na(padj_mouse) & !is.na(padj_human) &
          padj_mouse < 0.05 & padj_human < 0.05 & NES_mouse < 0 & NES_human < 0 ~ "both_down",
        !is.na(padj_mouse) & !is.na(padj_human) &
          padj_mouse < 0.05 & padj_human < 0.05 & NES_mouse > 0 & NES_human < 0 ~ "mouse_up_human_down",
        !is.na(padj_mouse) & !is.na(padj_human) &
          padj_mouse < 0.05 & padj_human < 0.05 & NES_mouse < 0 & NES_human > 0 ~ "mouse_down_human_up",
        !is.na(padj_mouse) & is.na(padj_human) &
          padj_mouse < 0.05 ~ "mouse_specific",
        is.na(padj_mouse) & !is.na(padj_human) &
          padj_human < 0.05 ~ "human_specific",
        TRUE ~ "not_significant_or_missing"
      )
    )
  
  # 4) Species-specific functional terms (optional)
  mouse_specific_terms <- gsea_both %>%
    filter(pattern == "mouse_specific")
  
  human_specific_terms <- gsea_both %>%
    filter(pattern == "human_specific")
  
  shared_terms <- gsea_both %>%
    filter(pattern %in% c("both_up", "both_down",
                          "mouse_up_human_down", "mouse_down_human_up"))
  
  # Conserved activated pathways
  conserved_up <- gsea_both %>%
    filter(pattern == "both_up") %>%
    arrange(padj_mouse, padj_human)
  
  # Conserved repressed pathways
  conserved_down <- gsea_both %>%
    filter(pattern == "both_down") %>%
    arrange(padj_mouse, padj_human)
  
  # Opposite response
  opposite <- gsea_both %>%
    filter(pattern %in% c("mouse_up_human_down", "mouse_down_human_up"))
  
  
  # --------------------
  # Save data
  # --------------------
  dir.create(file.path(DATA_OUT_DIR, genotype_human), showWarnings = FALSE, recursive = TRUE)
  
  # Save overlaps table
  write.table(
    gsea_both,
    file = file.path(DATA_OUT_DIR, genotype_human,
                     paste0("mouse_human_GSEA_overlaps_", 
                            group_mouse, "_AND_", group_human, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  
  # --------------------
  # 6) Barplot: counts per pattern
  # --------------------
  tag <- paste0(group_mouse, "_AND_", group_human)
  
  dir.create(file.path(OUT_DIR, genotype_human), showWarnings = FALSE, recursive = TRUE)
  
  if (nrow(gsea_both) > 0) {
    pattern_counts <- gsea_both %>%
      group_by(pattern) %>%
      summarise(n_terms = n(), .groups = "drop")
    
    p_bar_pattern <- ggplot(pattern_counts,
                            aes(x = pattern, y = n_terms, fill = pattern)) +
      geom_col(width = 0.5) +
      theme_minimal() +
      labs(
        title = paste0("GO patterns: ", group_mouse, " and ", group_human),
        x = "Pattern",
        y = "Number of GO terms"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 15),
            legend.position = "none")
    
    ggsave(
      filename = file.path(OUT_DIR, genotype_human,
                           paste0("barplot_GO_patterns_", tag, ".png")),
      plot = p_bar_pattern,
      width = 7, height = 4, dpi = 300
      )
    }
  }