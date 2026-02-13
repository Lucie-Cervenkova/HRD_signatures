#################
### Hallmarks ###
#################
library(tidyr)
library(stats)

merged_data <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")

signatures <- c("BRCA1ness", "HRD_signature_Peng", "HRD_signature_Walens", "HRD_Beinse_scaled", "HRD_Zhuang_scaled")
hallmark_cols <- grep("^HALLMARK", names(merged_data), value = TRUE)  # Detect hallmark columns
# hallmark_cols <- signatures

cor_results <- data.frame()
pval_results <- data.frame()

# stats_corr_df <- data.frame(hallmark = character(), HRD_signature = character(), correlation = numeric(), p_value = numeric(), CI.low = numeric(), CI.high = numeric())

for (hallmark in hallmark_cols) {
  for (signature in signatures) {
    cor_test <- cor.test(merged_data[[signature]], merged_data[[hallmark]], method = "pearson")  # Change to "spearman" if needed
    cor_results <- rbind(cor_results, data.frame(Hallmark = hallmark, Signature = signature, Correlation = cor_test$estimate))
    pval_results <- rbind(pval_results, data.frame(Hallmark = hallmark, Signature = signature, PValue = cor_test$p.value))

    # stats_corr_df <- rbind(stats_corr_df, data.frame(hallmark = hallmark, HRD_signature = signature, correlation = cor_test$estimate, p_value = cor_test$p.value, CI.low = cor_test$conf.int[1], CI.high = cor_test$conf.int[2]))
  }
}

# corr_df <- stats_corr_df %>% filter(abs(correlation) >= 0.5 & p_value < 0.05)

# altto_df <- corr_df

# ### Save to Excel file
# library(openxlsx)

# # create new workbook
# wb <- createWorkbook()

# addWorksheet(wb, "ALTTO")
# writeData(wb, sheet = "ALTTO", x = altto_df)


# saveWorkbook(wb, "results/HRD/correlations.xlsx", overwrite = TRUE)



pval_results <- pval_results %>%
  mutate(Adj_PValue = p.adjust(PValue, method = "fdr"))

df_plot <- cor_results %>%
  left_join(pval_results, by = c("Hallmark", "Signature")) %>%
  mutate(Size = -log10(Adj_PValue + 1e-10))  # Avoid log(0)

df_plot <- df_plot %>% mutate(Hallmark = gsub("^HALLMARK_", "", Hallmark))  # Remove "HALLMARK_" prefix

ordered_pathways <- c(
  "E2F_TARGETS", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "G2M_CHECKPOINT", "MITOTIC_SPINDLE", "KRAS_SIGNALING",
  "ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE",
  "P53_PATHWAY",
  "APOPTOSIS", "TNFA_SIGNALING_VIA_NFKB",
  "ANGIOGENESIS", "HYPOXIA",
  "EPITHELIAL_MESENCHYMAL_TRANSITION", "TGF_BETA_SIGNALING",
  "ALLOGRAFT_REJECTION", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", 
  "IL2_STAT5_SIGNALING", "IL6_JAK_STAT3_SIGNALING", "INFLAMMATORY_RESPONSE",
  "GLYCOLYSIS", "OXIDATIVE_PHOSPHORYLATION", "FATTY_ACID_METABOLISM", "CHOLESTEROL_HOMEOSTASIS",
  "DNA_REPAIR", "UNFOLDED_PROTEIN_RESPONSE", "REACTIVE_OXIGEN_SPECIES_PATHWAY",
  "COMPLEMENT", "COAGULATION",
  "NOTCH_SIGNALING", "WNT_BETA_CATENIN_SIGNALING", "HEDGEHOG_SIGNALING", "PI3K_AKT_MTOR_SIGNALING", "MTORC1_SIGNALING",
  "ADIPOGENESIS", "APICAL_JUNCTION", "APICAL_SURFACE", "PROTEIN_SECRETION", "MYOGENESIS"
)

df_plot <- df_plot %>% filter(Hallmark %in% ordered_pathways)

df_plot <- df_plot %>% mutate(Signature = ifelse(Signature == "HRD_signature_Peng", "Peng et al.",
                                            ifelse(Signature == "HRD_signature_Walens", "Walens et al.",
                                              ifelse(Signature == "HRD_Beinse_scaled", "Beinse et al.",
                                                ifelse(Signature == "HRD_Zhuang_scaled", "Zhuang et al.", "BRCA1ness")))))
lvls.sigs <- c("BRCA1ness", "Peng et al.", "Walens et al.", "Beinse et al.", "Zhuang et al.")

p <- ggplot(df_plot, aes(x = Signature, y = Hallmark, size = Size, color = Correlation)) +
  geom_point() +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  scale_size(range = c(2, 10)) +  
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.box = "vertical") +
  guides(
    size = guide_legend(title = "-log10(p-adj)", order = 2),
    color = guide_colorbar(title = "Correlation", order = 3, barwidth = 12, barheight = 1) 
  )      
p
ggsave(p, file = "results/figs/HRD/CALGB/hallmarks_dotplot.png", width = 7, height = 15)

######################
## Other signatures ##
######################
## Fix for each cohort
sigs_cols <- colnames(merged_data)[107:155]

cor_results <- data.frame()
pval_results <- data.frame()

for (sig in sigs_cols) {
  for (signature in signatures) {
    cor_test <- cor.test(merged_data[[signature]], merged_data[[sig]], method = "pearson")  # Change to "spearman" if needed
    cor_results <- rbind(cor_results, data.frame(Signature = sig, HRD = signature, Correlation = cor_test$estimate))
    pval_results <- rbind(pval_results, data.frame(Signature = sig, HRD = signature, PValue = cor_test$p.value))
  }
}

pval_results <- pval_results %>%
  mutate(Adj_PValue = p.adjust(PValue, method = "fdr"))

df_plot <- cor_results %>%
  left_join(pval_results, by = c("Signature", "HRD")) %>%
  mutate(Size = -log10(Adj_PValue + 1e-10))  # Avoid log(0)

df_plot <- df_plot %>% mutate(HRD = ifelse(HRD == "HRD_signature_Peng", "Peng et al.",
                                            ifelse(HRD == "HRD_signature_Walens", "Walens et al.",
                                              ifelse(HRD == "HRD_Beinse_scaled", "Beinse et al.",
                                                ifelse(HRD == "HRD_Zhuang_scaled", "Zhuang et al.", "BRCA1ness")))))
lvls.sigs <- c("BRCA1ness", "Peng et al.", "Walens et al.", "Beinse et al.", "Zhuang et al.")
df_plot <- df_plot %>% mutate(Signature = gsub("_PMID.*", "", Signature))

p <- ggplot(df_plot, aes(x = factor(HRD, lvls.sigs), y = factor(Signature, levels = unique(df_plot$Signature)), size = Size, color = Correlation)) +
  geom_point() +
  scale_color_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  scale_size(range = c(2, 8)) +  
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.box = "vertical") +
  guides(
    size = guide_legend(title = "-log10(p-adj)", order = 2),
    color = guide_colorbar(title = "Correlation", order = 3, barwidth = 12, barheight = 1) 
  )      
p
ggsave(p, file = "results/figs/HRD/CALGB/signatures_dotplot.png", width = 7, height = 15)

