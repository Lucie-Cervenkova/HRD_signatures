#################
### Hallmarks ###
#################
library(tidyr)
library(stats)
library(dplyr)
library(ggplot2)
library(ggpubr)

sigs_HRD <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Beinse" =  "HRD_Beinse_scaled", "Zhuang" = "HRD_Zhuang_scaled")

neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")
tnbc <- readRDS("data/HRD_signatures/TNBC_meta.rds")

hallmark_cols <- grep("^HALLMARK", names(neoaltto), value = TRUE)  # Detect hallmark columns
other_sig_cols <- c("ESR1_gene", "ERBB2_gene",grep("*PMID*", names(neoaltto), value = TRUE))

merged_data <- bind_rows(
  altto %>% mutate(study = "ALTTO") %>% select(study, all_of(c(hallmark_cols, other_sig_cols, unlist(sigs_HRD)))),
  neoaltto %>% mutate(study = "NeoALTTO") %>% select(study, all_of(c(hallmark_cols, other_sig_cols, unlist(sigs_HRD)))),
  calgb %>% mutate(study = "CALGB") %>% select(study, all_of(c(hallmark_cols, other_sig_cols, unlist(sigs_HRD)))),
)

#################
### Hallmarks ###
#################

cor_results <- data.frame()
pval_results <- data.frame()

for (hallmark in hallmark_cols) {
  for (signature in names(sigs_HRD)) {
    cor_test <- cor.test(merged_data[[signature]], merged_data[[hallmark]], method = "pearson")  # Change to "spearman" if needed
    cor_results <- rbind(cor_results, data.frame(Signature = hallmark, HRD = signature, Correlation = cor_test$estimate))
    pval_results <- rbind(pval_results, data.frame(Signature = hallmark, HRD = signature, PValue = cor_test$p.value))
  }
}

pval_results <- pval_results %>%
  mutate(Adj_PValue = p.adjust(PValue, method = "fdr"))

df_plot <- cor_results %>%
  left_join(pval_results, by = c("Signature", "HRD")) %>%
  mutate(Size = -log10(Adj_PValue + 1e-10))  # Avoid log(0)

df_plot <- df_plot %>% mutate(Signature = gsub("^HALLMARK_", "", Signature))  # Remove "HALLMARK_" prefix

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

df_plot <- df_plot %>% filter(Signature %in% ordered_pathways)

p1 <- ggplot(df_plot, aes(x = factor(HRD, levels = c("BRCA1ness", "Walens", "Peng", "Beinse", "Zhuang")), y = factor(Signature, levels = ordered_pathways), size = Size, color = Correlation)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(1,8)) + 
  labs(size = "-log10(FDR)", color = "Pearson correlation") +
  theme_light() +
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right",
  ) 

ggsave(p1, file = "results/figs/HRD/Fig2a_hallmarks_dotplot.pdf", width = 7, height = 12)

######################
## Other signatures ##
######################
cor_results <- data.frame()
pval_results <- data.frame()

for (sig in other_sig_cols) {
  for (signature in names(sigs_HRD)) {
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

df_plot <- df_plot %>% mutate(Signature = gsub("_PMID.*", "", Signature))

p2 <- ggplot(df_plot, aes(x = factor(HRD, c("BRCA1ness", "Walens", "Peng", "Beinse", "Zhuang")), y = factor(Signature, levels = unique(df_plot$Signature)), size = Size, color = Correlation)) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(1, 8)) + 
  labs(size = "-log10(FDR)", color = "Pearson correlation") + 
  guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(order = 2)
  ) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "right")   

ggsave(p2, file = "results/figs/HRD/Fig2b_signatures_dotplot.pdf", width = 7, height = 12)

p_both <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")
ggsave(p_both, file = "results/figs/HRD/Fig2_gene_signatures_dotplots.pdf", width = 13, height = 12)