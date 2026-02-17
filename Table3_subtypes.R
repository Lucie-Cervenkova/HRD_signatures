library(dplyr)
library(ggplot2)
library(ggpubr)
library(xlsx)

sigs_HRD <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Beinse" =  "HRD_Beinse_scaled", "Zhuang" = "HRD_Zhuang_scaled")

neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")

compare_subtypes <- function(df, signatures, group_var, study = NULL) {
  results <- data.frame()
  groups <- unique(df[[group_var]])

for (grp in groups) {
    # Binary comparison: group A vs. the rest
    df$binary_group <- ifelse(df[[group_var]] == grp, 1, 0)

    for (sig in signatures) {
      w_test <- wilcox.test(df[[sig]] ~ df$binary_group)
      model <- glm(binary_group ~ df[[sig]], data = df, family = "binomial")
      OR <- exp(coef(model)[2]) 
      
      results <- rbind(results, data.frame(
        subtype = grp,
        signature = sig,
        p_value = w_test$p.value,
        odds_ratio = OR,
        study = study
      ))
    }
  }

  # FDR Correction PER group
  results <- results %>%
    group_by(subtype) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
    ungroup()
  return(results[,c("subtype", "signature", "odds_ratio", "p_value", "p_adj", "study")])
}

res_neoaltto <- compare_subtypes(neoaltto, sigs_HRD, "HER2_subtype", study = "neoALTTO")
res_altto <- compare_subtypes(altto, sigs_HRD, "HER2_subtype", study = "ALTTO")
res_calgb <- compare_subtypes(calgb, sigs_HRD, "HER2_subtype", study = "CALGB")

res_all <- rbind(res_neoaltto, res_altto, res_calgb)

write.xlsx(res_all, file = "results/figs/HRD/Table3_subtypes.xlsx", sheetName = "HER2 subtypes", append = FALSE)
write.csv(res_all, file = "results/figs/HRD/Table3_subtypes.csv", row.names = FALSE)

res_neoaltto <- compare_subtypes(neoaltto, sigs_HRD, "PAM50", study = "neoALTTO")
res_altto <- compare_subtypes(altto, sigs_HRD, "AIMS", study = "ALTTO")
res_calgb <- compare_subtypes(calgb, sigs_HRD, "PAM50", study = "CALGB")

res_all <- rbind(res_neoaltto, res_altto, res_calgb)
write.xlsx(res_all, file = "results/figs/HRD/Table3_subtypes.xlsx", sheetName = "PAM50", append = TRUE)
write.csv(res_all, file = "results/figs/HRD/Table3_PAM50.csv", row.names = FALSE)