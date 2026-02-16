library(dplyr)
library(ggplot2)
library(stringr)
library(rlang)
library(ggpubr)

sigs_cols <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Beinse" =  "HRD_Beinse_scaled",  "Zhuang" = "HRD_Zhuang_scaled")

neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")
head(calgb)


################################
### Kolmogorov-Smirnov tests ###
################################
res_ks <- data.frame(signature = character(), comparison = character(), p_value = numeric())

for(sig in names(sigs_cols)) {
  sig_col <- sigs_cols[[sig]]
  
  # Perform KS test between ALTTO and neoALTTO
  ks_altto_neoaltto <- ks.test(altto[[sig_col]], neoaltto[[sig_col]])
  res_ks <- rbind(res_ks, data.frame(signature = sig, comparison = "ALTTO vs neoALTTO", p_value = ks_altto_neoaltto$p.value))
  
  # Perform KS test between ALTTO and CALGB
  ks_altto_calgb <- ks.test(altto[[sig_col]], calgb[[sig_col]])
  res_ks <- rbind(res_ks, data.frame(signature = sig, comparison = "ALTTO vs CALGB", p_value = ks_altto_calgb$p.value))

    # Perform KS test between NeoALTTO and CALGB
  ks_neoaltto_calgb <- ks.test(neoaltto[[sig_col]], calgb[[sig_col]])
  res_ks <- rbind(res_ks, data.frame(signature = sig, comparison = "neoALTTO vs CALGB", p_value = ks_neoaltto_calgb$p.value))
}

res_ks %>% filter(p_value < 0.05)
#   signature        comparison      p_value
# 1 BRCA1ness neoALTTO vs CALGB 3.207871e-02
# 2      Peng ALTTO vs neoALTTO 1.904194e-02
# 3    Zhuang ALTTO vs neoALTTO 3.035461e-12
# 4    Zhuang    ALTTO vs CALGB 3.772077e-03
# 5    Zhuang neoALTTO vs CALGB 2.042810e-14
# 6    Beinse    ALTTO vs CALGB 8.405795e-03
# 7    Beinse neoALTTO vs CALGB 2.715143e-03

res_ks <- res_ks %>% mutate(study1 = str_split(comparison, " vs ", simplify = TRUE)[,1],
                              study2 = str_split(comparison, " vs ", simplify = TRUE)[,2])

merged_data <- bind_rows(
  altto %>% mutate(study = "ALTTO") %>% select(study, all_of(unlist(sigs_cols))),
  neoaltto %>% mutate(study = "neoALTTO") %>% select(study, all_of(unlist(sigs_cols))),
  calgb %>% mutate(study = "CALGB") %>% select(study, all_of(unlist(sigs_cols)))
)

study_colors <- c("ALTTO" = "deepskyblue3", "neoALTTO" = "chocolate2", "CALGB" = "green4")

dist_plots <- list()
for (i in seq_len(nrow(res_ks))) { # loop through KS test results df rows
  if(res_ks[i,"p_value"] >= 0.05) next

  print(paste0("Plotting ", res_ks[i,"signature"], " for comparison ", res_ks[i,"comparison"], " with p-value ", round(res_ks[i,"p_value"], 3)))
  df <- merged_data %>% filter(study %in% c(res_ks[i,"study1"], res_ks[i,"study2"]))
  sig_col <- res_ks[i,"signature"]

  p <- ggplot(df, aes(x = !!sym(sig_col), fill = study, color = study)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    labs(
        title = sig_col,
        subtitle = paste0("p = ", round(res_ks[i,"p_value"], 3))
    ) +
    scale_fill_manual(values = study_colors[unique(df$study)], drop = FALSE) +
    scale_color_manual(values = study_colors[unique(df$study)], drop = FALSE) +
    theme_classic(base_size = 14) +
    theme(axis.title = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16), 
          legend.position = "bottom",
          legend.title = element_blank())
  dist_plots[[paste0(sig_col, "_", res_ks[i,"comparison"])]] <- p
  # ggsave(p, file = paste0("results/figs/HRD/KS_test/", sig_col, "_", res_ks[i,"comparison"], ".pdf"), width = 8, height = 6)
}

p <- ggarrange(plotlist = dist_plots, ncol = 4, nrow = 2)
ggsave(p, file = "results/figs/HRD/SuppFig5_KS_test.pdf", width = 12, height = 8)
