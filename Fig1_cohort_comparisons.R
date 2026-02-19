library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggtext)
library(tidyr)
library(xlsx)

sigs_cols <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Beinse" =  "HRD_Beinse_scaled", "Zhuang" = "HRD_Zhuang_scaled")

neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")
tnbc <- readRDS("data/HRD_signatures/TNBC_meta.rds")
head(tnbc)

medians <- data.frame(signature = character(), tnbc_median = numeric())

for(sig in names(sigs_cols)){
    median_value <- median(tnbc[[sigs_cols[[sig]]]], na.rm = TRUE)
    medians <- rbind(medians, data.frame(signature = sig, tnbc_median = median_value))
}

merged_data <- bind_rows(
  altto %>% mutate(study = "ALTTO") %>% select(study, all_of(unlist(sigs_cols))),
  neoaltto %>% mutate(study = "NeoALTTO") %>% select(study, all_of(unlist(sigs_cols))),
  calgb %>% mutate(study = "CALGB") %>% select(study, all_of(unlist(sigs_cols))),
  tnbc %>% mutate(study = "TNBC") %>% select(study, all_of(unlist(sigs_cols)))
)


plots <- list()

for (sig in names(sigs_cols)) {
  sig_median <- medians %>%
    filter(signature == sig) %>%
    pull(tnbc_median)
  
  plot_data <- merged_data %>%
    mutate(
      value = .data[[sig]],
      above_tnbc_median = value > sig_median
    )

  p_val <- kruskal.test(value ~ study, data = plot_data)$p.value

  subtitle_text <- if (p_val < 0.001) {
  "<b>p < 0.001</b>"
  } else if (p_val >= 0.05) {
    paste0("p = ", round(p_val, 3))
  } else {
    paste0("<b>p = ", round(p_val, 3), "</b>")
  }

  main_color <- "steelblue"

  p <- ggplot(plot_data, aes(x = factor(study, levels = c("NeoALTTO", "ALTTO", "CALGB", "TNBC")), y = value)) +
    geom_boxplot(fill = main_color, alpha = 0.5, outlier.shape = NA) +
    geom_jitter(
        aes(fill = above_tnbc_median),
        shape = 21,
        color = main_color,
        width = 0.25,
        size = 2,
        alpha = 0.9
        ) +
    scale_fill_manual(values = c("FALSE" = "white",
                                 "TRUE" = main_color)) +
    geom_hline(yintercept = sig_median,
               linetype = "dashed", linewidth = 1, color = "red") +
    theme_classic(base_size = 16) +
    labs(title = sig,
         y = "Score",
         subtitle = subtitle_text) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 18), 
          plot.subtitle = element_markdown())
  plots[[sig]] <- p
}

plots <- ggarrange(plotlist = plots, ncol = 3, nrow = 2)
plots

ggsave(plots, file = "results/figs/HRD/Fig1_cohort_comparisons.pdf", width = 16, height = 10)


### Get number of samples above TNBC median for each signature
counts_list <- list()
for (sig in names(sigs_cols)) {
  
  sig_median <- medians %>%
    filter(signature == sig) %>%
    pull(tnbc_median)
  
  tmp <- merged_data %>%
    mutate(value = .data[[sig]]) %>%
    group_by(study) %>%
    summarise(
      signature = sig,
      n_above_tnbc_median = sum(value > sig_median, na.rm = TRUE),
      n_total = n(),
      frac_above_tnbc_median = n_above_tnbc_median / n_total,
      percent_above_tnbc_median = paste0(n_above_tnbc_median, " (", round(frac_above_tnbc_median * 100, 2), "%)"),
      .groups = "drop"
    ) 
  
  counts_list[[sig]] <- tmp
}

counts_df <- bind_rows(counts_list) %>% filter(study != "TNBC")

df <- counts_df %>% select(study, signature, percent_above_tnbc_median) %>%
  pivot_wider(names_from = "signature", values_from = "percent_above_tnbc_median") %>% as.data.frame()
rownames(df) <- df$study
df$study <- NULL

write.csv(df, "results/figs/HRD/Fig1_sample_above_TNBC_median.csv", row.names = TRUE)
write.xlsx(df, "results/figs/HRD/Fig1_sample_above_TNBC_median.xlsx", row.names = TRUE)
