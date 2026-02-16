library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

sigs_HRD <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Beinse" =  "HRD_Beinse_scaled", "Zhuang" = "HRD_Zhuang_scaled")

altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")

## 22 samples removed from CALGB due to missing tumor size
# > sum(is.na(calgb$tsizepe))
# [1] 22

altto <- altto %>% mutate(tsize_cm = str_remove(tumour_size, " cm")) %>% mutate(tsize_cm = ifelse(tsize_cm == "=<2", "<=2", ifelse(tsize_cm == ">2 to =<5", ">2 to <=5", tsize_cm)))
calgb_subset <- calgb %>% filter(!is.na(tsizepe)) %>% mutate(tsize_cm = ifelse(tsizepe <= 2, "<=2", ifelse(tsizepe > 5, ">5", ">2 to <=5")))

merged_data <- bind_rows(
  altto %>% mutate(study = "ALTTO") %>% select(study, c("tsize_cm", all_of(unlist(sigs_HRD)))),
  calgb_subset %>% mutate(study = "CALGB") %>% select(study, c("tsize_cm", all_of(unlist(sigs_HRD))))
)

tsize_colors <- c("<=2" = "salmon1", ">2 to <=5" = "salmon3", ">5" = "salmon4")

sig_names <- names(sigs_HRD)
df_long <- merged_data %>%
  pivot_longer(
    cols = all_of(sig_names),
    names_to = "signature",
    values_to = "value"
  ) %>%
  mutate(
    tsize_cm  = factor(tsize_cm, levels = c("<=2", ">2 to <=5", ">5")),
    study     = factor(study, levels = c("ALTTO", "CALGB")),
    signature = factor(signature, levels = sig_names)
  )

p <- ggplot(df_long, aes(x = tsize_cm, y = value, color = tsize_cm)) +
  geom_boxplot(aes(fill = tsize_cm), linewidth = 1.25, alpha = 0.4, outlier.shape = NA) +
  geom_jitter(position = position_jitter(width = 0.2), size = 2.5, alpha = 0.9) +
  scale_color_manual(values = tsize_colors, drop = TRUE) +
  scale_fill_manual(values = tsize_colors, drop = TRUE) +
  facet_grid(
    rows = vars(study),
    cols = vars(signature),
    scales = "free_y",
    switch = "y"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_blank(),
    legend.position = "none",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    strip.text.x = element_text(face = "bold")
  ) +
#   stat_compare_means(label = "p.format", label.y = 1.1)
stat_compare_means(
  aes(
    label = sprintf("p = %.3f", after_stat(p)),
    fontface = after_stat(ifelse(p < 0.05, "bold", "plain"))
  )
)

ggsave(p, file = "results/figs/HRD/SuppFig2_tumor_size.pdf", width = 14, height = 7)
