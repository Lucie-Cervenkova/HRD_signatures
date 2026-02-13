library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)

sigs_HRD <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Beinse" =  "HRD_Beinse_scaled", "Zhuang" = "HRD_Zhuang_scaled")

neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")

merged_data <- bind_rows(
  altto %>% mutate(study = "ALTTO") %>% select(study, all_of(unlist(sigs_HRD))),
  neoaltto %>% mutate(study = "NeoALTTO") %>% select(study, all_of(unlist(sigs_HRD))),
  calgb %>% mutate(study = "CALGB") %>% select(study, all_of(unlist(sigs_HRD))))

df_long <- merged_data %>%
  pivot_longer(
    cols = all_of(names(sigs_HRD)),
    names_to = "signature",
    values_to = "value"
  ) %>% mutate(study = factor(study, levels = c("NeoALTTO", "ALTTO", "CALGB")), signature = factor(signature, levels = names(sigs_HRD)))


### All signatures
p <- ggplot(df_long, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40,
                 fill = "red",
                 color = "black",
                 alpha = 0.6) +
  geom_density(color = "red",
               linewidth = 0.8,
               alpha = 0.5) +
  facet_grid(
    rows = vars(study),
    cols = vars(signature),
    scales = "free_x",
    switch = "y",
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold")
  )

ggsave(p, file = "results/figs/HRD/SuppFig1_distributions.pdf", width = 12, height = 8)

### Exclude Zhuang
df_long_nozhuang <- df_long %>% filter(signature != "Zhuang")

p <- ggplot(df_long_nozhuang, aes(x = value)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40,
                 fill = "red",
                 color = "black",
                 alpha = 0.6) +
  geom_density(color = "red",
               linewidth = 0.8,
               alpha = 0.5) +
  facet_grid(
    rows = vars(study),
    cols = vars(signature),
    scales = "free_x",
    switch = "y",
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.title = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold")
  )

ggsave(p, file = "results/figs/HRD/SuppFig1_distributions_no_Zhuang.pdf", width = 12, height = 8)
