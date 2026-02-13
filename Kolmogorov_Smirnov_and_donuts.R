library(dplyr)
library(ggplot2)
library(stringr)
library(rlang)
library(ggpubr)

sigs_cols <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Zhuang" = "HRD_Zhuang_scaled", "Beinse" =  "HRD_Beinse_scaled")

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

for (i in seq_len(nrow(res_ks))) { # loop through KS test results df rows
  if(res_ks[i,"p_value"] >= 0.05) next

  print(paste0("Plotting ", res_ks[i,"signature"], " for comparison ", res_ks[i,"comparison"], " with p-value ", round(res_ks[i,"p_value"], 3)))
  df <- merged_data %>% filter(study %in% c(res_ks[i,"study1"], res_ks[i,"study2"]))
  sig_col <- res_ks[i,"signature"]

  p <- ggplot(df, aes(x = df[[sig_col]], fill = study, color = study)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    labs(
        x = sig_col,
        y = "Density",
        color = "Study",
        fill = "Study",
        subtitle = paste0("Kolmogorov-Smirnov test, p = ", round(res_ks[i,"p_value"], 4))
    ) +
    scale_fill_manual(values = study_colors) +
    scale_color_manual(values = study_colors) +
    theme_classic(base_size = 14) 

  ggsave(p, file = paste0("results/figs/HRD/KS_test/", sig_col, "_", res_ks[i,"comparison"], ".png"), width = 8, height = 6)
}

################
### Outliers ###
################

### Get outlier patients
df_quantiles <- merged_data %>% group_by(study) %>%
                    summarise(across(all_of(names(sigs_cols)), ~ quantile(.x, 0.9, na.rm = TRUE)))

mtx_quantiles <- as.matrix(df_quantiles)
rownames(mtx_quantiles) <- mtx_quantiles[, "study"]
mtx_quantiles <- mtx_quantiles[, colnames(mtx_quantiles) != "study"]

outliers_list <- list()
for(sig in names(sigs_cols)) {
    print(sig)
    sig_outliers <- c()
    for(study_id in unique(merged_data$study)) {
        threshold <- mtx_quantiles[study_id, sig]
        outlier_samples <- merged_data %>% filter(study == study_id, !!sym(sig) > threshold) %>% rownames()
        sig_outliers <- c(sig_outliers, outlier_samples)
    }
    outliers_list[[sig]] <- sig_outliers
}

## Donuts
er_colors <- readRDS("data/color_palettes/ER_colors.rds")
names(er_colors) <- c("HR-", "HR+")
subtype_colors <- readRDS("data/color_palettes/HER2_subtypes.rds")
nodes_colors <- c("N0" = "goldenrod1", "N+" = "goldenrod4")

make_donut <- function(data, group_var, pal, legend_title = NULL, plot_title = NULL, disable.legend = FALSE, title = NULL) {
  if(disable.legend) {
    legend_pos = "none"
  } else {
    legend_pos = "right"
  }

  group_var <- rlang::enquo(group_var)

  df_sum <- data %>%
    dplyr::count(!!group_var, name = "n") %>%
    mutate(label = paste0("n = ", n))

  df_sum <- data %>%
    dplyr::count(!!group_var, name = "n") %>%
    mutate(label = paste0("n = ", n))

  if(!is.null(plot_title)){
    n_samples <- sum(df_sum$n)  
    subtitle_lab <- paste0("n = ", n_samples)
  } else {
    subtitle_lab <- NULL
  }

  ggplot(df_sum, aes(x = 2, y = n, fill = !!group_var)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      color = "white",
      size = 4
    ) +
    labs(fill = legend_title, title = plot_title, subtitle = subtitle_lab) +
    scale_fill_manual(values = pal, drop = FALSE) +
    theme_void() +
    theme(
      legend.title = element_text(face = "bold"),
      legend.position = legend_pos, 
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), 
      plot.subtitle = element_text(hjust = 0.5, size = 12) 
    )
}

### Plot donuts for all signatures (90th quantile samples only)
for(sig in names(sigs_cols)) {
    print(sig)

    ### Prep data frames with outliers
    neoaltto_outliers <- neoaltto %>% filter(MATERIAL.ID %in% outliers_list[[sig]]) %>% select(!!sigs_cols[[sig]], erpgrstatus, HER2_subtype, nstage) 
    colnames(neoaltto_outliers) <- c(sig,"HR_status","HER2_subtype", "Nodal_status")
    neoaltto_outliers <- neoaltto_outliers %>% mutate(HR_status = ifelse(HR_status == "POSITIVE", "HR+", "HR-"), Nodal_status = ifelse(Nodal_status == "N0", "N0", "N+"))

    altto_outliers <- altto %>% filter(ID_brightcore %in% outliers_list[[sig]]) %>% select(!!sigs_cols[[sig]], hr, HER2_subtype, nodal_status)
    altto_outliers <- altto_outliers %>% mutate(hr = ifelse(hr == "Positive", "HR+", "HR-"), nodal_status = ifelse(nodal_status == "Node Negative", "N0", "N+"))
    colnames(altto_outliers) <- c(sig,"HR_status","HER2_subtype", "Nodal_status")

    # TO DO: add missing HER2 subtypes to the saved calgb object
    calgb_outliers <- calgb %>% filter(fq %in% outliers_list[[sig]]) %>% select(!!sigs_cols[[sig]], HR_reviewed, HER2_subtype, Clinical_N_Stage)
    calgb_outliers <- calgb_outliers %>% mutate(HR_reviewed = ifelse(HR_reviewed == "pos", "HR+", "HR-"), Clinical_N_Stage = ifelse(Clinical_N_Stage == "N0", "N0", "N+"))
    colnames(calgb_outliers) <- c(sig,"HR_status","HER2_subtype", "Nodal_status")

    ### Plotting
    p_neoaltto_er   <- make_donut(neoaltto_outliers, HR_status, pal = er_colors, legend_title = "Hormone receptor status")
    legend_er <- get_legend(p_neoaltto_er) %>% as_ggplot()
    p_neoaltto_nodes   <- make_donut(neoaltto_outliers, Nodal_status, pal = nodes_colors, legend_title = "Nodal status")
    legend_nodes <- get_legend(p_neoaltto_nodes) %>% as_ggplot()

    p_neoaltto_er   <- make_donut(neoaltto_outliers, HR_status, pal = er_colors, disable.legend = TRUE, plot_title = "NeoALTTO")
    p_neoaltto_sub  <- make_donut(neoaltto_outliers, HER2_subtype, pal = subtype_colors, disable.legend = TRUE)
    p_neoaltto_nodes   <- make_donut(neoaltto_outliers, Nodal_status, pal = nodes_colors, disable.legend = TRUE)
    p_altto_er  <- make_donut(altto_outliers, HR_status, pal = er_colors, disable.legend = TRUE, plot_title = "ALTTO")

    p_altto_sub <- make_donut(altto_outliers, HER2_subtype, pal = subtype_colors, legend_title = "HER2 subtype")
    legend_sub <- get_legend(p_altto_sub)

    p_altto_sub <- make_donut(altto_outliers, HER2_subtype, pal = subtype_colors, disable.legend = TRUE)
    p_altto_nodes   <- make_donut(altto_outliers, Nodal_status, pal = nodes_colors, disable.legend = TRUE)

    p_calgb_er <- make_donut(calgb_outliers, HR_status, pal = er_colors, disable.legend = TRUE, plot_title = "CALGB")
    p_calgb_sub  <- make_donut(calgb_outliers, HER2_subtype, pal = subtype_colors, disable.legend = TRUE)
    p_calgb_nodes   <- make_donut(calgb_outliers, Nodal_status, pal = nodes_colors, disable.legend = TRUE)

    p_final_donuts <- cowplot::plot_grid(p_neoaltto_er, p_altto_er, p_calgb_er,legend_er, 
                                    p_neoaltto_sub, p_altto_sub, p_calgb_sub, legend_sub, 
                                    p_neoaltto_nodes, p_altto_nodes, p_calgb_nodes, legend_nodes,
                    nrow = 3, ncol = 4, rel_widths = c(0.35, 0.35, 0.35, 0.3))

    ggsave(p_final_donuts, file = paste0("results/figs/HRD/donut_charts_outliers_", sig, ".png"), width = 12, height = 9)
}
