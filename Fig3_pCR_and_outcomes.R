library(dplyr)
library(ggplot2)
library(ggpubr)

sigs_HRD <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Beinse" =  "HRD_Beinse_scaled", "Zhuang" = "HRD_Zhuang_scaled")

neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")

source("scripts/HRD/forest_plots_pCR_and_survival.R")

###########
### pCR ###
###########

##  NeoALTTO
neoaltto <- neoaltto %>% mutate(pCR_numeric = ifelse(pCR == "No", 0, 1))
neoaltto <- neoaltto %>% mutate(HR_pos = ifelse(erpgrstatus == "NEGATIVE", 0, 1))

p1 <- make_forest_plot(
  data = neoaltto,
  signatures = sigs_HRD,
  observed_variable = "pCR_numeric",
  covariates = c("ArmCD", "HR_pos"), ## adjustment of HR status and treatment arm (T, L, T+L)
  signatures_labs = rev(names(sigs_HRD)),
  binary = TRUE, # for pCR
  title = "NeoALTTO"
)

## CALGB
calgb <- calgb %>% mutate(HR_pos = ifelse(HR_reviewed == "neg", 0, 1))

# 1. pCR
p2 <- make_forest_plot(
  data = calgb,
  signatures = sigs_HRD,
  observed_variable = "pCR",
  covariates = c("Tx_arm", "HR_pos"),
  signatures_labs = rev(names(sigs_HRD)),
  binary = TRUE, # for pCR
  title = "CALGB"
)

p_pcr <- ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave(p_pcr, file = "results/figs/HRD/Fig3a_pCR_forest_plots.pdf", width = 13, height = 7)

### Supplementary figure: split into HR groups
## NeoALTTO

tmp_data <- neoaltto %>% filter(HR_pos == 0)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs_HRD,
  observed_variable = "pCR_numeric",
  covariates = c("GGI_PMID_16478745", "IgG_HER2DX_PMID_34990895"),
  signatures_labs = rev(names(sigs_HRD)),
  binary = TRUE, # for pCR
  title = "Hormone negative (n = 117)"
)

tmp_data <- neoaltto %>% filter(HR_pos == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs_HRD,
  observed_variable = "pCR_numeric",
  covariates = c("GGI_PMID_16478745", "IgG_HER2DX_PMID_34990895"),
  signatures_labs = rev(names(sigs_HRD)),
  binary = TRUE, # for pCR
  title = "Hormone positive (n = 137)"
)

p <- ggarrange(p1,p2)
ggsave(p, file = "results/figs/HRD/SuppFig3a_pCR_in_HR_groups_NeoALTTO.pdf", width = 13, height = 7)

## CALGB
tmp_data <- calgb %>% filter(HR_pos == 0)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs_HRD,
  observed_variable = "pCR",
  covariates = c("GGI_PMID_16478745", "IgG_HER2DX_PMID_34990895"),
  signatures_labs = rev(names(sigs_HRD)),
  binary = TRUE, # for pCR
  title = "Hormone negative (n = 110)"
)

### Supp. HR groups
tmp_data <- calgb %>% filter(HR_pos == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs_HRD,
  observed_variable = "pCR",
  covariates = c("GGI_PMID_16478745", "IgG_HER2DX_PMID_34990895"),
  signatures_labs = rev(names(sigs_HRD)),
  binary = TRUE, # for pCR
  title = "Hormone positive (n = 154)"
)

p <- ggarrange(p1,p2)
ggsave(p, file = "results/figs/HRD/SuppFig3b_pCR_in_HR_groups_CALGB.pdf", width = 13, height = 7)

### Supplementary figure: treatment arms comparison
## 1. NeoALTTO
treatment_arms <- c("TRASTUZUMAB ALONE", "LAPATINIB ALONE", "LAPATINIB IN COMBINATION WITH TRASTUZUMAB")

p1 <- make_multi_group_forest_plot(
  data           = neoaltto,
  signatures     = sigs_HRD,
  group_var      = "ArmDesc",
  group_levels   = treatment_arms,
  observed_variable = "pCR_numeric",
  covariates     = c("erpgrstatus"),
  signatures_labs= rev(names(sigs_HRD)),
  pos_lab        = "Favors pCR",
  neg_lab        = "Favors non-pCR",
  binary         = TRUE,
  title         = "NeoALTTO"
)

# ggsave(p, file = "results/figs/HRD/SuppFig4a_pCR_forest_plot_treatment_arms.png", width = 8, height = 10)

## CALGB
calgb <- calgb %>% mutate(ArmDesc = ifelse(ArmCD == "A", "LAPATINIB ALONE", ifelse(ArmCD == "B", "TRASTUZUMAB ALONE", "LAPATINIB IN COMBINATION WITH TRASTUZUMAB")))

p2 <- make_multi_group_forest_plot(
  data           = calgb,
  signatures     = sigs_HRD,
  group_var      = "ArmDesc",
  group_levels   = treatment_arms,
  observed_variable = "pCR",
  covariates     = c("HR_pos"),
  signatures_labs= rev(names(sigs_HRD)),
  pos_lab        = "Favors pCR",
  neg_lab        = "Favors non-pCR",
  binary         = TRUE,
  title         = "CALGB"
)

# ggsave(p, file = "results/figs/HRD/CALGB/pCR_forest_plot_treatment_arms_HR.png", width = 8, height = 10)
p <- ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom")
ggsave(p, file = "results/figs/HRD/SuppFig4_pCR_treatment_arms.pdf", width = 13, height = 7)
###############
### EFS/DFS ###
###############

## NeoALTTO
neoaltto <- neoaltto %>% mutate(nstage_numeric = ifelse(nstage %in% c("NX", "N0"), 0, 1))
## !!!! Questionable units in tumor size column !!!!
# merged_data <- merged_data %>% mutate(tsize_cm = ifelse(tsize > 10, tsize/10, tsize))

p1 <- make_forest_plot(
  data = neoaltto,
  signatures = sigs_HRD,
  surv_obj = neoaltto$efs,
  covariates = c("Age", "nstage_numeric", "pCR_numeric"),
  signatures_labs = list(rev(names(sigs_HRD))), # must be list
  binary = FALSE,
  title = "NeoALTTO (EFS)"
)

## CALGB
p2 <- make_forest_plot(
  data = calgb,
  signatures = sigs_HRD,
  time = "efsYears",
  event = "efs",
  covariates = c("Age", "tsizepe", "N_STAGE_GROUPS", "pCR"),
  signatures_labs = list(rev(names(sigs_HRD))), 
  binary = FALSE, 
  title = "CALGB (EFS)"
)

p_efs <- ggarrange(p1, p2, ncol = 2, nrow = 1)
ggsave(p_efs, file = "results/figs/HRD/Fig3b_EFS_forest_plots.pdf", width = 13, height = 7)

## ALTTO
altto$tumour_size <- factor(altto$tumour_size, levels = c("=<2 cm", ">2 to =<5 cm", ">5 cm"))
altto$nstatus <- ifelse(altto$nodal_status == "Node Negative", 0, 1)

# 1. DFS
p3 <- make_forest_plot(
  data = altto,
  signatures = sigs_HRD,
  time = "dfs_years",
  event = "dfs_event",
  covariates = c("age", "tumour_size", "nstatus"),
  signatures_labs = list(rev(names(sigs_HRD))), 
  binary = FALSE,
  title = "ALTTO (DFS)"
)
ggsave(p3, file = "results/figs/HRD/Fig3c_DFS_ALTTO_forest_plot.pdf", width = 7, height = 7)
