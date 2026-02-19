### TO DO: make some functions to easily load data etc.

library(dplyr)
library(GSVA)
library(DESeq2)
library(stringr)
library(ggplot2)
library(readxl)
library(rstatix)
library(ggpubr)
library(oddsratio)
library(pracma)
######################
### Get signatures ###
######################


genesets <- list()

# BRCA1ness
brca1ness_sig <- read_excel("/mnt/bctl/gene_signatures/BRCA1ness.xlsx")
brca1ness <- c(brca1ness_sig$`77_genes`)
genesets$BRCA1ness <- brca1ness

# HRD Walens et al.
HRD_Walens <- readRDS("/mnt/bctl/gene_signatures/HRD_Walens.rds")
genesets$HRD_signature_Walens <- HRD_Walens

# Peng et al. HRD microarray signature
# Creating an R vector with the genes from Supplementary Table 1
HRD_Peng <- readRDS("/mnt/bctl/gene_signatures/HRD_Peng.rds")
genesets$HRD_signature_Peng <- HRD_Peng

genesets %>% lapply(., length) %>% unlist

###################
### Load inputs ###
###################

#################
## 1. neoALTTO ##
#################

# ## Load data matrix
# dir <- "/mnt/bctl/mattia_data/neoALTTO/"
# txi <- readRDS(paste0(dir, "txi_neoaltto.rds"))
# mtx <- txi$counts # raw counts
# tpm <- txi$abundance # TPM

# ## Load metadata
# clinical <- load("/mnt/bctl/davidV/neoALTTO/Clinical.RData") %>% get()
# subtype_data <- read_excel("/mnt/bctl/mattia_data/neoALTTO/source_data_NeoALTTO_GEX.xlsx")

# clinical <- clinical %>% filter(MATERIAL.ID %in% subtype_data$MATERIAL.ID)
# rownames(clinical) <- clinical$MATERIAL.ID

# meta <- left_join(clinical, subtype_data, by = "MATERIAL.ID")
# rownames(meta) <- meta$MATERIAL.ID

# ## VST normalization
# meta <- meta[colnames(txi$counts), ]
# dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ ArmCD)
# keep <- rowMeans(counts(dds)) >= 10 # only keep gened where rowMeans >= 10
# dds_filtered <- dds[keep, ]
# vsd <- vst(dds_filtered, blind = TRUE)  # blind=TRUE means no prior knowledge of conditions
# vst_counts <- assay(vsd)
# mtx <- vst_counts

# # Prep input for HRD200 = Pan et al. classifier
# tpm_df <- t(tpm) %>% as.data.frame
# tpm_df$id <- rownames(tpm_df)
# tpm_df <- tpm_df[,c("id",colnames(tpm_df)[-length(colnames(tpm_df))])]

# write.table(tpm_df, file="/mnt/bctl/lcer0007/HRD200/data/neoALTTO.csv", sep=",", row.names = FALSE, quote=FALSE)

##############
## 2. ALTTO ##
##############

## Load data matrix
# dir <- "/mnt/bctl/mattia_data/ALTTO/"
# txi <- readRDS(paste0(dir, "txi_altto.RDS"))
# mtx <- txi$counts # raw counts
# tpm <- txi$abundance # TPM
# vst <- readRDS(paste0(dir, "rna_norm_vst.RDS")) ## filtered out genes with counts rowMeans < 10 & high quality patients
# ## only keep the genes and samples in filtered VST counts
# tpm <- tpm[rownames(vst), colnames(vst)]
# mtx <- mtx[rownames(mtx), colnames(mtx)]

# ## Load metadata
# clinical <- read_excel(paste0(dir, "ALTTO_FULL_CLINICAL_CASE_CONTROL_COHORT.xlsx")) ## clinical data
# subtype_data <- read_excel(paste0(dir,"source_data_ALTTO_GEX.xlsx")) ## subtypes + signatures, n = 386
# subtype_data <- as.data.frame(subtype_data)
# rownames(subtype_data) <- subtype_data$ID_brightcore

# clinical <- clinical %>% filter(ID_brightcore %in% colnames(vst))
# rownames(clinical) <- clinical$ID_brightcore

# meta <- left_join(clinical, subtype_data, by = "ID_brightcore")
# rownames(meta) <- meta$ID_brightcore

# mtx <- vst

# ## Prep input for HRD200 = Pan et al. classifier
# tpm_df <- t(tpm) %>% as.data.frame
# tpm_df$id <- rownames(tpm_df)
# tpm_df <- tpm_df[,c("id",colnames(tpm_df)[-length(colnames(tpm_df))])]

# write.table(tpm_df, file="/mnt/bctl/lcer0007/HRD200/data/ALTTO.csv", sep=",", row.names = FALSE, quote=FALSE)

##############
## 3. CALGB ##
##############

# Load data matrix
# dir <- "/mnt/bctl/mattia_data/CALGB/"
# txi <- readRDS(paste0(dir, "txi_calgb.RDS"))
# mtx <- txi$counts # raw counts
# tpm <- txi$abundance # TPM

# ## Load metadata
# # Clinical
# clinical <- read_excel(paste0(dir, "CALGB_clinical_bcr_tcr.xlsx")) ## Note: use updated EFS/OS/pCR
# # PAM50
# pam50_meta <- read_excel(paste0(dir, "CALGB_PAM50_paper.xlsx"))
# # Hallmarks
# hallmarks_df <- read_excel(paste0(dir, "CALGB_GSVA_H_paper.xlsx"))

# key <- read_excel(paste0(dir, "paper_KEY_CALGB.xlsx")) # convert IDs to fq in clinical data
# pam50_meta <- left_join(pam50_meta, key, by = "PT_ID_calgb")
# hallmarks_df <- left_join(hallmarks_df, key, by = "PT_ID_calgb")
# pam50_meta$PT_ID_calgb <- NULL
# hallmarks_df$PT_ID_calgb <- NULL

# # Join all into metadata
# meta <- left_join(clinical, pam50_meta, by = "fq")
# meta <- left_join(meta, hallmarks_df, by = "fq")
# meta <- as.data.frame(meta)
# rownames(meta) <- meta$fq

# ## VST normalization
# meta <- meta[colnames(txi$counts), ]
# dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ Tx_arm)
# keep <- rowMeans(counts(dds)) >= 10 # only keep gened where rowMeans >= 10
# dds_filtered <- dds[keep, ]
# vsd <- vst(dds_filtered, blind = TRUE)  # blind=TRUE means no prior knowledge of conditions
# vst_counts <- assay(vsd)
# mtx <- vst_counts

# ## Compute signatures
# sigs <- readRDS("data/signatures_REF.RDS")
# sigs_subset <- sigs[!str_detect(names(sigs), "HALLMARK")]
# source("/mnt/bctl/Gene_signatures/calcSig.R")
# res <- calcSig(d = mtx, sig = sigs_subset) # matrix must be normalized first (mtx = VST counts)
# # Add to metadata
# meta <- cbind(meta, res)

# ## Prep input for HRD200 = Pan et al. classifier
# tpm_df <- t(tpm) %>% as.data.frame
# tpm_df$id <- rownames(tpm_df)
# tpm_df <- tpm_df[,c("id",colnames(tpm_df)[-length(colnames(tpm_df))])]

# write.table(tpm_df, file="/mnt/bctl/lcer0007/HRD200/data/CALGB.csv", sep=",", row.names = FALSE, quote=FALSE)

# tpm_df[1:5, 1:5]

#############
## 4. TCGA ##
#############

## Load data matrix
load("/mnt/bctl/TCGA/bctl_backup/rnaseqv2.RData")
mtx <- res

## Load metadata
load("/mnt/bctl/TCGA/bctl_backup/Clinical.RData")
clinical <- res

## Subset TNBCs
names(clinical) <- make.unique(names(clinical))
meta <- clinical %>% filter(type == "TNBC")

## Rename colnames in matrix
short_ids <- sub("^(([^-]+-){2}[^-]+).*", "\\1", colnames(mtx))
# rownames(meta)[!(rownames(meta) %in% short_ids)]
colnames(mtx) <- short_ids

mtx <- mtx[,intersect(rownames(meta), short_ids)]
meta <- meta[colnames(mtx), ]

##############
## Run GSVA ##
##############

gsva_results <- gsva(
  expr = mtx,      # Your normalized expression matrix
  gset.idx.list = genesets, # Your list of gene sets
  method = "gsva",        # Choose the method: "gsva", "ssgsea", "zscore", or "plage"
  kcdf = "Gaussian",      # Use "Gaussian" for continuous data like RNA-seq VST counts
  min.sz = 20,            # Minimum gene set size
  max.sz = 10000,           # Maximum gene set size
  verbose = TRUE
)

gsva_df <- as.data.frame(t(gsva_results)) 

# rescale
gsva_df <- as.data.frame(apply(gsva_df, 2, function(col) {
  2 * (col - min(col)) / (max(col) - min(col)) - 1
}))

head(gsva_df)

# ## neoALTTO
# gsva_df$MATERIAL.ID <- colnames(gsva_results)
# gsva_df$MATERIAL.ID <- as.integer(gsva_df$MATERIAL.ID)
# merged_data <- gsva_df %>%
#   left_join(meta, by = "MATERIAL.ID")
# rownames(merged_data) <- merged_data$MATERIAL.ID

# ## ALTTO
# gsva_df$ID_brightcore <- colnames(gsva_results)
# merged_data <- gsva_df %>%
#   left_join(meta, by = "ID_brightcore")
# rownames(merged_data) <- merged_data$ID_brightcore

# ## CALGB
# gsva_df$fq <- rownames(gsva_df)
# merged_data <- gsva_df %>%
#   left_join(meta, by = "fq")
# rownames(merged_data) <- merged_data$fq

## TCGA
merged_data <- cbind(meta, gsva_df)

## HRD: Zhuang et al. gene-pair signature 
HRD_Zhuang <- data.frame(A = c('CDT1', 'DEPDC1', 'DEPDC1', 'CEP55', 'CDC45', 'CENPA', 'E2F1', 'PLK4', 
                               'CDC45','CENPE', 'CKAP2L', 'BUB1', 'DEPDC1', 'KPNA2', 'CENPA', 'PLK4', 
                               'CDT1', 'GINS2', 'ORC6', 'ANLN', 'BIRC5', 'EXO1', 'FOXM1', 'CDC25C'),
                         B = c('POU4F1', 'POU4F1', 'MS4A1', 'TMEM176A', 'CREB3L3', 'UGT1A1', 'S100A1', 'UGT1A5',
                              'POU2AF1', 'GBA3', 'EFCAB1', 'HNF4G', 'SPATA22', 'SWAP70','TMPRSS11BNL', 'ZBP1', 
                                'KY', 'STK4', 'ZC3H7A', 'HNF4G', 'SLC43A1', 'SCYL3', 'ZNF124', 'EIF3K'))

computeZhuangSig <- function(tpm, sig_df) {
    HRD_genes <- unique(c(sig_df$A, sig_df$B))

    missing_genes <- setdiff(HRD_genes, rownames(tpm))
    print("Missing genes from Zhuang et al. signature:")
    print(missing_genes)

    HRD_genes <- intersect(HRD_genes, rownames(tpm))

    if (length(HRD_genes) == 0) {
        stop("Error: No valid genes found in TPM matrix.")
    }

    tpm_filtered <- tpm[HRD_genes, , drop = FALSE]  # Safe subsetting

    total_scores <- rep(NA, ncol(tpm_filtered))

    for (col_num in 1:ncol(tpm_filtered)) {
        score <- 0

        for (row_num in 1:nrow(sig_df)) {
            A_gene <- sig_df[row_num, "A"]
            B_gene <- sig_df[row_num, "B"]

            if (!(A_gene %in% rownames(tpm_filtered)) | !(B_gene %in% rownames(tpm_filtered))) {
                next
            }

            if (tpm_filtered[A_gene, col_num] > tpm_filtered[B_gene, col_num]) {
                score <- score + 1
            }
        }
        total_scores[col_num] <- score / 24
    }

    return(total_scores)
}

HRD_Zhuang_scores <- computeZhuangSig(mtx, HRD_Zhuang)

sum(HRD_Zhuang_scores >= 0.65)
sum(HRD_Zhuang_scores < 0.65)

merged_data$HRD_Zhuang <- HRD_Zhuang_scores
merged_data <- merged_data %>% mutate(HRD_Zhuang_groups = ifelse(HRD_Zhuang >= 0.65, "high HRD", "low HRD"))

### HRD Beinse et al
library(glmnet)

dir_Beinse <- "/mnt/bctl/gene_signatures/HRD_Beinse/"

cv.ridge.nuclear <- readRDS(paste0(dir_Beinse, "RNAseq_nuclear_model"))
cv.ridge.cyto <- readRDS(paste0(dir_Beinse, "RNAseq_cytoplasmic_model"))
HRD_model <- readRDS(paste0(dir_Beinse, "RNAseq_model"))

# log_tpm <- log2(tpm + 1)
log_tpm <- log2(mtx + 1) ## TNBC 

## Get ENSEMBL ID instead of HGNC gene names
library(biomaRt)

## Sometimes not working because of the DB connection :(
# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

mart <- useEnsembl(
  biomart = "ensembl",
  dataset = "hsapiens_gene_ensembl",
  mirror  = "asia"     # options: "www", "useast", "uswest", "asia"
)

converted_genes <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), 
                         filters="hgnc_symbol", 
                         values=rownames(log_tpm), 
                         mart=mart)

log_tpm <- log_tpm[converted_genes$hgnc_symbol,]
rownames(log_tpm) <- converted_genes$ensembl_gene_id

## See how many genes are present
# rownames(cv.ridge.nuclear$glmnet.fit$beta) %>% length
# intersect(rownames(cv.ridge.nuclear$glmnet.fit$beta), rownames(log_tpm)) %>% length

nuclear_genes <- intersect(rownames(cv.ridge.nuclear$glmnet.fit$beta), rownames(log_tpm))
missing_genes <- setdiff(rownames(cv.ridge.nuclear$glmnet.fit$beta), rownames(log_tpm))

zero_matrix <- matrix(0, nrow = length(missing_genes), ncol = ncol(log_tpm))
rownames(zero_matrix) <- missing_genes
colnames(zero_matrix) <- colnames(log_tpm)

# Append to the original log_tpm matrix
nuclear_log_tpm <- rbind(log_tpm, zero_matrix)

data_for_nuclear_score <- nuclear_log_tpm %>%
  as.data.frame()%>%
  t()%>%
  as.data.frame()

nuclear_score <- as.numeric(predict(object = cv.ridge.nuclear,
                                    as.matrix(data_for_nuclear_score[,rownames(cv.ridge.nuclear$glmnet.fit$beta)]),
                                    s = cv.ridge.nuclear$lambda.min,
                                    type="response"))

## Missing cytoplasmic genes
missing_genes <- setdiff(rownames(cv.ridge.cyto$glmnet.fit$beta), rownames(log_tpm))

zero_matrix <- matrix(0, nrow = length(missing_genes), ncol = ncol(log_tpm))
rownames(zero_matrix) <- missing_genes
colnames(zero_matrix) <- colnames(log_tpm)

# Append to the original log_tpm matrix
cyto_log_tpm <- rbind(log_tpm, zero_matrix)

data_for_cytoplasmic_score <- cyto_log_tpm %>%
  as.data.frame()%>%
  # dplyr::select(all_of(c(list_of_sample))) %>%
  t()%>%
  as.data.frame()


cytoplasmic_score <- as.numeric(predict(object = cv.ridge.cyto,
                                                          as.matrix(data_for_cytoplasmic_score[,rownames(cv.ridge.cyto$glmnet.fit$beta)]),
                                                          s = cv.ridge.cyto$lambda.min,
                                                          type="response"))

results <- data.frame(nuclear_score = as.numeric(nuclear_score),
                      cytoplasmic_score = as.numeric(cytoplasmic_score))

score <- predict(HRD_model, results[,c("cytoplasmic_score","nuclear_score")])

## > 42 = high HRD group (in Beinse et al. paper)
merged_data$HRD_Beinse <- score
merged_data <- merged_data %>% mutate(HRD_Beinse_groups = ifelse(HRD_Beinse > 42, "high HRD", "low HRD"))

table(merged_data$HRD_Beinse_groups)
table(merged_data$HRD_Zhuang_groups)

### Scale and save signatures

scaled_sigs_df <- as.data.frame(apply(merged_data[,c("HRD_Zhuang", "HRD_Beinse")], 2, function(col) {
  2 * (col - min(col)) / (max(col) - min(col)) - 1
}))
colnames(scaled_sigs_df) <- c("HRD_Zhuang_scaled", "HRD_Beinse_scaled")

merged_data <- cbind(merged_data, scaled_sigs_df)


saveRDS(merged_data, file = "data/HRD_signatures/TNBC_meta.rds")

#######################
### Load data frame ###
#######################
merged_tnbc <- readRDS("data/HRD_signatures/TNBC_meta.rds")
merged_calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")
merged_altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
merged_neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")

signatures <- c("BRCA1ness", "HRD_signature_Peng", "HRD_signature_Walens", "HRD_Beinse_scaled", "HRD_Zhuang_scaled")

cohorts <- list(
  CALGB = merged_calgb,
  ALTTO = merged_altto,
  neoALTTO = merged_neoaltto
)

# Function to compute overlap coefficient between two numeric vectors
compute_overlap <- function(vec1, vec2) {
  vec1 <- vec1[!is.na(vec1)]
  vec2 <- vec2[!is.na(vec2)]
  
  if (length(vec1) < 5 || length(vec2) < 5 || sd(vec1) == 0 || sd(vec2) == 0) {
    return(NA)
  }
  
  dens1 <- tryCatch(density(vec1), error = function(e) return(NULL))
  dens2 <- tryCatch(density(vec2), error = function(e) return(NULL))
  
  if (is.null(dens1) || is.null(dens2)) return(NA)

  x_common <- seq(min(dens1$x, dens2$x), max(dens1$x, dens2$x), length.out = 1000)
  y1 <- approx(dens1$x, dens1$y, xout = x_common, rule = 2)$y  # rule=2 extrapolates constant
  y2 <- approx(dens2$x, dens2$y, xout = x_common, rule = 2)$y

  y_min <- pmin(y1, y2)

  if (anyNA(y_min)) {
    valid <- !is.na(y_min)
    overlap <- trapz(x_common[valid], y_min[valid])
  } else {
    overlap <- trapz(x_common, y_min)
  }

  return(overlap)
}

# Initialize results matrix
results <- matrix(NA, nrow = length(cohorts), ncol = length(signatures),
                  dimnames = list(names(cohorts), signatures))

# Loop over cohorts and signatures
for (cohort_name in names(cohorts)) {
  cohort_data <- cohorts[[cohort_name]]
  cat("Checking cohort:", cohort_name, "\n")
  
  for (sig in signatures) {
    cat("  Signature:", sig, " ... ")

    if (sig %in% colnames(merged_tnbc) && sig %in% colnames(cohort_data)) {
      cat("exists — computing\n")
      
      vec_tnbc <- merged_tnbc[[sig]]
      vec_other <- cohort_data[[sig]]

      # Filter to Q4 (upper 25%)
      q3_tnbc <- quantile(vec_tnbc, 0.75, na.rm = TRUE)
      q3_other <- quantile(vec_other, 0.75, na.rm = TRUE)

      vec_tnbc_q4 <- vec_tnbc[vec_tnbc >= q3_tnbc]
      vec_other_q4 <- vec_other[vec_other >= q3_other]
      
      results[cohort_name, sig] <- compute_overlap(vec_tnbc_q4, vec_other_q4)
    } else {
      cat("NOT FOUND\n")
    }
  }
}


# Convert to data frame
overlap_df <- as.data.frame(results)
overlap_df

sum(is.na(merged_tnbc[,signatures]))
##########################################
### Signature distributions/histograms ###
##########################################

p1 <- ggplot(merged_data, aes(x = BRCA1ness)) +
  geom_histogram(aes(y = ..density..), fill = "red", bins = 50, alpha = 0.6, color = "black") +
  geom_density(size = 1, alpha = 0.4, color = "red") +
  labs(title = "BRCA1ness",
       x = "",
       y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom") 
  
p2 <- ggplot(merged_data, aes(x = HRD_signature_Peng)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "red", alpha = 0.6, color = "black") +
  geom_density(size = 1, alpha = 0.4, color = "red") +
  labs(title = "Peng et al.",
       x = "",
       y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

p3 <- ggplot(merged_data, aes(x = HRD_signature_Walens)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "red", alpha = 0.6, color = "black") +
  geom_density(size = 1, alpha = 0.4, color = "red") +
  labs(title = "Walens et al.",
       y = "Density",
       x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

p4 <- ggplot(merged_data, aes(x = HRD_Beinse_scaled)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "red", alpha = 0.6, color = "black") +
  geom_density(size = 1, alpha = 0.4, color = "red") +
  labs(title = "Beinse et al.",
       x = "",
       y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

p5 <- ggplot(merged_data, aes(x = HRD_Zhuang_scaled)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "red", alpha = 0.6, color = "black") +
  geom_density(size = 1, alpha = 0.4, color = "red") +
  labs(title = "Zhuang et al.",
       y = "Density",
       x = "") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

p <- ggarrange(p1, p2, p3, p4, p5, ncol = 3, nrow = 2)

## replace with appropriate cohort name
ggsave(p, file = "results/figs/HRD/TCGA_histograms_sigs_distribution.png")

######################### 
### Density in groups ###
#########################
library(plyr)

grouped_density_plot <- function(df, group_var, group_var_lab = NULL, obs_var, obs_var_lab = NULL, custom_colors = NULL){
  
  if(is.null(group_var_lab)){
    group_var_lab = group_var
  }

  if(is.null(obs_var_lab)){
    obs_var_lab = obs_var
  }
  print(obs_var_lab)
  # mu <- ddply(merged_data, obs_var, summarise, grp.mean=mean(!!sym(obs_var)))
  
  p <- ggplot(merged_data, aes(x = !!sym(obs_var), fill = !!sym(group_var))) +
      # geom_histogram(aes(y = ..density.., fill = !!sym(group_var)), bins = 50, alpha = 0.4, color = "black") +
       geom_density(alpha= 0.5) +
      #  geom_vline(data=mu, aes(xintercept=grp.mean, color=hr), linetype="dashed") +
      #  scale_fill_manual(values = subtype_colors) +
      #  scale_color_manual(values = subtype_colors) +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(title = obs_var_lab, y = "", x = "",
            fill = group_var_lab, color = group_var_lab)
  return(p)
}

grouped_histogram <- function(df, group_var, group_var_lab = NULL, obs_var, obs_var_lab = NULL, custom_colors = NULL){
  
  if(is.null(group_var_lab)){
    group_var_lab = group_var
  }

  if(is.null(obs_var_lab)){
    obs_var_lab = obs_var
  }
  print(obs_var_lab)
  # mu <- ddply(merged_data, obs_var, summarise, grp.mean=mean(!!sym(obs_var)))
  
  p <- ggplot(merged_data, aes(x = !!sym(obs_var), fill = !!sym(group_var))) +
      geom_histogram(aes(y = ..density.., fill = !!sym(group_var)), bins = 50, alpha = 0.4, color = "black") +
      #  geom_density(alpha= 0.5) +
      #  geom_vline(data=mu, aes(xintercept=grp.mean, color=hr), linetype="dashed") +
      #  scale_fill_manual(values = subtype_colors) +
      #  scale_color_manual(values = subtype_colors) +
      scale_y_continuous(expand = expansion(mult = c(0, 0))) +
      theme_classic() +
      theme(legend.position = "bottom",
            plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(title = obs_var_lab, y = "", x = "",
            fill = group_var_lab, color = group_var_lab)
  return(p)
}

split_histogram <- function(df, group_var, group_var_lab = NULL, obs_var, obs_var_lab = NULL, custom_colors = NULL){
  
  if(is.null(group_var_lab)){
    group_var_lab = group_var
  }

  if(is.null(obs_var_lab)){
    obs_var_lab = obs_var
  }

  print(obs_var_lab)

  # mu <- ddply(df, group_var, summarise, grp.mean=mean(!!sym(obs_var)))

  mu <- df %>%
    dplyr::group_by(!!sym(group_var)) %>%
    dplyr::summarise(grp.mean = mean(!!sym(obs_var)), na.rm = TRUE)

  p <- ggplot(df, aes(x = !!sym(obs_var), fill = !!sym(group_var))) +
    geom_histogram(aes(y = ..density..), 
                   bins = 50, 
                   alpha = 0.4, 
                   color = "black") +
    geom_density(alpha= 0.5) +
    geom_vline(data = mu, aes(xintercept = grp.mean), color = "black", linetype = "dashed") +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.text = element_text(size = 12)) +
    labs(title = obs_var_lab, y = "", x = "",
         fill = group_var_lab, color = group_var_lab) +
    facet_wrap(vars(.data[[group_var]]), ncol = 1, scales = "free_y")

  return(p)
}

ER_lab <- "Estrogen receptor status"
p1 <- split_histogram(merged_data, group_var = "hr", group_var_lab = ER_lab, obs_var = "BRCA1ness")
p2 <- split_histogram(merged_data, group_var = "hr", group_var_lab = ER_lab, obs_var = "HRD_Beinse_scaled", obs_var_lab = "Beinse et al.")
p3 <- split_histogram(merged_data, group_var = "hr", group_var_lab = ER_lab, obs_var = "HRD_Zhuang_scaled", obs_var_lab = "Zhuang et al.")
p4 <- split_histogram(merged_data, group_var = "hr", group_var_lab = ER_lab, obs_var = "HRD_signature_Peng", obs_var_lab = "Peng et al.")
p5 <- split_histogram(merged_data, group_var = "hr", group_var_lab = ER_lab, obs_var = "HRD_signature_Walens", obs_var_lab = "Walens et al.")


merged_data <- merged_data %>% mutate(age_group = ifelse(age < 40, "<40", ">=40"))
age_lab <- "Age"

table(merged_data$age_group)

p1 <- split_histogram(merged_data, group_var = "age_group", group_var_lab = age_lab, obs_var = "BRCA1ness")
p2 <- split_histogram(merged_data, group_var = "age_group", group_var_lab = age_lab, obs_var = "HRD_Beinse_scaled", obs_var_lab = "Beinse et al.")
p3 <- split_histogram(merged_data, group_var = "age_group", group_var_lab = age_lab, obs_var = "HRD_Zhuang_scaled", obs_var_lab = "Zhuang et al.")
p4 <- split_histogram(merged_data, group_var = "age_group", group_var_lab = age_lab, obs_var = "HRD_signature_Peng", obs_var_lab = "Peng et al.")
p5 <- split_histogram(merged_data, group_var = "age_group", group_var_lab = age_lab, obs_var = "HRD_signature_Walens", obs_var_lab = "Walens et al.")

###################################
### Clinical variables boxplots ###
###################################
library(RColorBrewer)
library(ggpubr)
library(rlang)

clinical_variable_boxplot <- function(merged_data, clinical_variable, variable_lab = NULL, signature, signature_lab = NULL, 
                                      cohort_name, save_plot = TRUE, custom_colors = NULL, custom_levels = NULL, y_lab = NULL, return_plot = FALSE){
  ##### WARNING: This way we do NOT perform multiple testing correction across all HRD signatures (potentially fine, but beware if p value = 0.01 etc.)
  ##### use forest plots instead

  if (length(unique(merged_data[[clinical_variable]])) < 2) {
    message(paste("Skipping", cohort_name, "- Not enough unique values for statistical testing"))
    return(NULL)
  }

  if(length(unique(merged_data[[clinical_variable]])) > 2){
    # Kruskal-Wallis test
    kruskal_res <- kruskal.test(merged_data[[signature]], merged_data[[clinical_variable]])
    print(paste0("Kruskal-Wallis p-val = ", kruskal_res$p.value))

    ### Skip plotting if Kruskal-Wallis is insignificant
    # if (kruskal_res$p.value >= 0.05) {
    #   message(paste("Skipping - Kruskal-Wallis p-value:", kruskal_res$p.value))
    #   return(NULL)
    # }

    subtitle_lab = paste0("Kruskal-Wallis p-value: ", signif(kruskal_res$p.value, 3))
  } 

  ### new
  if(is.null(custom_levels)){
    merged_data[[clinical_variable]] <- factor(merged_data[[clinical_variable]], levels = unique(merged_data[[clinical_variable]]))
  }
  else {
    merged_data[[clinical_variable]] <- factor(merged_data[[clinical_variable]], levels = custom_levels)
  }
  


  pairwise_wilcox <- tryCatch({
      pairwise.wilcox.test(merged_data[[signature]], merged_data[[clinical_variable]], p.adj = "BH") %>% 
        broom::tidy() %>%
        mutate(sig_label = paste0("p.adj = ", signif(p.value, 3)))
    }, error = function(e) {
      message(paste("Skipping - Wilcoxon test failed:", e$message))
      return(NULL)
    })

    if (is.null(pairwise_wilcox)) return(NULL)  # Skip plotting if Wilcoxon test failed
    
    print(pairwise_wilcox)

    if(length(unique(merged_data[[clinical_variable]])) == 2) {
      subtitle_lab = paste0("Wilcoxon p-value: ", signif(pairwise_wilcox$p.value, 3))
    }
    comparisons_list <- combn(as.character(unique(merged_data[[clinical_variable]])), 2, simplify = FALSE)

    # comparisons_list <- combn(unique(merged_data[[clinical_variable]]), 2, simplify = FALSE)

  if(is.null(variable_lab)){
    variable_lab = clinical_variable
  }

  if(is.null(signature_lab)){
    signature_lab = signature
  }
  
  if(is.null(y_lab)){
    y_lab = "Signature score"
  }

  p <- ggplot(merged_data, aes(x = !!sym(clinical_variable), y = !!sym(signature), fill = !!sym(clinical_variable))) +
    geom_boxplot() +
    geom_jitter(shape = 16, position = position_jitter(0.2)) +
    theme_classic() +
    labs(title = signature_lab,
        subtitle = subtitle_lab,
        x = variable_lab,
        y = y_lab) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin(t = 5, b = 5)),
          plot.subtitle = element_text(hjust = 0.5, size = 12),
          legend.position = "none") +
    # stat_compare_means(comparisons = comparisons_list) 
    stat_compare_means(comparisons = comparisons_list, method = "wilcox")
 

  if (!is.null(custom_colors)) {
    p <- p + scale_fill_manual(values = custom_colors)
  }

  if(save_plot){
    ggsave(p, file = paste0("results/figs/HRD/", cohort_name, "/", signature,"_boxplot_", clinical_variable,".png"), width = 7, height = 8)
  }

  if(return_plot){
    return(p)
  }

} 

## neoALTTO
# Hormone receptors
unique(merged_data$erpgrstatus)
hr_levels <- c("NEGATIVE", "POSITIVE")

clinical_variable_boxplot(merged_data, clinical_variable = "erstatus", variable_lab = "Estrogen receptor status", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "erstatus", variable_lab = "Estrogen receptor status", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "erstatus", variable_lab = "Estrogen receptor status", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "erstatus", variable_lab = "Estrogen receptor status", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "erstatus", variable_lab = "Estrogen receptor status", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)

clinical_variable_boxplot(merged_data, clinical_variable = "pgrstatus", variable_lab = "Progesteron receptor status", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "pgrstatus", variable_lab = "Progesteron receptor status", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "pgrstatus", variable_lab = "Progesteron receptor status", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "pgrstatus", variable_lab = "Progesteron receptor status", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "pgrstatus", variable_lab = "Progesteron receptor status", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)

clinical_variable_boxplot(merged_data, clinical_variable = "erpgrstatus", variable_lab = "Hormone receptor status", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "erpgrstatus", variable_lab = "Hormone receptor status", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "erpgrstatus", variable_lab = "Hormone receptor status", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "erpgrstatus", variable_lab = "Hormone receptor status", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "erpgrstatus", variable_lab = "Hormone receptor status", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "neoALTTO", custom_levels = hr_levels)

# Grade
histograde_data <- merged_data %>% filter(histograde %in% c("G1", "G2", "G3"))
grade_levels <- c("G1", "G2", "G3")

# ns for Walens
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = grade_levels)
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_levels = grade_levels)
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "neoALTTO", custom_levels = grade_levels)
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "neoALTTO", custom_levels = grade_levels)
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "neoALTTO", custom_levels = grade_levels)

# Nodal stage
nodal_data <- merged_data %>% filter(nstage != "NX")

unique(nodal_data$nstage)
nodes_levels <- c("N0","N1", "N2A", "N3B","N3C")

# nothing significant
clinical_variable_boxplot(nodal_data, clinical_variable = "nstage", variable_lab = "Nodal stage", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = nodes_levels)
clinical_variable_boxplot(nodal_data, clinical_variable = "nstage", variable_lab = "Nodal stage", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_levels = nodes_levels)
clinical_variable_boxplot(nodal_data, clinical_variable = "nstage", variable_lab = "Nodal stage", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "neoALTTO", custom_levels = nodes_levels)
clinical_variable_boxplot(nodal_data, clinical_variable = "nstage", variable_lab = "Nodal stage", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "neoALTTO", custom_levels = nodes_levels)
clinical_variable_boxplot(nodal_data, clinical_variable = "nstage", variable_lab = "Nodal stage", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "neoALTTO", custom_levels = nodes_levels)

# age
merged_data <- merged_data %>% mutate(age_group = ifelse(Age < 50, "<50", ">=50"))
age_levels <- c("<50", ">=50")

# ns
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = age_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_levels = age_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "neoALTTO", custom_levels = age_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "neoALTTO", custom_levels = age_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "neoALTTO", custom_levels = age_levels)


# menopause
menop_data <- merged_data %>% filter(menop %in% c("PREMENOPAUSAL", "POSTMENOPAUSAL"))
menop_levels <- c("PREMENOPAUSAL", "POSTMENOPAUSAL")

clinical_variable_boxplot(menop_data, clinical_variable = "menop", variable_lab = "", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = menop_levels)
clinical_variable_boxplot(menop_data, clinical_variable = "menop", variable_lab = "", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_levels = menop_levels)
clinical_variable_boxplot(menop_data, clinical_variable = "menop", variable_lab = "", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "neoALTTO", custom_levels = menop_levels)
clinical_variable_boxplot(menop_data, clinical_variable = "menop", variable_lab = "", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "neoALTTO", custom_levels = menop_levels)
clinical_variable_boxplot(menop_data, clinical_variable = "menop", variable_lab = "", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "neoALTTO", custom_levels = menop_levels)


# Treatment arm
unique(merged_data$ArmDesc)
## Key: A = Lapatinib alone, B = Trastuzumab alone, C = Combination
arms_levels <- c("A", "B", "C")

# ns
clinical_variable_boxplot(merged_data, clinical_variable = "ArmCD", variable_lab = "Treatment arm", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "ArmCD", variable_lab = "Treatment arm", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "ArmCD", variable_lab = "Treatment arm", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "neoALTTO", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "ArmCD", variable_lab = "Treatment arm", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "neoALTTO", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "ArmCD", variable_lab = "Treatment arm", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "neoALTTO", custom_levels = arms_levels)


# ### Continuous clinical variables
# cor.test(merged_data$tsize, merged_data$BRCA1ness, method = "spearman")
# cor.test(merged_data$tsize, merged_data$HRD_signature_Peng, method = "spearman")
# cor.test(merged_data$tsize, merged_data$HRD_signature_Walens, method = "spearman")

# cor.test(merged_data$height, merged_data$BRCA1ness, method = "spearman")
# cor.test(merged_data$height, merged_data$HRD_signature_Peng, method = "spearman")
# cor.test(merged_data$height, merged_data$HRD_signature_Walens, method = "spearman")

# cor.test(merged_data$weight, merged_data$BRCA1ness, method = "spearman")
# cor.test(merged_data$weight, merged_data$HRD_signature_Peng, method = "spearman")
# cor.test(merged_data$weight, merged_data$HRD_signature_Walens, method = "spearman")



#### ALTTO
hr_levels <- c("Negative", "Positive")
# Hormone receptors
clinical_variable_boxplot(merged_data, clinical_variable = "hr", variable_lab = "Hormone receptor status", signature = "BRCA1ness", cohort_name = "ALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "hr", variable_lab = "Hormone receptor status", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "ALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "hr", variable_lab = "Hormone receptor status", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "ALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "hr", variable_lab = "Hormone receptor status", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "ALTTO", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "hr", variable_lab = "Hormone receptor status", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "ALTTO", custom_levels = hr_levels)

# Grade
merged_data <- merged_data %>%
  mutate(histograde = as.factor(hgrade))

histograde_data <- merged_data %>% filter(histograde %in% c("G1: Well differentiated", "G2: Moderately differentiated", "G3: Poorly differentiated/Undifferentiated"))
histograde_data <- histograde_data %>% mutate(histograde = ifelse(histograde == "G1: Well differentiated", "G1",
                                                                    ifelse(histograde == "G2: Moderately differentiated", "G2", "G3")))
grade_levels <- c("G1", "G2", "G3")

clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "BRCA1ness", cohort_name = "ALTTO", custom_levels = grade_levels)
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "ALTTO", custom_levels = grade_levels)
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "ALTTO", custom_levels = grade_levels)
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "ALTTO", custom_levels = grade_levels)
clinical_variable_boxplot(histograde_data, clinical_variable = "histograde", variable_lab = "Histological grade", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "ALTTO", custom_levels = grade_levels)

# Nodal stage
# merged_data <- merged_data %>%
#   mutate(nodal_status = as.factor(nodal_status))
unique(merged_data$nodal_status)
nodes_levels <- c("Node Negative", "1-3 Positive Nodes", ">=4 Positive Nodes")

# nothing significant
clinical_variable_boxplot(merged_data, clinical_variable = "nodal_status", variable_lab = "Nodal stage", signature = "BRCA1ness", cohort_name = "ALTTO", custom_levels = nodes_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "nodal_status", variable_lab = "Nodal stage", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "ALTTO", custom_levels = nodes_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "nodal_status", variable_lab = "Nodal stage", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "ALTTO", custom_levels = nodes_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "nodal_status", variable_lab = "Nodal stage", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "ALTTO", custom_levels = nodes_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "nodal_status", variable_lab = "Nodal stage", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "ALTTO", custom_levels = nodes_levels)

# age
merged_data <- merged_data %>% mutate(age_group = ifelse(age < 40, "<40", ">=40"))
merged_data$age_group <- factor(merged_data$age_group, levels = c("<40", ">=40"))
table(merged_data$age_group)

levels(merged_data$age_group)
# ns
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "BRCA1ness", cohort_name = "ALTTO", custom_levels = levels(merged_data$age_group))
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "ALTTO", custom_levels = levels(merged_data$age_group))
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "ALTTO", custom_levels = levels(merged_data$age_group))
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "ALTTO", custom_levels = levels(merged_data$age_group))
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "ALTTO", custom_levels = levels(merged_data$age_group))

# Tumour size
unique(merged_data$tumour_size)
tsize_levels <- c("=<2 cm", ">2 to =<5 cm", ">5 cm")

# only Peng significant - basically the same as grade relationship but opposite?
clinical_variable_boxplot(merged_data, clinical_variable = "tumour_size", variable_lab = "Tumor size", signature = "BRCA1ness", cohort_name = "ALTTO", custom_levels = tsize_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "tumour_size", variable_lab = "Tumor size", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "ALTTO", custom_levels = tsize_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "tumour_size", variable_lab = "Tumor size", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "ALTTO", custom_levels = tsize_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "tumour_size", variable_lab = "Tumor size", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "ALTTO", custom_levels = tsize_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "tumour_size", variable_lab = "Tumor size", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "ALTTO", custom_levels = tsize_levels)

# Treatment arm
table(merged_data$chemo)
arms_levels <- c("Concurrent", "Sequential")

# not significant
clinical_variable_boxplot(merged_data, clinical_variable = "chemo", variable_lab = "Chemotherapy timing", signature = "BRCA1ness", cohort_name = "ALTTO", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "chemo", variable_lab = "Chemotherapy timing", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "ALTTO", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "chemo", variable_lab = "Chemotherapy timing", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "ALTTO", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "chemo", variable_lab = "Chemotherapy timing", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "ALTTO", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "chemo", variable_lab = "Chemotherapy timing", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "ALTTO", custom_levels = arms_levels)

#### CALGB
# Hormone receptors
unique(merged_data$PR_reviewed)
hr_levels <- c("neg", "pos")

# as expected, except Walens is switched ??
clinical_variable_boxplot(merged_data, clinical_variable = "ER_reviewed", variable_lab = "Estrogen receptor status", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "ER_reviewed", variable_lab = "Estrogen receptor status", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "ER_reviewed", variable_lab = "Estrogen receptor status", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "ER_reviewed", variable_lab = "Estrogen receptor status", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "ER_reviewed", variable_lab = "Estrogen receptor status", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = hr_levels)

clinical_variable_boxplot(merged_data, clinical_variable = "PR_reviewed", variable_lab = "Progesteron receptor status", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PR_reviewed", variable_lab = "Progesteron receptor status", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PR_reviewed", variable_lab = "Progesteron receptor status", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PR_reviewed", variable_lab = "Progesteron receptor status", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PR_reviewed", variable_lab = "Progesteron receptor status", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = hr_levels)

clinical_variable_boxplot(merged_data, clinical_variable = "HR_reviewed", variable_lab = "Hormone receptor status", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HR_reviewed", variable_lab = "Hormone receptor status", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HR_reviewed", variable_lab = "Hormone receptor status", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HR_reviewed", variable_lab = "Hormone receptor status", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = hr_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HR_reviewed", variable_lab = "Hormone receptor status", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = hr_levels)

# Nodal stage
nodal_data <- merged_data %>% filter(!is.na(Clinical_N_Stage) & Clinical_N_Stage != "NA")
unique(nodal_data$Clinical_N_Stage)

nodes_levels <- c("N0","N1","N2","N3")
# nothing significant
clinical_variable_boxplot(nodal_data, clinical_variable = "Clinical_N_Stage", variable_lab = "Nodal stage", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = nodes_levels)
clinical_variable_boxplot(nodal_data, clinical_variable = "Clinical_N_Stage", variable_lab = "Nodal stage", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = nodes_levels)
clinical_variable_boxplot(nodal_data, clinical_variable = "Clinical_N_Stage", variable_lab = "Nodal stage", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = nodes_levels)
clinical_variable_boxplot(nodal_data, clinical_variable = "Clinical_N_Stage", variable_lab = "Nodal stage", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = nodes_levels)
clinical_variable_boxplot(nodal_data, clinical_variable = "Clinical_N_Stage", variable_lab = "Nodal stage", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = nodes_levels)

# age
merged_data <- merged_data %>% mutate(age_group = ifelse(Age < 50, "<50", ">=50"))
table(merged_data$age_group)
# ns
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = age_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = age_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = age_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = age_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "age_group", variable_lab = "Age", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = age_levels)

# Tumour size ## TO DO!!! levels
tstage_data <- merged_data %>% filter(!is.na(T_stage))
# tstage_data <- tstage_data %>% mutate(tstage = ifelse(T_stage == ">=T3", "T3/4", T_stage))
# head(tstage_data)
unique(tstage_data$T_stage)
tstage_levels <- c("T1", "T2", ">=T3")

# only BRCA1ness
clinical_variable_boxplot(tstage_data, clinical_variable = "T_stage", variable_lab = "Tumor stage", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = tstage_levels)
clinical_variable_boxplot(tstage_data, clinical_variable = "T_stage", variable_lab = "Tumor stage", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = tstage_levels)
clinical_variable_boxplot(tstage_data, clinical_variable = "T_stage", variable_lab = "Tumor stage", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = tstage_levels)
clinical_variable_boxplot(tstage_data, clinical_variable = "T_stage", variable_lab = "Tumor stage", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = tstage_levels)
clinical_variable_boxplot(tstage_data, clinical_variable = "T_stage", variable_lab = "Tumor stage", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = tstage_levels)

# Treatment arm
unique(merged_data$Arm_combi_single)
arms_levels <- c("single", "combi")

# not significant
clinical_variable_boxplot(merged_data, clinical_variable = "Arm_combi_single", variable_lab = "Treatment arm", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "Arm_combi_single", variable_lab = "Treatment arm", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "Arm_combi_single", variable_lab = "Treatment arm", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "Arm_combi_single", variable_lab = "Treatment arm", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = arms_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "Arm_combi_single", variable_lab = "Treatment arm", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = arms_levels)

# PAM50
pam50_levels <- c("Basal", "Her2", "LumB", "LumA", "Normal")

clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = pam50_levels)

## Split by pCR
tmp_data <- merged_data %>% filter(pCR == "Yes") 

clinical_variable_boxplot(tmp_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = pam50_levels)
clinical_variable_boxplot(tmp_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = pam50_levels)
clinical_variable_boxplot(tmp_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = pam50_levels)
clinical_variable_boxplot(tmp_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = pam50_levels)
clinical_variable_boxplot(tmp_data, clinical_variable = "PAM50", variable_lab = "PAM50", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = pam50_levels)

## Tumor size
merged_data <- merged_data %>% mutate(tsize_groups = ifelse(is.na(tsizepe), NA ,ifelse(tsizepe <= 2, "<=2", ifelse(tsizepe <= 5, "2-5", ">5"))))
tsize_levels <- c("<=2", "2-5", ">5")
tsize_data <- merged_data %>% filter(!is.na(tsize_groups))

clinical_variable_boxplot(tsize_data, clinical_variable = "tsize_groups", variable_lab = "Tumor size (cm)", signature = "BRCA1ness", cohort_name = "CALGB", custom_levels = tsize_levels)
clinical_variable_boxplot(tsize_data, clinical_variable = "tsize_groups", variable_lab = "Tumor size (cm)", signature = "HRD_signature_Walens", signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_levels = tsize_levels)
clinical_variable_boxplot(tsize_data, clinical_variable = "tsize_groups", variable_lab = "Tumor size (cm)", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)", cohort_name = "CALGB", custom_levels = tsize_levels)
clinical_variable_boxplot(tsize_data, clinical_variable = "tsize_groups", variable_lab = "Tumor size (cm)", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)", cohort_name = "CALGB", custom_levels = tsize_levels)
clinical_variable_boxplot(tsize_data, clinical_variable = "tsize_groups", variable_lab = "Tumor size (cm)", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)", cohort_name = "CALGB", custom_levels = tsize_levels)


#####################
#### Forest plots ###
#####################
source("/mnt/bctl/lcer0007/her2/scripts/neoALTTO/my_forest_plots.R")

sigs <- c("BRCA1ness", "HRD_signature_Peng", "HRD_signature_Walens", "HRD_Beinse_scaled", "HRD_Zhuang_scaled")
sigs_labs <- c("BRCA1ness", "Peng et al.", "Walens et al.", "Beinse et al.", "Zhuang et al.")

#### neoALTTO
merged_data <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")

# 1. pCR
merged_data <- merged_data %>% mutate(pCR_numeric = ifelse(pCR == "No", 0, 1))
merged_data <- merged_data %>% mutate(HR_pos = ifelse(erpgrstatus == "NEGATIVE", 0, 1))

p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "pCR_numeric",
  covariates = c("ArmCD", "HR_pos"),
  signatures_labs = sigs_labs,
  binary = TRUE # for pCR
)
p
ggsave(p, file = "results/figs/HRD/neoALTTO/pCR_forest_HR_treatment.png", width = 6, height = 6)


### HR status groups split:
tmp_data <- merged_data %>% filter(HR_pos == 0)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  observed_variable = "pCR_numeric",
  covariates = c("GGI_PMID_16478745", "IgG_HER2DX_PMID_34990895"),
  signatures_labs = sigs_labs,
  binary = TRUE, # for pCR
  title = "Hormone negative (n = 117)"
)

tmp_data <- merged_data %>% filter(HR_pos == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  observed_variable = "pCR_numeric",
  covariates = c("GGI_PMID_16478745", "IgG_HER2DX_PMID_34990895"),
  signatures_labs = sigs_labs,
  binary = TRUE, # for pCR
  title = "Hormone positive (n = 137)"
)

p <- ggarrange(p1,p2)
ggsave(p, file = "results/figs/HRD/neoALTTO/pCR_forest_plot_HR_effect_IgG_GGI.png", width = 12, height = 6)

# 2. EFS
merged_data <- merged_data %>% mutate(nstage_numeric = ifelse(nstage %in% c("NX", "N0"), 0, 1))

## !!!! Questionable units in tumor size column !!!!
# merged_data <- merged_data %>% mutate(tsize_cm = ifelse(tsize > 10, tsize/10, tsize))

p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  surv_obj = merged_data$efs,
  covariates = c("Age", "nstage_numeric", "pCR_numeric"),
  signatures_labs = sigs_labs, # must be list
  binary = FALSE
)
p
ggsave(p, file = "results/figs/HRD/neoALTTO/EFS_forest_plot_age_tsize_nstatus_pCR.png", width = 6, height = 6)

### HR status groups:

tmp_data <- merged_data %>% filter(HR_pos == 0)
dim(tmp_data)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  surv_obj = tmp_data$efs,
  covariates = c("Age", "nstage_numeric", "pCR_numeric"),
  signatures_labs = sigs_labs, # must be list
  binary = FALSE,
  title = "Hormone negative (n = 117)"
)

tmp_data <- merged_data %>% filter(HR_pos == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  surv_obj = tmp_data$efs,
  covariates = c("Age", "nstage_numeric", "pCR_numeric"),
  signatures_labs = sigs_labs,
  binary = FALSE, 
  title = "Hormone positive (n = 170)"
)

p <- ggarrange(p1,p2)
ggsave(p, file = "results/figs/HRD/neoALTTO/EFS_forest_plot_HR_effect_age_nstage_pCR.png", width = 12, height = 6)

### Nodal stage groups:
tmp_data <- merged_data %>% filter(nstage_numeric == 0)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  surv_obj = tmp_data$efs,
  covariates = c("Age", "pCR_numeric"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Node negative (n = 68)"
)

tmp_data <- merged_data %>% filter(nstage_numeric == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  surv_obj = tmp_data$efs,
  covariates = c("Age", "pCR_numeric"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Node positive (n = 186)"
)

p <- ggarrange(p1, p2)
ggsave(p, file = "results/figs/HRD/neoALTTO/EFS_forest_plot_node_effect_age_pCR.png", width = 12, height = 6)


# 3. OS
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  surv_obj = merged_data$os,
  covariates = c("Age", "nstage_numeric", "pCR_numeric"),
  signatures_labs = sigs_labs,
  binary = FALSE
)
p
ggsave(p, file = "results/figs/HRD/neoALTTO/OS_forest_plot_age_tsize_nstatus_pCR.png", width = 6, height = 6)

#### CALGB
merged_data <- readRDS("data/HRD_signatures/CALGB_meta.rds")
merged_data <- merged_data %>% mutate(HR_pos = ifelse(HR_reviewed == "neg", 0, 1))

# 1. pCR
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "pCR",
  covariates = c("Tx_arm", "HR_pos"),
  signatures_labs = sigs_labs,
  binary = TRUE # for pCR
)
p
ggsave(p, file = "results/figs/HRD/CALGB/pCR_forest_plot_HR_treatement.png", width = 6, height = 6)


### HR status groups split:
merged_data <- merged_data %>% mutate(HR_pos = ifelse(HR_reviewed == "neg", 0, 1))

tmp_data <- merged_data %>% filter(HR_pos == 0)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  observed_variable = "pCR",
  covariates = c("GGI_PMID_16478745", "IgG_HER2DX_PMID_34990895"),
  signatures_labs = sigs_labs,
  binary = TRUE, # for pCR
  title = "Hormone negative (n = 110)"
)

tmp_data <- merged_data %>% filter(HR_pos == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  observed_variable = "pCR",
  covariates = c("GGI_PMID_16478745", "IgG_HER2DX_PMID_34990895"),
  signatures_labs = sigs_labs,
  binary = TRUE, # for pCR
  title = "Hormone positive (n = 154)"
)

p <- ggarrange(p1,p2)
ggsave(p, file = "results/figs/HRD/CALGB/pCR_forest_plot_HR_effect_IgG_GGI.png", width = 12, height = 6)


# 2. EFS
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  time = "efsYears",
  event = "efs",
  covariates = c("Age", "tsizepe", "N_STAGE_GROUPS", "pCR"),
  signatures_labs = sigs_labs, 
  binary = FALSE
)
p
ggsave(p, file = "results/figs/HRD/CALGB/EFS_forest_plot_age_tsize_nstatus_pCR.png", width = 6, height = 6)

### HR groups:
tmp_data <- merged_data %>% filter(HR_pos == 0)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  time = "efsYears",
  event = "efs",
  covariates = c("Age", "tsizepe", "N_STAGE_GROUPS", "pCR"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Hormone negative (n = 110)"
)

tmp_data <- merged_data %>% filter(HR_pos == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  time = "efsYears",
  event = "efs",
  covariates = c("Age", "tsizepe", "N_STAGE_GROUPS", "pCR"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Hormone positive (n = 154)"
)

p <- ggarrange(p1, p2)
ggsave(p, file = "results/figs/HRD/CALGB/EFS_forest_plot_HR_effect_age_tsize_nstatus_pCR.png", width = 12, height = 6)

### nodal stage groups:
tmp_data <- merged_data %>% filter(N_STAGE_GROUPS == 0)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  time = "efsYears",
  event = "efs",
  covariates = c("Age", "tsizepe", "N_STAGE_GROUPS", "pCR"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Node negative (n = 115)"
)

tmp_data <- merged_data %>% filter(N_STAGE_GROUPS == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  time = "efsYears",
  event = "efs",
  covariates = c("Age", "tsizepe", "N_STAGE_GROUPS", "pCR"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Node positive (n = 113)"
)

p <- ggarrange(p1, p2)
ggsave(p, file = "results/figs/HRD/CALGB/EFS_forest_plot_node_effect_age_tsize_nstatus_pCR.png", width = 12, height = 6)



# 3. OS
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  time = "osYears_2",
  event = "os",
  covariates = c("Age", "tsizepe", "N_STAGE_GROUPS", "pCR"),
  signatures_labs = sigs_labs,
  binary = FALSE
)
p
ggsave(p, file = "results/figs/HRD/CALGB/OS_forest_plot_age_tsize_nstatus_pCR.png", width = 6, height = 6)

#### ALTTO
merged_data <- readRDS("data/HRD_signatures/ALTTO_meta.rds")

merged_data$tumour_size <- factor(merged_data$tumour_size, levels = c("=<2 cm", ">2 to =<5 cm", ">5 cm"))
merged_data <- merged_data %>% mutate(nstatus = ifelse(nodal_status == "Node Negative", 0, 1))

# 1. DFS
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  time = "dfs_years",
  event = "dfs_event",
  covariates = c("age", "tumour_size", "nstatus"),
  signatures_labs = sigs_labs, 
  binary = FALSE
)
p
ggsave(p, file = "results/figs/HRD/ALTTO/DFS_forest_plot_age_tsize_nstatus.png", width = 6, height = 6)

### HR groups:
tmp_data <- merged_data %>% filter(hr == "Negative")

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  time = "dfs_years",
  event = "dfs_event",
  covariates = c("age", "tumour_size", "nstatus"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Hormone negative (n = 174)"
)

tmp_data <- merged_data %>% filter(hr == "Positive")

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  time = "dfs_years",
  event = "dfs_event",
  covariates = c("age", "tumour_size", "nstatus"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Hormone positive (n = 212)"
)

p <- ggarrange(p1, p2)
ggsave(p, file = "results/figs/HRD/ALTTO/EFS_forest_plot_HR_effect_age_tsize_nstatus.png", width = 12, height = 6)

### nodal stage groups:
tmp_data <- merged_data %>% filter(nstatus == 0)

p1 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  time = "dfs_years",
  event = "dfs_event",
  covariates = c("age", "tumour_size", "nstatus"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Node negative (n = 93)"
)

tmp_data <- merged_data %>% filter(nstatus == 1)

p2 <- make_forest_plot(
  data = tmp_data,
  signatures = sigs,
  time = "dfs_years",
  event = "dfs_event",
  covariates = c("age", "tumour_size", "nstatus"),
  signatures_labs = sigs_labs, 
  binary = FALSE,
    title = "Node positive (n = 293)"
)

p <- ggarrange(p1, p2)
ggsave(p, file = "results/figs/HRD/ALTTO/EFS_forest_plot_node_effect_age_tsize_nstatus.png", width = 12, height = 6)



# 2. OS
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  time = "os_years",
  event = "os_event",
  covariates = c("age", "tumour_size", "nstatus"),
  signatures_labs = sigs_labs,
  binary = FALSE
)
p
ggsave(p, file = "results/figs/HRD/ALTTO/OS_forest_plot_age_tsize_nstatus.png", width = 6, height = 6)



####################################### 
## Clinical variables - forest plots ##
#######################################

#### neoALTTO
merged_data <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")

merged_data <- merged_data %>% mutate(HR_pos = ifelse(erpgrstatus == "NEGATIVE", 0, 1))
merged_data <- merged_data %>% mutate(ER_pos = ifelse(erstatus == "NEGATIVE", 0, 1))
merged_data <- merged_data %>% mutate(PR_pos = ifelse(pgrstatus == "NEGATIVE", 0, 1))

# 1. HR
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "HR_pos",
  neg_lab = "HR-",
  pos_lab = "HR+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/neoALTTO/forest_plots/HR.png", width = 5, height = 7)

# 2. ER
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "ER_pos",
  neg_lab = "ER-",
  pos_lab = "ER+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/neoALTTO/forest_plots/ER.png", width = 5, height = 7)

# 3. PR
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "PR_pos",
  neg_lab = "PR-",
  pos_lab = "PR+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/neoALTTO/forest_plots/PR.png", width = 5, height = 7)

# 4. Nodes
nodes_data <- merged_data %>% filter(nstage != "NX")
nodes_data <- nodes_data %>% mutate(nodes = ifelse(nstage == "N0", 0, 1))

p <- make_forest_plot(
  data = nodes_data,
  signatures = sigs,
  observed_variable = "nodes",
  neg_lab = "N0",
  pos_lab = "N+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/neoALTTO/forest_plots/nodal_status.png", width = 5, height = 7)

# 5. Grade
histo_data <- merged_data %>% filter(!(histograde %in% c("GX","NK")))
histo_data <- histo_data %>% mutate(grade = ifelse(histograde == "G3", 1, 0))

p <- make_forest_plot(
  data = histo_data,
  signatures = sigs,
  observed_variable = "grade",
  neg_lab = "G1-2",
  pos_lab = "G3",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/neoALTTO/forest_plots/grade_binary.png", width = 5, height = 7)

# 6. Age
merged_data <- merged_data %>% mutate(age_group = ifelse(Age >= 40, 1, 0)) # v2: 40, v1: 50
table(merged_data$age_group)

p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "age_group",
  neg_lab = "<40",
  pos_lab = ">=40",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/neoALTTO/forest_plots/age_binary_v2.png", width = 5, height = 7)

# 7. Menopause
menop_data <- merged_data %>% filter(menop %in% c("PREMENOPAUSAL", "POSTMENOPAUSAL"))
menop_data <- menop_data %>% mutate(menop_numeric = ifelse(menop ==  "POSTMENOPAUSAL", 1, 0))

p <- make_forest_plot(
  data = menop_data,
  signatures = sigs,
  observed_variable = "menop_numeric",
  neg_lab = "Premenopausal",
  pos_lab = "Postmenopausal",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/neoALTTO/forest_plots/menop_binary.png", width = 5, height = 7)

#### CALGB
merged_data <- readRDS("data/HRD_signatures/CALGB_meta.rds")

merged_data <- merged_data %>% mutate(HR_pos = ifelse(HR_reviewed == "neg", 0, 1))
merged_data <- merged_data %>% mutate(ER_pos = ifelse(ER_reviewed == "neg", 0, 1))
merged_data <- merged_data %>% mutate(PR_pos = ifelse(PR_reviewed == "neg", 0, 1))

# 1. HR
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "HR_pos",
  neg_lab = "HR-",
  pos_lab = "HR+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/CALGB/forest_plots/HR.png", width = 5, height = 7)

# 2. ER
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "ER_pos",
  neg_lab = "ER-",
  pos_lab = "ER+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/CALGB/forest_plots/ER.png", width = 5, height = 7)

# 3. PR
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "PR_pos",
  neg_lab = "PR-",
  pos_lab = "PR+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/CALGB/forest_plots/PR.png", width = 5, height = 7)

# 4. Nodes
nodes_data <- merged_data %>% filter(!(Clinical_N_Stage %in% c("NA", NA)))
nodes_data <- nodes_data %>% mutate(nodes = ifelse(Clinical_N_Stage == "N0", 0, 1))

p <- make_forest_plot(
  data = nodes_data,
  signatures = sigs,
  observed_variable = "nodes",
  neg_lab = "N0",
  pos_lab = "N+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/CALGB/forest_plots/nodal_status.png", width = 5, height = 7)


# 5. Age
merged_data <- merged_data %>% mutate(age_group = ifelse(Age_on_Study >= 40, 1, 0))
table(merged_data$age_group)

p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "age_group",
  neg_lab = "<40",
  pos_lab = ">=40",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/CALGB/forest_plots/age_binary_v2.png", width = 5, height = 7)

# # 6. Tumor size
# p <- make_forest_plot(
#   data = merged_data,
#   signatures = sigs,
#   observed_variable = "age_group",
#   neg_lab = "<40",
#   pos_lab = ">=40",
#   signatures_labs = sigs_labs,
#   binary = TRUE # logistic regression
# )

### ALTTO
merged_data <- readRDS("data/HRD_signatures/ALTTO_meta.rds")

merged_data <- merged_data %>% mutate(HR_pos = ifelse(hr == "Negative", 0, 1))

# 1. HR
p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "HR_pos",
  neg_lab = "HR-",
  pos_lab = "HR+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/ALTTO/forest_plots/HR.png", width = 5, height = 7)

# 2. Grade
histo_data <- merged_data %>% filter(hgrade != "GX: Differentiation cannot be assessed")
histo_data <- histo_data %>% mutate(grade = ifelse(hgrade == "G3: Poorly differentiated/Undifferentiated", 1, 0))

p <- make_forest_plot(
  data = histo_data,
  signatures = sigs,
  observed_variable = "grade",
  neg_lab = "G1-2",
  pos_lab = "G3",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/ALTTO/forest_plots/grade_binary.png", width = 5, height = 7)


# 4. Nodes
merged_data <- merged_data %>% mutate(nodes = ifelse(nodal_status == "Node Negative", 0, 1))

p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "nodes",
  neg_lab = "N0",
  pos_lab = "N+",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/ALTTO/forest_plots/nodal_status.png", width = 5, height = 7)


# 5. Age
merged_data <- merged_data %>% mutate(age_group = ifelse(age >= 40, 1, 0))

p <- make_forest_plot(
  data = merged_data,
  signatures = sigs,
  observed_variable = "age_group",
  neg_lab = "<40",
  pos_lab = ">=40",
  signatures_labs = sigs_labs,
  binary = TRUE # logistic regression
)

ggsave(p, file = "results/figs/HRD/ALTTO/forest_plots/age_binary_v2.png", width = 5, height = 7)

####################
### Bubble plots ###
####################
source("/mnt/bctl/lcer0007/her2/scripts/neoALTTO/bubble_plot.R")

sigs <- c("BRCA1ness", "HRD_signature_Peng", "HRD_signature_Walens", "HRD_Beinse_scaled", "HRD_Zhuang_scaled")
sigs_labs <- c("BRCA1ness", "Peng et al.", "Walens et al.", "Beinse et al.", "Zhuang et al.")

### neoALTTO
merged_data <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")

subtypes_order <- c("IM", "P/Met", "Mes/S", "LUM", "ERBB2-D")
p <- signature_dotplot(df = merged_data, signatures = sigs, signatures_labs = sigs_labs, group_var = "HER2_subtype", group_order = subtypes_order)
p

ggsave(p, file = "results/figs/HRD/neoALTTO/dotplot_HER2_subtypes.png", width = 6, height = 6)

pam50_levels <- c("Basal", "Her2", "LumB", "LumA", "Normal")
p <- signature_dotplot(df = merged_data, signatures = sigs, signatures_labs = sigs_labs, group_var = "PAM50", group_order = pam50_levels)
ggsave(p, file = "results/figs/HRD/neoALTTO/dotplot_PAM50_subtypes.png", width = 6, height = 6)

### CALGB
merged_data <- readRDS("data/HRD_signatures/CALGB_meta.rds")

subtypes_order <- c("IM", "P/Met", "Mes/S", "LUM", "ERBB2-D")
p <- signature_dotplot(df = merged_data, signatures = sigs, signatures_labs = sigs_labs, group_var = "HER2_subtype", group_order = subtypes_order)
p

ggsave(p, file = "results/figs/HRD/CALGB/dotplot_HER2_subtypes.png", width = 6, height = 6)

pam50_levels <- c("Basal", "Her2", "LumB", "LumA", "Normal")
p <- signature_dotplot(df = merged_data, signatures = sigs, signatures_labs = sigs_labs, group_var = "PAM50", group_order = pam50_levels) ## ! AIMS
ggsave(p, file = "results/figs/HRD/CALGB/dotplot_PAM50_subtypes.png", width = 6, height = 6)

### ALTTO
merged_data <- readRDS("data/HRD_signatures/ALTTO_meta.rds")

subtypes_order <- c("IM", "P/Met", "Mes/S", "LUM", "ERBB2-D")
p <- signature_dotplot(df = merged_data, signatures = sigs, signatures_labs = sigs_labs, group_var = "HER2_subtype", group_order = subtypes_order)
p

ggsave(p, file = "results/figs/HRD/ALTTO/dotplot_HER2_subtypes.png", width = 6, height = 6)

pam50_levels <- c("Basal", "Her2", "LumB", "LumA", "Normal")
p <- signature_dotplot(df = merged_data, signatures = sigs, signatures_labs = sigs_labs, group_var = "AIMS", group_order = pam50_levels) ## ! AIMS
ggsave(p, file = "results/figs/HRD/ALTTO/dotplot_PAM50_subtypes.png", width = 6, height = 6)


#####################
### HER2 subtypes ###
#####################

## Computation pipeline - CALGB
dir <- "/mnt/bctl/lcer0007/HER2_subtypes/"

sigs_groups_class_tot <- readRDS(paste0(dir, "input_files/sigs_groups_class_final.RDS"))
median_genes <- readRDS(paste0(dir, "input_files/median_genes.RDS"))
xcenter_genes <- readRDS(paste0(dir, "input_files/x_mean_genes.RDS"))
xcenter_sd_genes <- readRDS(paste0(dir, "input_files/x_sd_genes.RDS"))

library(matrixStats)

source(paste0(dir, "calc_HER2_groups_function.R"))

HER2_subtypes_obj <- calc_HER2_groups(tpm, type = "TPM", sig = sigs_groups_class_tot, median_genes = median_genes, center_genes = xcenter_genes, sd_genes = xcenter_sd_genes)

merged_data$HER2_subtype <- HER2_subtypes_obj$HER2_subtype 
# calgb$HER2_subtype <- HER2_subtypes_obj$HER2_subtype 

table(merged_data$HER2_subtype)

subtype_colors <- readRDS("data/color_palettes/HER2_subtypes.rds")
subtype_levels <- names(subtype_colors)

## neoALTTO
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "BRCA1ness", 
                         cohort_name = "neoALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_signature_Walens", 
                          signature_lab = "HRD signature (Walens et al.)",cohort_name = "neoALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_signature_Peng", 
                          signature_lab = "HRD signature (Peng et al.)",cohort_name = "neoALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_Zhuang", 
                          signature_lab = "HRD gene-pair signature (Zhuang et al.)",cohort_name = "neoALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_Beinse", 
                          signature_lab = "HRD score (Beinse et al.)",cohort_name = "neoALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)

## ALTTO
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "BRCA1ness", 
                         cohort_name = "ALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_signature_Walens", 
                          signature_lab = "HRD signature (Walens et al.)",cohort_name = "ALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_signature_Peng", 
                          signature_lab = "HRD signature (Peng et al.)",cohort_name = "ALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_Zhuang", 
                          signature_lab = "HRD gene-pair signature (Zhuang et al.)",cohort_name = "ALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_Beinse", 
                          signature_lab = "HRD score (Beinse et al.)",cohort_name = "ALTTO", custom_colors = subtype_colors, custom_levels = subtype_levels)

## CALGB
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "BRCA1ness", 
                         cohort_name = "CALGB", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_signature_Walens", 
                          signature_lab = "HRD signature (Walens et al.)",cohort_name = "CALGB", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_signature_Peng", 
                          signature_lab = "HRD signature (Peng et al.)",cohort_name = "CALGB", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_Zhuang", 
                          signature_lab = "HRD gene-pair signature (Zhuang et al.)",cohort_name = "CALGB", custom_colors = subtype_colors, custom_levels = subtype_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "HER2_subtype", variable_lab = "HER2 subtype", signature = "HRD_Beinse", 
                          signature_lab = "HRD score (Beinse et al.)",cohort_name = "CALGB", custom_colors = subtype_colors, custom_levels = subtype_levels)

#############
### PAM50 ###
#############
pam50_levels <- c("Basal", "Her2", "LumB", "LumA", "Normal")

## neoALTTO
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", signature = "BRCA1ness", cohort_name = "neoALTTO", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", signature = "HRD_signature_Walens",signature_lab = "HRD signature (Walens et al.)", cohort_name = "neoALTTO", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)",cohort_name = "neoALTTO", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)",cohort_name = "neoALTTO", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "PAM50", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)",cohort_name = "neoALTTO", custom_levels = pam50_levels)

### Split by pCR - TO DO?
pCR_data <- merged_data %>% filter(pCR == "Yes")
RD_data <- merged_data %>% filter(pCR == "No")

p1 <- clinical_variable_boxplot(pCR_data, clinical_variable = "PAM50", signature = "BRCA1ness", y_lab = "BRCA1ness", signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)
p2 <- clinical_variable_boxplot(RD_data, clinical_variable = "PAM50", signature = "BRCA1ness", y_lab = "BRCA1ness", signature_lab = "non-pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)
ggarrange(p1, p2)

clinical_variable_boxplot(pCR_data, clinical_variable = "PAM50", signature = "HRD_signature_Walens",y_lab = "HRD signature (Walens et al.)", signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)
# clinical_variable_boxplot(pCR_data, clinical_variable = "PAM50", signature = "HRD_signature_Peng", y_lab = "HRD signature (Peng et al.)", signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)
# clinical_variable_boxplot(pCR_data, clinical_variable = "PAM50", signature = "HRD_Zhuang", y_lab = "HRD gene-pair signature (Zhuang et al.)",signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)
# clinical_variable_boxplot(pCR_data, clinical_variable = "PAM50", signature = "HRD_Beinse", y_lab = "HRD score (Beinse et al.)", signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)


# clinical_variable_boxplot(RD_data, clinical_variable = "PAM50", signature = "HRD_signature_Walens",y_lab = "HRD signature (Walens et al.)", signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)
# clinical_variable_boxplot(RD_data, clinical_variable = "PAM50", signature = "HRD_signature_Peng", y_lab = "HRD signature (Peng et al.)", signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)
# clinical_variable_boxplot(RD_data, clinical_variable = "PAM50", signature = "HRD_Zhuang", y_lab = "HRD gene-pair signature (Zhuang et al.)",signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)
# clinical_variable_boxplot(RD_data, clinical_variable = "PAM50", signature = "HRD_Beinse", y_lab = "HRD score (Beinse et al.)", signature_lab = "pCR", cohort_name = "neoALTTO", custom_levels = pam50_levels, save_plot = FALSE, return_plot = TRUE)


## ALTTO
clinical_variable_boxplot(merged_data, clinical_variable = "AIMS", signature = "BRCA1ness", cohort_name = "ALTTO", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "AIMS", signature = "HRD_signature_Walens",signature_lab = "HRD signature (Walens et al.)", cohort_name = "ALTTO", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "AIMS", signature = "HRD_signature_Peng", signature_lab = "HRD signature (Peng et al.)",cohort_name = "ALTTO", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "AIMS", signature = "HRD_Zhuang", signature_lab = "HRD gene-pair signature (Zhuang et al.)",cohort_name = "ALTTO", custom_levels = pam50_levels)
clinical_variable_boxplot(merged_data, clinical_variable = "AIMS", signature = "HRD_Beinse", signature_lab = "HRD score (Beinse et al.)",cohort_name = "ALTTO", custom_levels = pam50_levels)

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


#################
### KM curves ###
#################
library(survival)
library(survminer)
library(ggpubr)

## CALGB
# Beinse et al.

## EFS
km_fit <- survfit(Surv(efsYears, efs) ~ HRD_Beinse_groups, data = merged_data)

surv_plot <- ggsurvplot(
  km_fit,
  data = merged_data,
  pval = TRUE,               # Show p-value for statistical significance
  conf.int = TRUE,           # Show confidence interval
  risk.table = TRUE,         # Add risk table below the plot
  ggtheme = theme_classic(), # Apply minimal theme
  legend.title = "HRD group (Beinse et al.)",
  legend.labs = c("High", "Low"),
  palette = c("blue", "red") # Custom colors for groups
) 
# Arrange the KM curve and risk table together
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))

# Save the full plot (KM curve + risk table)
ggsave("results/figs/HRD/CALGB/KM_EFS_Beinse.png", plot = full_km_plot, height = 7, width = 9)

## OS
km_fit <- survfit(Surv(osYears_2, os) ~ HRD_Beinse_groups, data = merged_data)

surv_plot <- ggsurvplot(
  km_fit,
  data = merged_data,
  pval = TRUE,               # Show p-value for statistical significance
  conf.int = TRUE,           # Show confidence interval
  risk.table = TRUE,         # Add risk table below the plot
  ggtheme = theme_classic(), # Apply minimal theme
  legend.title = "HRD group (Beinse et al.)",
  legend.labs = c("High", "Low"),
  palette = c("blue", "red") # Custom colors for groups
) 
# Arrange the KM curve and risk table together
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))

# Save the full plot (KM curve + risk table)
ggsave("results/figs/HRD/CALGB/KM_OS_Beinse.png", plot = full_km_plot, height = 7, width = 9)

## neoALTTO
# Beinse et al.

## EFS
# efs contains the Surv formula format already
km_fit <- survfit(efs ~ HRD_Beinse_groups, data = merged_data)

surv_plot <- ggsurvplot(
  km_fit,
  data = merged_data,
  pval = TRUE,               # Show p-value for statistical significance
  conf.int = TRUE,           # Show confidence interval
  risk.table = TRUE,         # Add risk table below the plot
  ggtheme = theme_classic(), # Apply minimal theme
  legend.title = "HRD group (Beinse et al.)",
  legend.labs = c("High", "Low"),
  palette = c("blue", "red") # Custom colors for groups
) 
# Arrange the KM curve and risk table together
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))

# Save the full plot (KM curve + risk table)
ggsave("results/figs/HRD/neoALTTO/KM_EFS_Beinse.png", plot = full_km_plot, height = 7, width = 9)

## OS
km_fit <- survfit(os ~ HRD_Beinse_groups, data = merged_data)

surv_plot <- ggsurvplot(
  km_fit,
  data = merged_data,
  pval = TRUE,               # Show p-value for statistical significance
  conf.int = TRUE,           # Show confidence interval
  risk.table = TRUE,         # Add risk table below the plot
  ggtheme = theme_classic(), # Apply minimal theme
  legend.title = "HRD group (Beinse et al.)", 
  legend.labs = c("High", "Low"),
  palette = c("blue", "red") # Custom colors for groups
) 
# Arrange the KM curve and risk table together
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))

# Save the full plot (KM curve + risk table)
ggsave("results/figs/HRD/neoALTTO/KM_OS_Beinse.png", plot = full_km_plot, height = 7, width = 9)


### Distribution metrics across datasets

merged_data <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
# median(merged_data$BRCA1ness)

# wb <- createWorkbook()

dataset <- "TNBC"

df <- data.frame(signature = character(), mean = numeric(), median = numeric(), sd = numeric())

for(sig in signatures){
  df <- rbind(df, data.frame(signature = sig, mean = mean(merged_data[[sig]]), median = median(merged_data[[sig]]), sd = sd(merged_data[[sig]])))
}
# ### Save to Excel file
library(openxlsx)

wb <- loadWorkbook("results/HRD/distributions_stats.xlsx")

addWorksheet(wb, dataset)
writeData(wb, sheet = dataset, x = df)

saveWorkbook(wb, "results/HRD/distributions_stats.xlsx", overwrite = TRUE)


##### Compare distributions of different datasets
tnbc <- merged_data

merged_data <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
altto <- merged_data

merged_data <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
neoaltto <- merged_data

merged_data <- readRDS("data/HRD_signatures/CALGB_meta.rds")
calgb <- merged_data

altto$study <- "ALTTO"
neoaltto$study <- "neoALTTO"
calgb$study <- "CALGB"
tnbc$study <- "TNBC"

altto <- altto[,c("study", signatures)]
neoaltto <- neoaltto[,c("study", signatures)]
calgb <- calgb[,c("study", signatures)]
tnbc <- tnbc[,c("study", signatures)]

merged_her2 <- rbind(altto, neoaltto, calgb)
merged_all <- rbind(merged_her2, tnbc)

kruskal.test(BRCA1ness ~ study, data = merged_all)

# ggplot(merged_all, aes(x = study, y = HRD_signature_Peng)) + geom_boxplot()

clinical_variable_boxplot(merged_all, clinical_variable = "study", variable_lab = "Dataset", signature = "BRCA1ness", cohort_name = "all")
clinical_variable_boxplot(merged_all, clinical_variable = "study", variable_lab = "Dataset", signature = "HRD_signature_Peng", signature_lab = "Peng et al.",cohort_name = "all")
clinical_variable_boxplot(merged_all, clinical_variable = "study", variable_lab = "Dataset", signature = "HRD_signature_Walens", signature_lab = "Walens et al.",cohort_name = "all")
clinical_variable_boxplot(merged_all, clinical_variable = "study", variable_lab = "Dataset", signature = "HRD_Beinse_scaled", signature_lab = "Beinse et al.",cohort_name = "all")
clinical_variable_boxplot(merged_all, clinical_variable = "study", variable_lab = "Dataset", signature = "HRD_Zhuang_scaled", signature_lab = "Zhuang et al.",cohort_name = "all")
