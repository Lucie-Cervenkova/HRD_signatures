library(dplyr)
library(ggplot2)
library(ggpubr)

sigs_HRD <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Beinse" =  "HRD_Beinse_scaled", "Zhuang" = "HRD_Zhuang_scaled")

neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")

source("scripts/HRD/forest_plots_pCR_and_survival.R")

####################################### 
## Clinical variables - forest plots ##
#######################################
clin_var_df <- data.frame()
neoaltto_clin_var <- data.frame()
altto_clin_var <- data.frame()
calgb_clin_var <- data.frame()

#### neoALTTO
neoaltto <- neoaltto %>% mutate(HR_positive = ifelse(erpgrstatus == "NEGATIVE", 0, 1))
neoaltto <- neoaltto %>% mutate(ER_positive = ifelse(erstatus == "NEGATIVE", 0, 1))
neoaltto <- neoaltto %>% mutate(PR_positive = ifelse(pgrstatus == "NEGATIVE", 0, 1))

# 1. HR
res <- make_forest_plot(
  data = neoaltto,
  signatures = sigs_HRD,
  observed_variable = "HR_positive",
  neg_lab = "HR-",
  pos_lab = "HR+",
  signatures_labs = names(sigs_HRD),
  binary = TRUE,
  return_object = TRUE,
)$data
neoaltto_clin_var <- rbind(neoaltto_clin_var, res)

# 2. ER
res <- make_forest_plot(
  data = neoaltto,
  signatures = sigs_HRD,
  observed_variable = "ER_positive",
  neg_lab = "ER-",
  pos_lab = "ER+",    
  signatures_labs = names(sigs_HRD),
  binary = TRUE,
  return_object = TRUE,
)$data
neoaltto_clin_var <- rbind(neoaltto_clin_var, res)

# 3. PR
res <- make_forest_plot(
  data = neoaltto,
  signatures = sigs_HRD,
  observed_variable = "PR_positive",
  neg_lab = "PR-",
  pos_lab = "PR+",    
  signatures_labs = names(sigs_HRD),
  binary = TRUE,
  return_object = TRUE,
)$data
neoaltto_clin_var <- rbind(neoaltto_clin_var, res)

# 4. Nodes
nodes_data <- neoaltto %>% filter(nstage != "NX")
nodes_data <- nodes_data %>% mutate(node_positive = ifelse(nstage == "N0", 0, 1))

res <- make_forest_plot(
  data = nodes_data,
  signatures = sigs_HRD,
  observed_variable = "node_positive",
  neg_lab = "N0",
  pos_lab = "N+",
  signatures_labs = names(sigs_HRD),
  binary = TRUE, # logistic regression
  return_object = TRUE
)$data
neoaltto_clin_var <- rbind(neoaltto_clin_var, res)

# 5. Grade
histo_data <- neoaltto %>% filter(!(histograde %in% c("GX","NK")))
histo_data <- histo_data %>% mutate(grade = ifelse(histograde == "G3", 1, 0))

res <- make_forest_plot(
  data = histo_data,
  signatures = sigs_HRD,
  observed_variable = "grade",
  neg_lab = "G1-2",
  pos_lab = "G3",
  signatures_labs = names(sigs_HRD),
  binary = TRUE, # logistic regression
  return_object = TRUE
)$data
neoaltto_clin_var <- rbind(neoaltto_clin_var, res)

# 6. Age
neoaltto <- neoaltto %>% mutate(age = ifelse(Age >= 40, 1, 0)) # v2: 40, v1: 50
table(neoaltto$age)

res <- make_forest_plot(
  data = neoaltto,
  signatures = sigs_HRD,
  observed_variable = "age",
  neg_lab = "<40",
  pos_lab = ">=40",
  signatures_labs = names(sigs_HRD),
  binary = TRUE,
  return_object = TRUE
)$data
neoaltto_clin_var <- rbind(neoaltto_clin_var, res)

# 7. Menopause
menop_data <- neoaltto %>% filter(menop %in% c("PREMENOPAUSAL", "POSTMENOPAUSAL"))
menop_data <- menop_data %>% mutate(menopause = ifelse(menop ==  "POSTMENOPAUSAL", 1, 0))

res <- make_forest_plot(
  data = menop_data,
  signatures = sigs_HRD,
  observed_variable = "menopause",
  neg_lab = "Premenopausal",
  pos_lab = "Postmenopausal",
  signatures_labs = names(sigs_HRD),
  binary = TRUE, 
  return_object = TRUE
)$data
neoaltto_clin_var <- rbind(neoaltto_clin_var, res)

neoaltto_clin_var$study <- "NeoALTTO"

### CALGB
calgb <- calgb %>% mutate(HR_positive = ifelse(HR_reviewed == "neg", 0, 1))
calgb <- calgb %>% mutate(ER_positive = ifelse(ER_reviewed == "neg", 0, 1))
calgb <- calgb %>% mutate(PR_positive = ifelse(PR_reviewed == "neg", 0, 1))

# 1. HR
res <- make_forest_plot(
  data = calgb,
  signatures = sigs_HRD,
  observed_variable = "HR_positive",
  neg_lab = "HR-",
  pos_lab = "HR+",
  signatures_labs = names(sigs_HRD),
  binary = TRUE, # logistic regression
    return_object = TRUE
)$data

calgb_clin_var <- rbind(calgb_clin_var, res)

# 2. ER
res <- make_forest_plot(
  data = calgb,
  signatures = sigs_HRD,
  observed_variable = "ER_positive",
  neg_lab = "ER-",
  pos_lab = "ER+",
  signatures_labs = names(sigs_HRD),
  binary = TRUE, # logistic regression
    return_object = TRUE
)$data

calgb_clin_var <- rbind(calgb_clin_var, res)

# 3. PR
res <- make_forest_plot(
  data = calgb,
  signatures = sigs_HRD,
  observed_variable = "PR_positive",
  neg_lab = "PR-",
  pos_lab = "PR+",
  signatures_labs = names(sigs_HRD),
  binary = TRUE, # logistic regression
    return_object = TRUE
)$data

calgb_clin_var <- rbind(calgb_clin_var, res)

# 4. Nodes
nodes_data <- calgb %>% filter(!(Clinical_N_Stage %in% c("NA", NA)))
nodes_data <- nodes_data %>% mutate(nodes = ifelse(Clinical_N_Stage == "N0", 0, 1))

res <- make_forest_plot(
  data = nodes_data,
  signatures = sigs_HRD,
  observed_variable = "nodes",
  neg_lab = "N0",
  pos_lab = "N+",
  signatures_labs = names(sigs_HRD),
  binary = TRUE, # logistic regression
  return_object = TRUE
)$data

calgb_clin_var <- rbind(calgb_clin_var, res)

# 5. Age
calgb <- calgb %>% mutate(age = ifelse(Age_on_Study >= 40, 1, 0))
table(calgb$age)

res <- make_forest_plot(
  data = calgb,
  signatures = sigs_HRD,
  observed_variable = "age",
  neg_lab = "<40",
  pos_lab = ">=40",
  signatures_labs = names(sigs_HRD),
  binary = TRUE,
  return_object = TRUE
)$data
calgb_clin_var <- rbind(calgb_clin_var, res)
calgb_clin_var$study <- "CALGB"

### ALTTO
altto <- altto %>% mutate(HR_positive = ifelse(hr == "Negative", 0, 1))

# 1. HR
res <- make_forest_plot(
  data = altto,
  signatures = sigs_HRD,
  observed_variable = "HR_positive",
  neg_lab = "HR-",
  pos_lab = "HR+",
  signatures_labs = names(sigs_HRD),
  binary = TRUE, 
  return_object = TRUE
)$data
altto_clin_var <- rbind(altto_clin_var, res)

# 2. Grade
histo_data <- altto %>% filter(hgrade != "GX: Differentiation cannot be assessed")
histo_data <- histo_data %>% mutate(grade = ifelse(hgrade == "G3: Poorly differentiated/Undifferentiated", 1, 0))

res <- make_forest_plot(
  data = histo_data,
  signatures = sigs_HRD,
  observed_variable = "grade",
  neg_lab = "G1-2",
  pos_lab = "G3",
  signatures_labs = names(sigs_HRD),
  binary = TRUE,
  return_object = TRUE
)$data

altto_clin_var <- rbind(altto_clin_var, res)

# 4. Nodes
altto <- altto %>% mutate(nodes = ifelse(nodal_status == "Node Negative", 0, 1))

res <- make_forest_plot(
  data = altto,
  signatures = sigs_HRD,
  observed_variable = "nodes",
  neg_lab = "N0",
  pos_lab = "N+",
  signatures_labs = names(sigs_HRD),
  binary = TRUE,
  return_object = TRUE
)$data

altto_clin_var <- rbind(altto_clin_var, res)

# 5. Age
altto <- altto %>% mutate(age = ifelse(age >= 40, 1, 0))

res <- make_forest_plot(
  data = altto,
  signatures = sigs_HRD,
  observed_variable = "age",
  neg_lab = "<40",
  pos_lab = ">=40",
  signatures_labs = names(sigs_HRD),
  binary = TRUE,
  return_object = TRUE
)$data

altto_clin_var <- rbind(altto_clin_var, res)
altto_clin_var$study <- "ALTTO"


clin_var_df <- rbind(neoaltto_clin_var, calgb_clin_var, altto_clin_var)
View(clin_var_df)
write.csv(clin_var_df, "results/figs/HRD/Table2_clinical_variables.csv", row.names = FALSE)
