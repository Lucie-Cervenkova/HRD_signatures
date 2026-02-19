### This is in separate script because I'm unable to pass different variable as string in function to survfit and I already wasted too much time fixing it. 
library(dplyr)
library(ggplot2)
# library(stringr)
# library(rlang)
library(ggpubr)
library(survival)
library(survminer)

## 1. Load data
sigs_cols <- list("BRCA1ness" = "BRCA1ness", "Walens" = "HRD_signature_Walens", "Peng" = "HRD_signature_Peng", "Zhuang" = "HRD_Zhuang_scaled", "Beinse" =  "HRD_Beinse_scaled")

neoaltto <- readRDS("data/HRD_signatures/neoALTTO_meta.rds")
altto <- readRDS("data/HRD_signatures/ALTTO_meta.rds")
calgb <- readRDS("data/HRD_signatures/CALGB_meta.rds")

## 2. Split by median and or quartiles
neoaltto <- neoaltto %>% mutate(HRD_Beinse_median_groups = ifelse(HRD_Beinse_scaled >= median(HRD_Beinse_scaled, na.rm = TRUE), "High", "Low")) %>% 
                         mutate(HRD_Beinse_quartile_groups = ifelse(HRD_Beinse_scaled >= quantile(HRD_Beinse_scaled,0.75, na.rm = TRUE), "High", ifelse(HRD_Beinse_scaled < quantile(HRD_Beinse_scaled,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(BRCA1ness_median_groups = ifelse(BRCA1ness >= median(BRCA1ness, na.rm = TRUE), "High", "Low")) %>%
                         mutate(BRCA1ness_quartile_groups = ifelse(BRCA1ness >= quantile(BRCA1ness,0.75, na.rm = TRUE), "High", ifelse(BRCA1ness < quantile(BRCA1ness,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Walens_median_groups = ifelse(HRD_signature_Walens >= median(HRD_signature_Walens, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Walens_quartile_groups = ifelse(HRD_signature_Walens >= quantile(HRD_signature_Walens,0.75, na.rm = TRUE), "High", ifelse(HRD_signature_Walens < quantile(HRD_signature_Walens,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Peng_median_groups = ifelse(HRD_signature_Peng >= median(HRD_signature_Peng, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Peng_quartile_groups = ifelse(HRD_signature_Peng >= quantile(HRD_signature_Peng,0.75, na.rm = TRUE), "High", ifelse(HRD_signature_Peng < quantile(HRD_signature_Peng,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Zhuang_median_groups = ifelse(HRD_Zhuang_scaled >= median(HRD_Zhuang_scaled, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Zhuang_quartile_groups = ifelse(HRD_Zhuang_scaled >= quantile(HRD_Zhuang_scaled,0.75, na.rm = TRUE), "High", ifelse(HRD_Zhuang_scaled < quantile(HRD_Zhuang_scaled,0.25, na.rm = TRUE), "Low", "Medium")))

altto <- altto %>% mutate(HRD_Beinse_median_groups = ifelse(HRD_Beinse_scaled >= median(HRD_Beinse_scaled, na.rm = TRUE), "High", "Low")) %>% 
                         mutate(HRD_Beinse_quartile_groups = ifelse(HRD_Beinse_scaled >= quantile(HRD_Beinse_scaled,0.75, na.rm = TRUE), "High", ifelse(HRD_Beinse_scaled < quantile(HRD_Beinse_scaled,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(BRCA1ness_median_groups = ifelse(BRCA1ness >= median(BRCA1ness, na.rm = TRUE), "High", "Low")) %>%
                         mutate(BRCA1ness_quartile_groups = ifelse(BRCA1ness >= quantile(BRCA1ness,0.75, na.rm = TRUE), "High", ifelse(BRCA1ness < quantile(BRCA1ness,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Walens_median_groups = ifelse(HRD_signature_Walens >= median(HRD_signature_Walens, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Walens_quartile_groups = ifelse(HRD_signature_Walens >= quantile(HRD_signature_Walens,0.75, na.rm = TRUE), "High", ifelse(HRD_signature_Walens < quantile(HRD_signature_Walens,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Peng_median_groups = ifelse(HRD_signature_Peng >= median(HRD_signature_Peng, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Peng_quartile_groups = ifelse(HRD_signature_Peng >= quantile(HRD_signature_Peng,0.75, na.rm = TRUE), "High", ifelse(HRD_signature_Peng < quantile(HRD_signature_Peng,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Zhuang_median_groups = ifelse(HRD_Zhuang_scaled >= median(HRD_Zhuang_scaled, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Zhuang_quartile_groups = ifelse(HRD_Zhuang_scaled >= quantile(HRD_Zhuang_scaled,0.75, na.rm = TRUE), "High", ifelse(HRD_Zhuang_scaled < quantile(HRD_Zhuang_scaled,0.25, na.rm = TRUE), "Low", "Medium")))

calgb <- calgb %>% mutate(HRD_Beinse_median_groups = ifelse(HRD_Beinse_scaled >= median(HRD_Beinse_scaled, na.rm = TRUE), "High", "Low")) %>% 
                         mutate(HRD_Beinse_quartile_groups = ifelse(HRD_Beinse_scaled >= quantile(HRD_Beinse_scaled,0.75, na.rm = TRUE), "High", ifelse(HRD_Beinse_scaled < quantile(HRD_Beinse_scaled,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(BRCA1ness_median_groups = ifelse(BRCA1ness >= median(BRCA1ness, na.rm = TRUE), "High", "Low")) %>%
                         mutate(BRCA1ness_quartile_groups = ifelse(BRCA1ness >= quantile(BRCA1ness,0.75, na.rm = TRUE), "High", ifelse(BRCA1ness < quantile(BRCA1ness,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Walens_median_groups = ifelse(HRD_signature_Walens >= median(HRD_signature_Walens, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Walens_quartile_groups = ifelse(HRD_signature_Walens >= quantile(HRD_signature_Walens,0.75, na.rm = TRUE), "High", ifelse(HRD_signature_Walens < quantile(HRD_signature_Walens,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Peng_median_groups = ifelse(HRD_signature_Peng >= median(HRD_signature_Peng, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Peng_quartile_groups = ifelse(HRD_signature_Peng >= quantile(HRD_signature_Peng,0.75, na.rm = TRUE), "High", ifelse(HRD_signature_Peng < quantile(HRD_signature_Peng,0.25, na.rm = TRUE), "Low", "Medium"))) %>% 
                         mutate(Zhuang_median_groups = ifelse(HRD_Zhuang_scaled >= median(HRD_Zhuang_scaled, na.rm = TRUE), "High", "Low")) %>%
                         mutate(Zhuang_quartile_groups = ifelse(HRD_Zhuang_scaled >= quantile(HRD_Zhuang_scaled,0.75, na.rm = TRUE), "High", ifelse(HRD_Zhuang_scaled < quantile(HRD_Zhuang_scaled,0.25, na.rm = TRUE), "Low", "Medium")))

## 3. KM plots

## a) CALGB
### Peng
# Median
km_fit <- survfit(Surv(efsYears, efs) ~ Peng_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_Peng_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ Peng_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_Peng_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
calgb_subset <- calgb %>% filter(Peng_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(efsYears, efs) ~ Peng_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_Peng_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ Peng_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_Peng_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### BRCA1ness
# Median
km_fit <- survfit(Surv(efsYears, efs) ~ BRCA1ness_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_BRCA1ness_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ BRCA1ness_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_BRCA1ness_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
calgb_subset <- calgb %>% filter(BRCA1ness_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(efsYears, efs) ~ BRCA1ness_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_BRCA1ness_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ BRCA1ness_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_BRCA1ness_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Walens
# Median
km_fit <- survfit(Surv(efsYears, efs) ~ Walens_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_Walens_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ Walens_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_Walens_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
calgb_subset <- calgb %>% filter(Walens_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(efsYears, efs) ~ Walens_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_Walens_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ Walens_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_Walens_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Beinse
# Median
km_fit <- survfit(Surv(efsYears, efs) ~ Beinse_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_Beinse_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ Beinse_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_Beinse_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
calgb_subset <- calgb %>% filter(Beinse_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(efsYears, efs) ~ Beinse_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_Beinse_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ Beinse_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_Beinse_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Zhuang
# Median
km_fit <- survfit(Surv(efsYears, efs) ~ Zhuang_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_Zhuang_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ Zhuang_median_groups, data = calgb)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_Zhuang_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
calgb_subset <- calgb %>% filter(Zhuang_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(efsYears, efs) ~ Zhuang_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_EFS_Zhuang_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(osYears_2, os) ~ Zhuang_quartile_groups, data = calgb_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = calgb_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/CALGB_OS_Zhuang_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)


## b) ALTTO
### Peng
# Median
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ Peng_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_Peng_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(os_years, os_event) ~ Peng_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_Peng_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
altto_subset <- altto %>% filter(Peng_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ Peng_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_Peng_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(Surv(os_years, os_event) ~ Peng_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_Peng_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### BRCA1ness
# Median
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ BRCA1ness_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_BRCA1ness_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(os_years, os_event) ~ BRCA1ness_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_BRCA1ness_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
altto_subset <- altto %>% filter(BRCA1ness_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ BRCA1ness_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_BRCA1ness_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(Surv(os_years, os_event) ~ BRCA1ness_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_BRCA1ness_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Walens
# Median
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ Walens_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_Walens_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(os_years, os_event) ~ Walens_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_Walens_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
altto_subset <- altto %>% filter(Walens_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ Walens_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_Walens_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(Surv(os_years, os_event) ~ Walens_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_Walens_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Beinse
# Median
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ Beinse_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_Beinse_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(os_years, os_event) ~ Beinse_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_Beinse_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
altto_subset <- altto %>% filter(Beinse_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ Beinse_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_Beinse_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(Surv(os_years, os_event) ~ Beinse_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_Beinse_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Zhuang
# Median
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ Zhuang_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_Zhuang_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(Surv(os_years, os_event) ~ Zhuang_median_groups, data = altto)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_Zhuang_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
altto_subset <- altto %>% filter(Zhuang_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(Surv(dfs_years, dfs_event) ~ Zhuang_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_DFS_Zhuang_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(Surv(os_years, os_event) ~ Zhuang_quartile_groups, data = altto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = altto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/ALTTO_OS_Zhuang_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)


## a) NeoALTTO
### Peng
# Median
km_fit <- survfit(efs ~ Peng_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_Peng_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(os ~ Peng_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_Peng_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
neoaltto_subset <- neoaltto %>% filter(Peng_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(efs ~ Peng_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_Peng_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(os ~ Peng_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Peng quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_Peng_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### BRCA1ness
# Median
km_fit <- survfit(efs ~ BRCA1ness_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_BRCA1ness_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(os ~ BRCA1ness_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_BRCA1ness_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
neoaltto_subset <- neoaltto %>% filter(BRCA1ness_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(efs ~ BRCA1ness_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_BRCA1ness_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(os ~ BRCA1ness_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "BRCA1ness quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_BRCA1ness_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Walens
# Median
km_fit <- survfit(efs ~ Walens_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_Walens_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(os ~ Walens_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_Walens_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
neoaltto_subset <- neoaltto %>% filter(Walens_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(efs ~ Walens_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_Walens_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(os ~ Walens_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Walens quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_Walens_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Beinse
# Median
km_fit <- survfit(efs ~ Beinse_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_Beinse_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(os ~ Beinse_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_Beinse_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
neoaltto_subset <- neoaltto %>% filter(Beinse_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(efs ~ Beinse_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_Beinse_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(os ~ Beinse_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Beinse quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_Beinse_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)

### Zhuang
# Median
km_fit <- survfit(efs ~ Zhuang_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_Zhuang_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

km_fit <- survfit(os ~ Zhuang_median_groups, data = neoaltto)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang median groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_Zhuang_median_groups.png"), plot = full_km_plot, height = 7, width = 9)

# Quartiles
neoaltto_subset <- neoaltto %>% filter(Zhuang_quartile_groups != "Medium") # Exclude medium group for clearer comparison between high and low
km_fit <- survfit(efs ~ Zhuang_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_EFS_Zhuang_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
    
km_fit <- survfit(os ~ Zhuang_quartile_groups, data = neoaltto_subset)

surv_plot <- ggsurvplot(
            km_fit,
            data = neoaltto_subset,
            pval = TRUE,               # Show p-value for statistical significance
            conf.int = TRUE,           # Show confidence interval
            risk.table = TRUE,         # Add risk table below the plot
            ggtheme = theme_classic(), # Apply minimal theme
            legend.title = "Zhuang quartile groups",
            legend.labs = c("High", "Low"),
            palette = c("blue", "red") # Custom colors for groups
        ) 
        
full_km_plot <- ggarrange(surv_plot$plot, surv_plot$table, ncol = 1, heights = c(2, 1))
ggsave(paste0("results/figs/HRD/KM_plots/NeoALTTO_OS_Zhuang_quartile_groups.png"), plot = full_km_plot, height = 7, width = 9)
