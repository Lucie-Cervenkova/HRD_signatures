library(survival)
library(broom)
library(ggplot2)
library(dplyr)

make_forest_plot <- function(data, signatures, surv_obj = NULL, time = NULL, event = NULL, covariates = NULL, signatures_labs = NULL, pos_lab = NULL, neg_lab = NULL, 
                             binary = FALSE, observed_variable = NULL, return_object = FALSE, title = "") {
  results <- data.frame()
  
  for (sig in signatures) {
    if (binary) { # Logistic regression
      formula_str <- paste(observed_variable, "~", sig)
      if (!is.null(covariates)) {
        formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
      }
      model <- glm(as.formula(formula_str), data = data, family = binomial)
      tidy_model <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
      sig_row <- tidy_model %>% filter(term == sig)
      sig_row$signature <- sig
      results <- rbind(results, sig_row)

      if(is.null(neg_lab) & is.null(pos_lab)){
        pos_lab <- "Favors pCR"
        neg_lab <- "Favors non-pCR"
      }
    } else { # Cox regression
      if (!is.null(time) & !is.null(event)) {
        surv_obj <- Surv(time = data[[time]], event = data[[event]])
      }
      formula_str <- paste("surv_obj ~", sig)
      if (!is.null(covariates)) {
        formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
      }
      model <- coxph(as.formula(formula_str), data = data)
      tidy_model <- broom::tidy(model, exponentiate = TRUE, conf.int = TRUE)
      sig_row <- tidy_model %>% filter(term == sig)
      sig_row$signature <- sig
      results <- rbind(results, sig_row)

      neg_lab <- "Better outcome"
      pos_lab <- "Worse outcome"
    }
  }
  
  results <- results %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
    mutate(color = case_when(
      p.adj >= 0.05 ~ "grey",
      !binary & p.adj < 0.05 & estimate < 1 ~ "blue",
      !binary & p.adj < 0.05 & estimate > 1 ~ "red",
      binary & p.adj < 0.05 & estimate > 1 ~ "blue",
      binary & p.adj < 0.05 & estimate < 1 ~ "red"
    ))

  # View(results)

  min_HR <- min(results$conf.low, na.rm = TRUE)
  max_HR <- max(results$conf.high, na.rm = TRUE)

  # Positions for labels
  y_top <- length(signatures) + 0.5   # Put text above the top row
  left_x <- 1 / 1.15                   # Left of OR = 1
  right_x <- 1 * 1.15                  # Right of OR = 1
  
  p <- ggplot(results, aes(x = signature, y = estimate, ymin = conf.low, ymax = conf.high)) +
        geom_pointrange(aes(color = color), shape = 15, size = 1) +
        geom_hline(yintercept = 1, linetype = "dashed") +
        coord_flip() +
        scale_y_log10(limits = c(min_HR * 0.8, max_HR * 1.2)) +
        scale_color_identity() +
        scale_x_discrete(labels = signatures_labs) +
        # Here's the fixed annotation perfectly at the top, around the OR = 1
        annotate("text", x = y_top, y = left_x, label = neg_lab, hjust = 1, size = 4) +
        annotate("text", x = y_top, y = right_x, label = pos_lab, hjust = 0, size = 4) +
        labs(y = ifelse(binary, "Odds Ratio", "Hazard Ratio"), x = "", color = "FDR < 0.05") +
        theme_bw(base_size = 14) +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  
  if(return_object){
    df <- results %>% select(signature, estimate, p.value, p.adj)
    df$clinical_variable <- observed_variable
    df <- df %>% select(clinical_variable, everything())
    df$signature <- signatures_labs
    colnames(df) <- c("clinical_variable", "signature", "odds_ratio", "p_value", "p_adj")

    return(list(plot = p, data = df))
  }
  else {
    return(p)
  }
}

make_multi_group_forest_plot <- function(data,
                             signatures,
                             group_var,
                             group_levels,
                             surv_time    = NULL,
                             surv_event   = NULL,
                             covariates   = NULL,
                             signatures_labs = NULL,
                             pos_lab      = NULL,
                             neg_lab      = NULL,
                             binary       = FALSE,
                             observed_variable = NULL,
                             title        = "") {
  
  results <- lapply(group_levels, function(g) {
    df_g <- filter(data, !!sym(group_var) == g)
    
    lapply(signatures, function(sig) {
      if (binary) {
        # logistic
        f_str <- paste(observed_variable, "~", sig,
                       if (!is.null(covariates)) paste("+", paste(covariates, collapse = " + ")) else "")
        mdl <- glm(as.formula(f_str), data = df_g, family = binomial)
      } else {
        # Cox
        if (!is.null(surv_time) && !is.null(surv_event)) {
          surv_obj <- Surv(time = df_g[[surv_time]], event = df_g[[surv_event]])
        }
        f_str <- paste("surv_obj ~", sig,
                       if (!is.null(covariates)) paste("+", paste(covariates, collapse = " + ")) else "")
        mdl <- coxph(as.formula(f_str), data = df_g)
      }
      
      broom::tidy(mdl, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term == sig) %>%
        transmute(
          signature = sig,
          estimate  = estimate,
          conf.low  = conf.low,
          conf.high = conf.high,
          p.value   = p.value,
          cohort    = g
        )
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  dodge_width <- 0.7
  n_groups <- length(group_levels)
  palette  <- RColorBrewer::brewer.pal(n = n_groups, name = "Accent")
  
  if (is.null(pos_lab)) pos_lab <- if (binary) "Favors case" else "Worse outcome"
  if (is.null(neg_lab)) neg_lab <- if (binary) "Favors control" else "Better outcome"
  
  ggplot(results,
         aes(x = signature,
             y = estimate,
             ymin = conf.low,
             ymax = conf.high,
             color = factor(cohort, levels = group_levels))) +
    geom_pointrange(position = position_dodge(dodge_width), size = 1, shape = 15) +
    geom_vline(
       xintercept = seq(1.5, length(signatures) - 0.5, by = 1),
       color      = "grey",
       size       = 0.3
     ) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    coord_flip() +
    scale_y_log10() +
    scale_color_manual(values = palette, name = "Treatment arm") +
    # scale_shape_manual(values = shapes, name = group_var) +
    scale_x_discrete(labels = signatures_labs) +
    # annotations at the top for directionality
    annotate("text",
             x = length(signatures) + 0.5,
             y = 1 / 1.15,
             label = neg_lab,
             hjust = 1, size = 4) +
    annotate("text",
             x = length(signatures) + 0.5,
             y = 1 * 1.15,
             label = pos_lab,
             hjust = 0, size = 4) +
    labs(
      x     = NULL,
      y     = ifelse(binary, "Odds Ratios", "Hazard Ratios"),
      title = title
    ) +
    guides(
    color = guide_legend(
      title.position = "top",
      title.hjust     = 0.5)) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical"
    )
}