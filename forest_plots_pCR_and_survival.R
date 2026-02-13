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
    return(list(plot = p, data = results))
  }
  else {
    return(p)
  }
}
