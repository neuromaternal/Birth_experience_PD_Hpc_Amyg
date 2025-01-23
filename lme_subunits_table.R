#' lme_subunits_table
#'
#' This function performs an analysis of brain subunits using linear mixed-effects (LME) models.
#'
#' @param C A matrix of contrasts, where each row represents a contrast to be tested.
#' @param C_names A vector of contrast names corresponding to the rows of `C`.
#' @param X A design matrix for the fixed effects in the LME model.
#' @param df A data frame containing the input data with columns: `ID`, `subunits`, `value`, `ses`, and `Group`.
#' @param region_names A vector of subunit or region names to be analyzed.
#' @param descriptives Logical, if `TRUE`, descriptive statistics (mean and SD) will be computed and merged with the results.
#' 
#' @return A data frame summarizing the analysis, including F-statistics, p-values, signed effect sizes, and optionally, descriptive statistics.
#'
#' @details The function fits an LME model for each subunit and contrast, computes F-statistics, and adjusts p-values using the FDR method. If `descriptives` is `TRUE`, the function also computes mean and standard deviation for each subunit and merges these with the statistical results.
#'
#' @export
lme_subunits_table <- function(C, C_names, X, df, region_names, descriptives = TRUE,  group_var = NULL) {
  # Initialize summary table
  sum_table <- data.frame(matrix(NA, nrow = length(region_names) * length(C_names), ncol = 7))
  names(sum_table) <- c("subunits", "contrast", "Fstatistic", "uncP", "dfWithin", "dfBetween", "sgn")
  
  
  # Loop over subunits and contrasts
  row_idx <- 1
  for (sub in region_names) {
    # Prepare data for current subunit
    ni <- matrix(unname(table(df$ID[df$subunits == sub])), ncol = 1)
    Y <- df %>% filter(subunits == sub) %>% pull(value) %>% as.matrix()
    
    # Fit model and calculate F-statistics for each contrast
    stats <- lme_fit_FS(X, Zcols = c(1), Y, ni)
    
    for (con in 1:length(C_names)) {
      F_C <- lme_F(stats, matrix(C[con, ], nrow = 1))
      sum_table[row_idx, ] <- c(sub, C_names[con], F_C$F[1], F_C$pval[1], F_C$df[1:2], F_C$sgn)
      row_idx <- row_idx + 1
    }
  }
  
  # Adjust p-values using FDR
  sum_table <- sum_table %>%
    mutate(across(c(Fstatistic, uncP, dfWithin, dfBetween, sgn), as.numeric)) %>%
    group_by(contrast) %>%
    mutate(fdr = p.adjust(uncP, method = "fdr")) %>%
    ungroup()
  
  # Calculate partial eta squared, fdr-corrected p-values and format output
  sum_table <- sum_table %>%
    mutate(
      partial_etaSq = Fstatistic * dfWithin / (Fstatistic * dfWithin + dfBetween),
      fdr_signif = case_when(
        uncP < 0.05 & fdr > 0.05 ~ "*", 
        uncP < 0.05 & fdr < 0.05 ~ "**", 
        TRUE ~ ""),
      Fstatistic = round(Fstatistic, 2),
      sgn_effect = gsub("0*$", "", format(signif(sgn * partial_etaSq, 2), scientific = FALSE)),
      uncP = paste0(gsub("0*$", "", format(signif(uncP, 2), scientific = FALSE)), fdr_signif),
      df = paste0(dfWithin, ", ", round(dfBetween, 2))
    )
  
  # Reshape data for final table
  stats_long <- sum_table %>%
    select(subunits, contrast, Fstatistic, df, uncP, sgn_effect) %>%
    melt(id.vars = c("subunits", "contrast"), variable.name = "statistic") %>%
    mutate(statistic = as.character(statistic)) %>%
    dcast(subunits + statistic ~ contrast, value.var = "value")
  
  # Compute descriptive statistics if descriptives is TRUE
  descriptives_long <- NULL  
  if (descriptives) {
    df$Group2 <- df[[group_var]]
    descriptives_long <- df %>%
      filter(subunits %in% region_names) %>%
      group_by(Group2, subunits, ses) %>%
      summarize(Mean = round(mean(value), 1), SD = round(sd(value), 1), .groups = 'drop') %>%
      melt(id.vars = c("Group2", "subunits", "ses"), variable.name = "statistic") %>%
      mutate(statistic = as.character(statistic)) %>%
      dcast(subunits + statistic ~ Group2 * ses, value.var = "value")
  }
  
  # Merge descriptives with stats_long if descriptives_long is not NULL
  if (!is.null(descriptives_long)) {
    smtable <- merge(descriptives_long, stats_long, by = c("subunits", "statistic"), all = TRUE)
  } else {
    smtable <- stats_long
  }
  
  # Final formatting of the table
  smtable <- smtable %>%
    arrange(match(statistic, c("Mean", "SD", "Fstatistic", "df", "uncP", "sgn_effect"))) %>%
    arrange(match(subunits, c("hippocampus_L", "hippocampus_R", hippo_subunit_names,
                              "amygdala_L", "amygdala_R", amyg_subunit_names))) %>%
    separate(subunits, c("subunits", "hemi"), c("_")) %>%
    mutate(hemi = ifelse(hemi == "L", "Left", "Right")) %>%
    mutate(subunits = case_when(subunits == "SubComplex" ~ "Subicular Complex",
                           subunits == "DentateGyrus" ~ "Dentate Gyrus",
                           subunits == "Tail" ~ "Hippocampal Tail",
                           subunits == "CentroMedial" ~ "Centromedial",
                           subunits == "LateroBasal" ~ "Laterobasal",
                           subunits == "hippocampus" ~ "Whole Hippocampus",
                           subunits == "amygdala" ~ "Whole Amygdala",
                           TRUE ~ subunits),
      statistic = case_when(statistic == "Fstatistic" ~ "F-stat.",
                            statistic == "df" ~ "DoF",
                            statistic == "uncP" ~ "unc. P-value",
                            statistic == "sgn_effect" ~ "Signed effect size",
                            TRUE ~ statistic)
    ) %>%
    mutate(subunits = paste(hemi, subunits)) %>%
    select(-hemi)
  
  # Return the summary table
  return(smtable)
}

