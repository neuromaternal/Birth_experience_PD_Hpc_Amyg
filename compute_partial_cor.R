compute_partial_cor <- function(x, y, covariates, corr_method = "pearson") {
  # Combine the variables into a single matrix
  data_matrix <- cbind(x, y, covariates)
  # Compute the partial correlation
  pcor_result <- ppcor::pcor(data_matrix, method = corr_method)
  # Extract the partial correlation coefficient and p-value
  list(cor = pcor_result$estimate[2, 1], p.value = pcor_result$p.value[2, 1])
}