perform_factorization <- function(method, data_matrix, rank, alpha) {
  if (method == "stICA") {
    result <- unbiased_stICA(data_matrix, rank, alpha = alpha)
  } else if (method == "SVD") {
    result <- svd(data_matrix, nu = rank, nv = rank)
    result$A <- result$u %*% diag(result$d[1:rank], rank)
    result$B <- result$v
  } else {
    stop("Unsupported factorization method")
  }
  return(result)
}

get_batch_indices <- function(batch_R2, threshold) {
  if (threshold < 0 || threshold > 1) {
    stop("batch_threshold not in [0, 1]")
  }
  batch_indices <- which(batch_R2$allpv < threshold)
  return(batch_indices)
}

get_group_indices <- function(group_R2, group_types, group_threshold) {
  if (any(group_threshold < 0 | group_threshold > 1)) {
    stop("group_threshold not in [0, 1]")
  }
  
  group_indices <- integer(0)
  
  if (length(group_threshold) != length(group_types)) {
    if (length(group_threshold) == 1) {
      group_threshold <- rep(group_threshold, length(group_types))
    } else {
      stop("length(group_threshold) should be equal to 1 or length(group_types)")
    }
  }

  for (i in seq_along(group_types)) {
    group_indices <- c(group_indices, which(group_R2$allpv[, i] < group_threshold[i]))
  }

  return(group_indices)
}

remove_components <- function(data_matrix, factor_matrix_A, factor_matrix_B, batch_indices) {
  components_to_remove_A <- factor_matrix_A[, batch_indices]
  components_to_remove_B <- factor_matrix_B[, batch_indices]
  print(paste("Removing", length(batch_indices), "components with P value less than batch_threshold"))
  normalized_matrix <- data_matrix - components_to_remove_A %*% t(components_to_remove_B)
  return(normalized_matrix)
}

combine_R2_results <- function(batch_R2, group_matrix, group_R2) {
  combined_R2 <- batch_R2$allR2
  if (!is.null(group_matrix)) {
    combined_R2 <- cbind(combined_R2, group_R2$allR2)
  }
  return(combined_R2)
}

#' Normalize Data by Removing Batch Effects
#'
#' Function to normalize data matrix `X` by factorizing `X = A %*% t(B)` and removing components with a high R² value with respect to batch effects.
#' Components with high R² values with respect to group effects are retained.
#'
#' @param factorization_method Character. The factorization method, either 'SVD' or 'stICA'.
#' @param data_matrix Matrix (n x p). Samples x features matrix to normalize.
#' @param batch_vector Vector (n). Variable representing the batch information to remove from `data_matrix`.
#' @param batch_type Character. Type of `batch_vector`, either 'categorical' or 'continuous'.
#' @param rank Integer. Rank of the low-rank decomposition. Default is 20.
#' @param batch_threshold Numeric. Threshold in [0,1]; if R²(component, batch_vector) > batch_threshold, the component is removed. Default is 0.5.
#' @param group_matrix Matrix (n x l), optional. Each column represents a group variable to retain in `data_matrix`.
#' @param group_types Character vector. Each element corresponds to the type of `group_matrix` columns, either 'categorical' or 'continuous'.
#' @param group_threshold Numeric vector. Threshold(s) in [0,1]; if R²(component, group_vector) > group_threshold[i], the component is retained. If scalar, this threshold applies to all group variables. Default is 0.5.
#' @param alpha Numeric, optional. Parameter for `stICA` factorization.
#' @return A list with the following components:
#' \item{normalized_matrix}{Matrix (n x p). Normalized version of `data_matrix`.}
#' \item{R2_matrix}{Matrix. R² values between `B` components and batch/group variables.}
#' \item{removed_components}{Matrix. Components of `B` correlating with `batch_vector` but not with `group_matrix`, removed from `data_matrix`.}
#' \item{factor_matrix_A}{Matrix. `A` in the matrix factorization `X = A %*% t(B)`.}
#' \item{factor_matrix_B}{Matrix. `B` in the matrix factorization `X = A %*% t(B)`.}
#' @references Renard E., Branders S., Absil P.-A.: Independent Component Analysis to Remove Batch Effects from Merged Microarray Datasets (WABI2016).
normFact <- function(
  factorization_method,
  data_matrix,
  batch_vector,
  batch_type,
  rank = 20,
  batch_threshold = 0.5,
  group_matrix = NULL,
  group_types = NULL,
  group_threshold = 0.5,
  alpha = NULL
) {
  factorization_result <- perform_factorization(factorization_method, data_matrix, rank, alpha)
  factor_matrix_A <- factorization_result$A
  factor_matrix_B <- factorization_result$B

  batch_R2 <- R2(batch_vector, factor_matrix_B, batch_type, pval = TRUE)
  batch_indices <- get_batch_indices(batch_R2, batch_threshold)

  if (!is.null(group_matrix)) {
    group_R2 <- R2(group_matrix, factor_matrix_B, group_types, pval = TRUE)
    group_indices <- get_group_indices(group_R2, group_types, group_threshold)
    indices_to_keep <- intersect(batch_indices, group_indices)
    print(paste("Keeping", length(indices_to_keep), "components with P value less than group_threshold"))
    batch_indices <- setdiff(batch_indices, indices_to_keep)
  }

  normalized_matrix <- remove_components(data_matrix, factor_matrix_A, factor_matrix_B, batch_indices)

  combined_R2 <- combine_R2_results(batch_R2, group_matrix, group_R2)

  list(
    normalized_matrix = normalized_matrix,
    R2_matrix = combined_R2,
    removed_components = factor_matrix_B[, batch_indices],
    factor_matrix_A = factor_matrix_A,
    factor_matrix_B = factor_matrix_B
  )
}
