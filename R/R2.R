linear_regression_estimation <- function(component_values, poi_values) {
  coefs <- coef(lm(component_values ~ as.numeric(poi_values)))
  estimated_values <- coefs[2] * as.numeric(poi_values) + coefs[1]
  return(estimated_values)
}

class_mean_estimation <- function(component_values, poi_values) {
  classes <- unique(poi_values)
  estimated_values <- rep(NA, length(component_values))
  for (class in classes) {
    class_indices <- which(poi_values == class)
    estimated_values[class_indices] <- mean(component_values[class_indices])
  }
  return(estimated_values)
}

#' Compute R-squared and P-values between Phenotypes of Interest and Components
#'
#' This function computes the R-squared and optionally the p-values between phenotypes of interest (POI) and components.
#'
#' @param poi Matrix (p x l). Representing the phenotypes of interest.
#' @param components Matrix (p x k). The components where the rows correspond to the samples.
#' @param poi_types Character vector (length l). Types of POI: 'continuous' (for linear regression) or 'categorical' (for class mean).
#' @param pval Logical. If TRUE, compute p-values in addition to R-squared values. Default is TRUE.
#' @return A list containing:
#' \item{R2}{Vector (length l). Highest R-squared value between a column of `components` and each POI.}
#' \item{idxCorr}{Vector (length l). Index of the column of `components` giving the highest R-squared value.}
#' \item{allR2}{Matrix (k x l). R-squared values for each column of `components` with each POI.}
#' \item{pv}{Vector (length l). Smallest p-value association between a column of `components` and each POI (if `pval` is TRUE).}
#' \item{idxCorr2}{Vector (length l). Index of the column of `components` giving the smallest p-value (if `pval` is TRUE).}
#' \item{allpv}{Matrix (k x l). P-values for each column of `components` with each POI (if `pval` is TRUE).}
#' @export
R2 <- function(poi, components, poi_types, pval = TRUE) {
  if (is.vector(components)) {
    components <- matrix(components, ncol = 1)
  }
  if (is.vector(poi)) {
    poi <- matrix(poi, nrow = length(poi))
  }

  n_samples <- nrow(components)
  n_components <- ncol(components)
  n_poi <- length(poi_types)

  if (is.null(n_poi)) {
    stop("POI type(s) needed")
  }

  poi_rows <- nrow(poi)
  poi_cols <- ncol(poi)

  if (poi_cols != n_poi) {
    if (poi_rows == n_poi) {
      poi <- t(poi)
      warning("Transposing POI to match POI types dimension")
      poi_rows <- nrow(poi)
    } else {
      stop("POI dimensions do not match POI types dimension")
    }
  }

  if (n_samples != poi_rows) {
    if (poi_rows == n_components) {
      warning("Transposing components to match POI dimension")
      components <- t(components)
      n_components <- n_samples
      n_samples <- poi_rows
    } else {
      stop("POI and components dimensions incompatible")
    }
  }

  R2_values <- rep(-1, n_poi)
  names(R2_values) <- colnames(poi)
  idx_corr <- R2_values
  R2_tmp <- matrix(rep(-1, n_components * n_poi), n_components, n_poi, dimnames = list(colnames(components), colnames(poi)))

  if (pval) {
    p_values <- R2_values
    idx_corr2 <- R2_values
    p_values_tmp <- R2_tmp
  }

  for (component_idx in 1:n_components) {
    component_values <- components[, component_idx]
    for (poi_idx in 1:n_poi) {
      finite_indices <- is.finite(as.factor(poi[, poi_idx]))
      poi_values <- poi[finite_indices, poi_idx]
      finite_component_values <- component_values[finite_indices]

      if (poi_types[poi_idx] == "continuous") {
        estimated_values <- linear_regression_estimation(finite_component_values, poi_values)
        num_classes <- 2
      } else if (poi_types[poi_idx] == "categorical") {
        estimated_values <- class_mean_estimation(finite_component_values, poi_values)
        num_classes <- length(unique(poi_values))
      } else {
        stop("Incorrect poi_type. Select 'continuous' or 'categorical'.")
      }

      sse <- sum((finite_component_values - estimated_values)^2)
      sst <- sum((finite_component_values - mean(finite_component_values))^2)
      R2_tmp[component_idx, poi_idx] <- 1 - sse / sst

      if (pval) {
        F_value <- ((sst - sse) / (num_classes - 1)) / (sse / (n_samples - num_classes))
        p_values_tmp[component_idx, poi_idx] <- 1 - pf(F_value, num_classes - 1, n_samples - num_classes)
        if (!is.finite(p_values_tmp[component_idx, poi_idx])) {
          warning(sprintf("Non-finite p-value for component %d (pv=%g, F=%g), assigning NA", component_idx, p_values_tmp[component_idx, poi_idx], F_value))
          p_values_tmp[component_idx, poi_idx] <- NA
        }
      }
    }
  }

  for (poi_idx in 1:n_poi) {
    if (pval) {
      p_values[poi_idx] <- min(p_values_tmp[, poi_idx])
      idx_corr2[poi_idx] <- which(p_values_tmp[, poi_idx] == p_values[poi_idx])[1]
    }
    R2_values[poi_idx] <- max(R2_tmp[, poi_idx])
    idx_corr[poi_idx] <- which(R2_tmp[, poi_idx] == R2_values[poi_idx])[1]
  }

  if (pval) {
    return(list(R2 = R2_values, idxCorr = idx_corr, allR2 = R2_tmp, pv = p_values, idxCorr2 = idx_corr2, allpv = p_values_tmp))
  } else {
    return(list(R2 = R2_values, idxCorr = idx_corr, allR2 = R2_tmp))
  }
}
