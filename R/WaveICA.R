#' @description Wavelet decomposition.
#' @importFrom waveslim modwt
wt_decomposition <- function(data, level, wf) {
  coef <- list()
  for (k in 1:(level + 1)) {
    coef[[k]] <- matrix(NA, nrow(data), ncol(data))
  }
  cat(paste("Decomposing...\n"))
  for (j in 1:ncol(data)) {
    data_temp <- data[, j]
    x_modwt <- modwt(data_temp, wf = wf, n.levels = level)
    for (k in 1:(level + 1)) {
      coef[[k]][, j] <- x_modwt[[k]]
    }
  }
  return(coef)
}

#' @description Wavelet reconstruction.
#' @importFrom waveslim imodwt
wt_reconstruction <- function(data, data_wave_ICA, wf) {
  index <- ncol(data)
  index1 <- length(data_wave_ICA)
  data_coef <- matrix(NA, nrow(data_wave_ICA[[1]]), index1)
  data_wave <- matrix(NA, nrow(data_wave_ICA[[1]]), ncol(data_wave_ICA[[1]]))
  cat(paste("Reconstructing...\n"))
  for (i in 1:index) {
    for (j in 1:index1) {
      data_coef[, j] <- data_wave_ICA[[j]][, i]
    }
    data_temp <- data[, i]
    data_coef <- as.data.frame(data_coef)
    colnames(data_coef) <- c(paste("d", 1:(index1 - 1), sep = ""), paste("s", (index1 - 1), sep = ""))
    y <- as.list(data_coef)
    attributes(y)$class <- "modwt"
    attributes(y)$wavelet <- wf
    attributes(y)$boundary <- "periodic"
    data_wave[, i] <- imodwt(y) + mean(data_temp)
  }
  return(data_wave)
}

#' @importFrom mgcv gam
gam_wrapper <- function(x, injection_order) {
  corr <- gam(x ~ s(injection_order))
  corr_summary <- summary(corr)
  corr_r <- corr_summary$r.sq
  return(corr_r)
}

#' @title WaveICA_nonbatchwise
#' @description Removing batch effects for metabolomics data without batch information.
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @author RECETOX
#' \email{GalaxyToolsDevelopmentandDeployment@space.muni.cz}
#' @param data Dataframe. Feature table with intensities.
#' @param wf String. Wavelet function. The default is "haar".
#' @param injection_order Vector. Injection order of the samples.
#' @param alpha Float between 0 and 1. The trade-off value between the independence of samples and those
#' of variables. The default is 0.
#' @param cutoff Float between 0 and 1. The threshold of the variation explained by the injection order for
#' independent components.
#' @param K Integer. The maximal number of independent components. The default is 10.
#' @return Dataframe. Feature table with intensities corrected of batch effects.
#' @importFrom parallel mclapply
#' @export
waveica_nonbatchwise <- function(data, wf = "haar", injection_order, alpha = 0, cutoff, K = 20) {
  level <- floor(log(nrow(data), 2))
  coef <- wt_decomposition(data, level, wf)

  ### ICA
  index <- level + 1
  data_wave_ICA <- list()
  cat(paste("Performing ICA...\n"))
  for (i in (1:index)) {
    data_coef <- coef[[i]]
    data_coef_ICA <- unbiased_stICA(X = t(data_coef), k = K, alpha)
    B <- data_coef_ICA$B
    A <- data_coef_ICA$A
    B <- as.data.frame(B)

    ### Gam
    corr <- mclapply(B, gam_wrapper, injection_order)
    corr <- unlist(corr)
    label <- which(corr >= cutoff)
    B_new <- B[, label, drop = F]
    A_new <- A[, label, drop = F]
    Xn <- data_coef - t(A_new %*% t(B_new))

    data_wave_ICA[[i]] <- Xn
  }
  data_wave <- wt_reconstruction(data, data_wave_ICA, wf)

  rownames(data_wave) <- rownames(data)
  colnames(data_wave) <- colnames(data)
  return(data_wave)
}

#' @title WaveICA
#' @description Removing batch effects for metabolomics data.
#' @author Kui Deng
#' \email{dengkui_stat@163.com}
#' @author RECETOX
#' \email{GalaxyToolsDevelopmentandDeployment@space.muni.cz}
#' @param data Dataframe. Feature table with intensities.
#' @param wf String. Wavelet function, the default is "haar".
#' @param batch Vector. Batch number of each sample.
#' @param factorization String. Matrix factorization method, options are ["stICA", "SVD"]. The default is "stICA".
#' @param group Vector, optional. Type of a sample (blank, sample, QC, standard) numerically encoded to blank:0, sample:1, QC:2, standard:3.
#' @param K Integer. The maximal number of independent components (for ICA) or singular vectors (SVD). The default is 20.
#' @param batch_threshold Float between 0 and 1. The threshold to consider a component associate with the batch. The default is 0.05.
#' @param group_threshold Float between 0 and 1. The threshold to consider a component associate with the group. The default is 0.05.
#' @param alpha Float between 0 and 1. The trade-off value between the independence of samples and those
#' of variables. The default is 0.
#' @return Dataframe. Feature table with intensities corrected of batch effects.
#' @export
waveica <- function(data,
                    wf = "haar",
                    batch,
                    factorization = "stICA",
                    group = NULL,
                    K = 20,
                    batch_threshold = 0.05,
                    group_threshold = 0.05,
                    alpha = 0) {
  if (!factorization %in% c("stICA", "SVD")) {
    stop("The factorization method should be 'stICA' or 'SVD'.")
  }

  level <- floor(log(nrow(data), 2))
  coef <- wt_decomposition(data, level, wf)

  ### Factorization
  index <- level + 1
  data_wave_ICA <- list()
  cat(paste("Performing matrix factorization...\n"))
  for (i in (1:index)) {
    data_coef <- coef[[i]]
    data_coef_ICA <- normFact(factorization_method = factorization, data_matrix = t(data_coef), batch_vector = batch, batch_type = "categorical", rank = K, batch_threshold = batch_threshold, group_matrix = group, group_types = "categorical", group_threshold = group_threshold, alpha)
    data_wave_ICA[[i]] <- t(data_coef_ICA$normalized_matrix)
  }
  data_wave <- wt_reconstruction(data, data_wave_ICA, wf)

  rownames(data_wave) <- rownames(data)
  colnames(data_wave) <- colnames(data)
  return(data_wave)
}
