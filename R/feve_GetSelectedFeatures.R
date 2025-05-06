"Functions called to set-up the data of a fEVE clustering recursion.

	2025/03/18 @yanisaspic"

feve_GetSelectedFeatures <- function(dataset, params, features_percent=0.25) {
  #' Get the variable features in an -omics dataset.
  #'
  #' The variable features have a variance superior to 0.
  #'
  #' @param dataset an -omics dataset.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param features_percent the proportion of highly variable features to select
  #' (according to the initial pool of features).
  #'
  #' @return a vector of features.
  #'
  #' @import stats
  #'
  #' @export
  #'
  n_features_init <- nrow(dataset)
  n_features <- round(n_features_init * features_percent)
  variances <- apply(X=dataset, MARGIN=1, FUN=stats::var, na.rm=TRUE)
  means <- apply(X=dataset, MARGIN=1, FUN=mean, na.rm=TRUE)
  dispersions <- variances / means
  dispersions <- dispersions[(dispersions!=0) & (is.na(dispersions)==FALSE)]
  dispersions <- dispersions[order(dispersions, decreasing=TRUE)]
  upper_bound <- min(n_features, length(dispersions))
  selected_features <- names(dispersions[1:upper_bound])
  return(selected_features)
}
