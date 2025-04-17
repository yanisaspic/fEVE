"Functions called to set-up the data of a brEVE clustering recursion.

	2025/03/18 @yanisaspic"

breve_GetSelectedFeatures <- function(dataset, params, n_features=1000) {
  #' Get the variable features in an -omics dataset.
  #'
  #' The variable features have a variance superior to 0.
  #'
  #' @param dataset an -omics dataset.
  #' Its rows are features and its columns are samples.
  #' @param params a list of parameters (cf. `feve::get_parameters()`).
  #' @param n_features the number of highly variable features to sample.
  #'
  #' @return a vector of genes.
  #' 
  #' @import stats
  #' 
  #' @export
  #'
  variances <- apply(X=dataset, MARGIN=1, FUN=stats::var, na.rm=TRUE)
  means <- apply(X=dataset, MARGIN=1, FUN=mean, na.rm=TRUE)
  dispersions <- variances / means
  dispersions <- dispersions[(dispersions!=0) & (is.na(dispersions)==FALSE)]
  dispersions <- dispersions[order(dispersions, decreasing=TRUE)]
  upper_bound <- min(n_features, length(dispersions))
  selected_features <- names(dispersions[1:upper_bound])
  return(selected_features)
}