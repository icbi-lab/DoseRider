
#' Compare Two GAM Models
#'
#' This function compares two Generalized Additive Models (GAMs) by performing a Chi-square test on their scores.
#' The function first checks if both models are GAMs and if they are built using the same method.
#' It then computes the test statistic based on the difference of the models' Un-Biased Risk Estimators (UBREs)
#' and their degrees of freedom.
#'
#' @param model1 A gam or bam model object. The first model to compare.
#' @param model2 A gam or bam model object. The second model to compare.
#'
#' @return The p-value from the Chi-square test comparing the two models.
#'
#' @examples
#' data(mtcars)
#' fit1 <- mgcv::gam(mpg ~ s(hp), data = mtcars)
#' fit2 <- mgcv::gam(mpg ~ s(hp, k = 5), data = mtcars)
#' compareGAMM(fit1, fit2)
#'
#' @importFrom mgcv gam
#' @importFrom stats pchisq
#' @export
compareGAMM<- function(model1, model2) {
  # Check if models are gam or bam objects:
  if (!("gam" %in% class(model1)) || !("gam" %in% class(model2))) {
    stop("Models are not gam objects (i.e., built with bam()/gam()).")
  }

  # Check if models are comparable:
  if (model1$method != model2$method) {
    stop(sprintf("Models are incomparable: method model1 = %s, method model2 = %s", model1$method, model2$method))
  }

  # Extract the number of degrees of freedom (ndf) for each model:
  ndf1 <- length(model1$sp) + model1$nsdf + ifelse(length(model1$smooth) > 0,
                                                   sum(sapply(model1$smooth, FUN = function(x) {
                                                     x$null.space.dim
                                                   }, USE.NAMES = FALSE)), 0)
  ndf2 <- length(model2$sp) + model2$nsdf + ifelse(length(model2$smooth) > 0,
                                                   sum(sapply(model2$smooth, FUN = function(x) {
                                                     x$null.space.dim
                                                   }, USE.NAMES = FALSE)), 0)

  # Calculate AIC scores for the models:
  ml1 <- AIC(model1)
  ml2 <- AIC(model2)

  # Perform the chi-square test to obtain the p-value for model comparison:
  p_value <- pchisq(2 * (ml1 - ml2), abs(ndf1 - ndf2), lower.tail = FALSE)

  return(p_value)
}
