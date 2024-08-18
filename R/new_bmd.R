#' #' DoseRider: Dose-Response Analysis for Omics Data
#' #'
#' #' This package provides functions to fit generalized linear mixed models (GLMMs) for dose-response analysis of omics data and calculate the Benchmark Dose (BMD) values.
#' #'
#'
#' # Load necessary libraries
#' library(lme4)
#' library(splines)
#'
#' #' Define Linear Dose-Response Function
#' #'
#' #' Defines the linear dose-response function.
#' #'
#' #' @param x Dose value.
#' #' @param b Slope parameter.
#' #' @param d Intercept parameter.
#' #'
#' #' @return Response value.
#' flin <- function(x, b, d) {
#'   return(b * x + d)
#' }
#'
#' #' Inverse Linear Dose-Response Function
#' #'
#' #' Calculates the dose value given a response value for the linear dose-response function.
#' #'
#' #' @param y Response value.
#' #' @param b Slope parameter.
#' #' @param d Intercept parameter.
#' #'
#' #' @return Dose value.
#' invlin <- function(y, b, d) {
#'   return((y - d) / b)
#' }
#'
#' #' Define Cubic Splines Dose-Response Function
#' #'
#' #' Defines the cubic splines dose-response function.
#' #'
#' #' @param x Dose value.
#' #' @param b Vector of spline coefficients.
#' #' @param d Intercept parameter.
#' #' @param knots Knot locations for the spline.
#' #'
#' #' @return Response value.
#' fspline <- function(x, b, d, knots, df) {
#'   spline_basis <- bs(x, degree = 3, df = df, Boundary.knots = knots, intercept = TRUE)
#'   return(d + spline_basis %*% b)
#' }
#'
#'
#' #' Inverse Cubic Splines Dose-Response Function
#' #'
#' #' Calculates the dose value given a response value for the cubic splines dose-response function.
#' #'
#' #' @param y Response value.
#' #' @param b Vector of spline coefficients.
#' #' @param d Intercept parameter.
#' #' @param knots Knot locations for the spline.
#' #' @param interval Interval for the numerical solution.
#' #'
#' #' @return Dose value.
#' invSpline <- function(y, b, d, knots, interval = c(0, 100)) {
#'   spline_fun <- function(x) fspline(x, b, d, knots) - y
#'   root <- uniroot(spline_fun, interval = interval)$root
#'   return(root)
#' }
#'
#' #' Calculate BMD for GLMM
#' #'
#' #' Calculates the Benchmark Dose (BMD) values using parameters from a fitted GLMM model.
#' #'
#' #' @param model Fitted GLMM model.
#' #' @param z Number of standard deviations for BMR.
#' #' @param x Percent change for BMR.
#' #' @param minBMD Minimum BMD value.
#' #' @param ratio2switchinlog Ratio to switch to logarithmic space.
#' #'
#' #' @return List containing the BMD calculations.
#' #' @export
#' bmdcalc_glmer <- function(model, z = 1, x = 10, minBMD, ratio2switchinlog = 100, dose_col = "Dose", long_df = long_df) {
#'   fixed_effects <- fixef(model)
#'   random_effects <- ranef(model)$gene
#'   model_formula <- as.character(formula(model))
#'
#'   data <- long_df
#'   dose <- as.vector(data$Dose)
#'   dosemax <- max(dose)
#'   dosemin <- min(dose[dose != 0])
#'
#'   if (missing(minBMD)) { minBMD <- dosemin / 100 }
#'   if (minBMD > dosemin) { minBMD <- dosemin / 100 }
#'   if (minBMD < 0) { stop("minBMD should be a strictly positive value.") }
#'
#'   gene_names <- rownames(random_effects)
#'   dcalc <- data.frame(xextrem = rep(dosemax, length(gene_names)),
#'                       yextrem = NA,
#'                       y0 = NA,
#'                       ydosemax = NA,
#'                       yp = NA,
#'                       ysd = NA,
#'                       BMDp = NA,
#'                       BMDsd = NA,
#'                       row.names = gene_names)
#'
#'   xdiv100 <- x / 100
#'
#'   # Extract residuals and compute standard deviations per gene
#'   residuals <- residuals(model)
#'   gene_sd <- tapply(residuals, data$gene, sd)
#'
#'   for (i in 1:nrow(random_effects)) {
#'     gene <- rownames(random_effects)[i]
#'     y0 <- random_effects[gene, "(Intercept)"] + fixed_effects["(Intercept)"]
#'
#'     if (grepl("Dose", model_formula)) {
#'       b <- fixed_effects["Dose"]
#'       d <- fixed_effects["(Intercept)"]
#'
#'       ydosemax <- dcalc$ydosemax[i] <- flin(dosemax, b, d)
#'       dcalc[gene,"y0"] <- y0
#'       dcalc[gene,"yp"] <- y0 * (1 + xdiv100)
#'       #dcalc[gene,"BMDp"] <- invlin(dcalc$yp[i], b, d)
#'       dcalc[gene,"ysd"] <- y0 + z * gene_sd[gene]
#'       dcalc[gene,"BMDsd"] <- invlin(dcalc$ysd[i], b, d)
#'
#'     } else if (grepl("bs\\(Dose\\)", model_formula)) {
#'
#'       knots <- attr(bs(data$Dose, degree = 3), "Boundary.knots")
#'       b <- fixed_effects[grepl("bs\\(Dose\\)", names(fixed_effects))]
#'       d <- fixed_effects["(Intercept)"]
#'
#'       ydosemax <- dcalc$ydosemax[i] <- fspline(dosemax, b, d, knots, df = 2)
#'       dcalc[gene,"y0"] <- y0
#'       dcalc[gene,"yp"] <- y0 * (1 + xdiv100)
#'       dcalc[gene,"BMDp"] <- invSpline(dcalc$yp[i], b, d, knots, interval = c(minBMD, dosemax))
#'       dcalc[gene,"ysd"] <- y0 + z * gene_sd[gene]
#'       dcalc[gene,"BMDsd"] <- invSpline(dcalc$ysd[i], b, d, knots, interval = c(minBMD, dosemax))
#'     }
#'
#'     dcalc$BMDsd[i] <- max(dcalc$BMDsd[i], minBMD)
#'     dcalc$BMDp[i] <- max(dcalc$BMDp[i], minBMD)
#'   }
#'
#'   dcalc$BMDp[dcalc$BMDp > dosemax] <- NA
#'   dcalc$BMDsd[dcalc$BMDsd > dosemax] <- NA
#'
#'   reslist <- list(res = dcalc, z = z, x = x, minBMD = minBMD, ratio2switchinlog = ratio2switchinlog)
#'
#'   return(structure(reslist, class = "bmdcalc"))
#' }
#'
#' # Assuming `predicted_values` is a data frame with columns `gene`, `Dose`, and `predictions`
#' # Example:
#' # predicted_values <- data.frame(gene = ..., Dose = ..., predictions = ...)
#' predicted_values <- dose_rider_results$BENPORATH_ES_CORE_NINE_CORRELATED$Smooth_Predictions[[1]]
#'
#' # Initialize dcalc with row names set to gene names
#' gene_names <- unique(predicted_values$gene)
#' dcalc <- data.frame(xextrem = rep(dosemax, length(gene_names)),
#'                     yextrem = NA,
#'                     y0 = NA,
#'                     ydosemax = NA,
#'                     yp = NA,
#'                     ysd = NA,
#'                     BMDp = NA,
#'                     BMDp_low = NA,
#'                     BMDp_high = NA,
#'                     BMDsd = NA,
#'                     BMDsd_low = NA,
#'                     BMDsd_high = NA,
#'                     row.names = gene_names)
#'
#' # Calculate y0 and standard deviations for each gene
#' for (gene in gene_names) {
#'   gene_data <- predicted_values[predicted_values$gene == gene, ]
#'   y0 <- mean(gene_data$predictions[gene_data$Dose == 0])  # Baseline response (Dose = 0)
#'   gene_sd <- sd(gene_data$predictions[gene_data$Dose == 0])  # Baseline response (Dose = 0)
#'
#'   dcalc[gene, "y0"] <- y0
#'   dcalc[gene, "yp"] <- y0 * (1 + xdiv100)
#'   dcalc[gene, "ysd"] <- y0 + z * gene_sd
#'
#'   # Calculate BMDp by finding the dose that corresponds to yp
#'   dose_yp <- approx(gene_data$predictions, gene_data$Dose, xout = dcalc[gene, "yp"])$y
#'   dcalc[gene, "BMDp"] <- dose_yp
#'
#'   # Calculate BMDsd by finding the dose that corresponds to ysd
#'   dose_ysd <- approx(gene_data$predictions, gene_data$Dose, xout = dcalc[gene, "ysd"])$y
#'   dcalc[gene, "BMDsd"] <- dose_ysd
#' }
#'
#' # Number of bootstrap samples
#' n_boot <- 1000
#'
#' # Initialize matrix to store bootstrap results
#' bmd_boot <- matrix(NA, nrow = n_boot, ncol = length(gene_names) * 2)
#'
#' # Define function to compute BMD for a bootstrap sample
#' compute_bmd <- function(data, indices, xdiv100, z, gene_sd) {
#'   resampled_data <- data[indices, ]
#'   bmd
#'
