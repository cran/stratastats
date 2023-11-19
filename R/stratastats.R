#' Stratified Analysis of 2x2 Contingency Tables
#'
#' This function performs a comprehensive stratified analysis of 2x2 contingency tables. It calculates odds ratios (OR),
#' 95 percent confidence intervals (CI), and conducts chi-squared tests, Cochran-Mantel-Haenszel (CMH) tests, and
#' Breslow-Day-Tarone tests for homogeneity of the odds ratios across strata. Additionally, it produces a
#' nicely-formatted table using the \code{\link[gt]{gt}} package, which includes both the results and their
#' interpretation, facilitating easier understanding and presentation of the analysis. The function is designed
#' to work with either a list of 2x2 tables or a 3-dimensional array representing stratified tables.\cr
#' Visit this \href{https://drive.google.com/file/d/1ofS6COi6RFzNm-0K5-uZqeRQN4gu0T-i/view?usp=sharing}{LINK} to access
#' the package's vignette.\cr
#'
#' @param input.data A list of 2x2 contingency tables or a 3-dimensional array where each 2x2 slice along the third
#'                   dimension represents a stratum.
#' @param variable.names Optional. A character vector containing the names of the two variables being cross-tabulated.
#'                      If not provided for a 3D array, they are extracted from the array's dimnames.
#'                      For a list of tables, these must be provided.
#' @param stratifying.variable.name Optional. A character vector containing the name of the stratifying variable.
#'                                  If not provided for a 3D array, it is extracted from the array's dimnames.
#' @param flip.or Optional. Logical value indicating whether to flip the odds ratio (default is FALSE).
#'                If TRUE, the reciprocal of the OR is calculated.
#' @param table.font.size Optional. Font size for the tables in the output (default is 13).
#' @param source.note.font.size Optional. Font size for the source notes in the output (default is 11).
#'
#' @details
#' The function employs statistical techniques appropriate for stratified data analysis. The odds ratio for each stratum
#' and the combined strata (marginal table) are calculated. Confidence intervals are derived based on the standard
#' error of the log odds ratio. Chi-squared tests are conducted for each table to assess the association between the
#' variables at each level of the stratifying variable. The CMH test assesses the overall association while accounting for
#' stratification. The Breslow-Day-Tarone test, implemented from Michael Hoehle's function (see References), evaluates
#' the homogeneity of odds ratios across strata.\cr
#'
#' The output includes a detailed breakdown of the odds ratios, confidence intervals, and test statistics for each
#' stratum and the marginal table. It also presents combined results with annotations explaining the significance and
#' implications of the tests.\cr
#'
#' Interpretational Scenarios:\cr
#'
#' 1. Significant CMH test with homogeneity (non-significant Breslow-Day-Tarone test): Indicates conditional dependence and
#' consistent association across strata. The common odds ratio (from CMH test) is a reliable summary of the association.\cr
#'
#' 2. Significant CMH test with heterogeneity (significant Breslow-Day-Tarone test): Suggests conditional dependence and
#' varying strength or direction of association across strata (interaction), cautioning against a simple summary
#' of the association.\cr
#'
#' 3. Non-significant CMH test: Implies conditional independence.\cr
#'
#' Note that the interpretation guidelines provided by the function are suggested based on the statistical tests' outcomes
#' and should be further evaluated within the context of your study.
#'
#' @return
#' A list containing the following elements:
#' \itemize{
#'   \item \code{odds_ratios}: A data frame of odds ratios and their 95 percent confidence intervals for each stratum and the marginal table.
#'    \item \code{partial_chi_sq_results}: A character vector storing the result of the chi-squared test for each partial table.
#'    \item \code{marginal_chi_sq_result}: A character vector storing the result of the chi-squared test for the marginal table.
#'    \item \code{cmh_test_result}: Result of the Cochran-Mantel-Haenszel test, providing an overall measure of association
#'     while accounting for stratification.
#'    \item \code{breslow_day_test_result}: Result of the Breslow-Day-Tarone test, used to assess the homogeneity of the odds ratios
#'     across different strata.
#'   }
#'
#'
#' @seealso
#' Refer to \code{\link{chisq.test}} for chi-squared tests, and to
#' \code{\link[stats]{mantelhaen.test}} for the Cochran-Mantel-Haenszel test.
#'
#'
#' @references Azen, R., & Walker, C. M. (2021). Categorical data analysis for the behavioral and social sciences (2nd ed.).
#' New York: Routledge.
#'
#' @references Breslow, N. E., & Day, N. E. (1980). Statistical methods in cancer research. Volume I -
#' The analysis of case-control studies. IARC Scientific Publications.
#'
#' @references Hoehle, M. (2000). Breslow-Day-Tarone Test.
#' Retrieved from https://online.stat.psu.edu/onlinecourses/sites/stat504/files/lesson04/breslowday.test_.R
#'
#' @references Lachin, J. M. (2000). Biostatistical methods: The assessment of relative risks. Wiley.
#'
#'
#' @import gt
#' @importFrom stats chisq.test mantelhaen.test pchisq
#' @importFrom abind abind
#'
#' @keywords stratastats
#'
#' @export
#'
#' @examples
#'
#' # EXAMPLE 1
#' # Survival on the Titanic
#'
#' # create three individual partial tables
#'
#' table1 <- matrix(c(118,5,61,139), byrow = TRUE, ncol=2)
#' table2 <- matrix(c(146,12,25,94), byrow = TRUE, ncol=2)
#' table3 <- matrix(c(418,110,75,106), byrow = TRUE, ncol=2)
#'
#' # make a list
#'
#' tables <- list(table1, table2, table3)
#'
#' # specify the variable names
#'
#' varnames <- c("Survival", "Gender")
#' stratvar <- "Class"
#'
#' # carry out the analysis
#' results <- stratastats(input.data = tables, variable.names = varnames,
#' stratifying.variable.name = stratvar)
#'
#'
#' # EXAMPLE 2
#' # Smoking status and breathing test results (after Azen-Walker 2021)
#'
#'
#' table1 <- matrix(c(577, 34, 682, 57), byrow = TRUE, ncol=2)
#' table2 <- matrix(c(164,4,245,74), byrow = TRUE, ncol=2)
#'
#' tables <- list(table1, table2)
#'
#' varnames <- c("Smoking Status", "Breathing Test Result")
#' stratvar <- "Age"
#'
#' results <-stratastats(input.data = tables, variable.names = varnames,
#' stratifying.variable.name = stratvar)
#'
#' # EXAMPLE 3
#' # Admission to graduate school at Berkeley in 1973 (3-dimensional array).
#' # Since the array contains variables name and the name of the stratifying variable,
#' # all we need to feed into the function is the dataset name 'UCBAdmissions'. However,
#' # the name of the variables can be customised using either the 'variable.names'
#' # or 'stratifying.variable.name' # parameter, or both.
#'
#' results <- stratastats(input.data = UCBAdmissions)
#'
#'
stratastats <- function(input.data, variable.names = NULL, stratifying.variable.name = NULL, flip.or = FALSE, table.font.size = 13, source.note.font.size = 11) {

  OR = CI_lower = CI_upper = NULL

  if (is.array(input.data) && length(dim(input.data)) == 3) {
    # Convert 3D array to a list of 2x2 tables
    list.of.tables <- lapply(1:dim(input.data)[3], function(i) input.data[,,i])

    # Use provided names or extract from array
    if (is.null(variable.names)) {
      variable.names <- names(dimnames(input.data))[1:2]
    }
    if (is.null(stratifying.variable.name)) {
      stratifying.variable.name <- names(dimnames(input.data))[3]
    }

    # Check if names are correctly specified
    if (length(variable.names) != 2 || is.null(stratifying.variable.name)) {
      stop("Could not correctly specify variable names and stratifying variable name.")
    }
  } else if (is.list(input.data)) {
    list.of.tables <- input.data

    # Require variable names and stratifying variable name for a list input
    if (is.null(variable.names) || is.null(stratifying.variable.name)) {
      stop("For a list of tables, please provide variable names and stratifying variable name.")
    }
  } else {
    stop("Input data must be either a list of 2x2 tables or a 3D array.")
  }

  # Check if the list is not empty and if variable names are provided
  if(length(list.of.tables) < 1) {
    stop("The list of tables is empty.")
  }
  if(length(variable.names) != 2) {
    stop("Please provide the names of the two variables being cross-tabulated.")
  }

  # Check if all tables are 2x2
  if(any(sapply(list.of.tables, function(tbl) any(dim(tbl) != c(2, 2))))) {
    stop("All tables must be 2x2.")
  }

  # Function to calculate OR and 95% CI
  calculate_or_and_ci <- function(tbl, flip.or) {
    a <- tbl[1, 1]
    b <- tbl[1, 2]
    c <- tbl[2, 1]
    d <- tbl[2, 2]

    or <- (a * d) / (b * c)

    if (flip.or) {
      or <- 1 / or
    }

    se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
    ci_lower <- exp(log(or) - 1.96 * se_log_or)
    ci_upper <- exp(log(or) + 1.96 * se_log_or)

    return(c(OR = or, CI_lower = ci_lower, CI_upper = ci_upper))
  }

  # Sum the tables to get the marginal table
  marginal_table <- Reduce("+", list.of.tables)

  # Calculate Odds Ratio and 95% CI for the marginal table
  marginal_or_and_ci <- calculate_or_and_ci(marginal_table, flip.or)

  # Calculate Odds Ratio and 95% CI for each partial table
  ors_and_cis <- lapply(list.of.tables, calculate_or_and_ci, flip.or)

  # Perform chi-squared tests for each partial table and create annotations with table numbers
  chi_sq_annotations <- sapply(seq_along(list.of.tables), function(i) {
    tbl <- list.of.tables[[i]]
    test_result <- chisq.test(tbl, simulate.p.value = FALSE, correct = FALSE)
    paste(sprintf("Table %d chi-sq:", i),
          format(test_result$statistic, digits = 2),
          "df:", test_result$parameter,
          "p-value:", signif(test_result$p.value, digits = 3))
  }, USE.NAMES = FALSE)

  # Combine all chi-squared annotations into one string with an initial newline and line breaks for readability
  chi_sq_annotation <- paste0("\n", paste(sapply(seq_along(list.of.tables), function(i) {
    tbl <- list.of.tables[[i]]
    test_result <- chisq.test(tbl, simulate.p.value = FALSE, correct = FALSE)
    sprintf("Partial Table %d chi-sq: %.2f, df: %d, p-value: %.3f", i,
            test_result$statistic, test_result$parameter,
            signif(test_result$p.value, digits = 3))
  }), collapse = "\n\n"))

  # Perform chi-squared test on the marginal table
  marginal_chi_sq_test <- chisq.test(marginal_table, simulate.p.value = FALSE, correct = FALSE)
  marginal_chi_sq_annotation <- sprintf("Marginal Table chi-sq: %.2f, df: %d, p-value: %.3f",
                                        marginal_chi_sq_test$statistic,
                                        marginal_chi_sq_test$parameter,
                                        signif(marginal_chi_sq_test$p.value, digits = 3))

  # Combine partial tables and marginal table chi-squared annotations
  chi_sq_annotation <- paste0(chi_sq_annotation, "\n\n", marginal_chi_sq_annotation)

  # Convert the list of tables to a 3D array
  array_of_tables <- abind::abind(list.of.tables, along = 3)

  # Cochran-Mantel-Haenszel test for conditional independence
  cmh_test <- mantelhaen.test(array_of_tables, correct = FALSE)

  breslowday.test <- function(x) {
    ######################################################################
    # Function to perform the Breslow and Day (1980) test including
    # the corrected test by Tarone
    # Uses the equations in Lachin (2000) p. 124-125.
    #
    # Programmed by Michael Hoehle <http://www-m4.ma.tum.de/pers/hoehle>
    # Note that the results of the Tarone corrected test do
    # not correspond to the numbers in the Lachin book...
    #
    # Params:
    #  x - a 2x2xK contingency table
    #
    # Returns:
    #  a vector with three values
    #   1st value is the Breslow and Day test statistic
    #   2nd value is the correct test by Tarone
    #   3rd value - p value based on the Tarone test statistic
    #               using a \chi^2(K-1) distribution
    ######################################################################

    #Find the common OR based on Mantel-Haenszel
    or.hat.mh <- mantelhaen.test(x)$estimate
    #Number of strata
    K <- dim(x)[3]
    #Value of the Statistic
    X2.HBD <- 0
    #Value of aj, tildeaj and Var.aj
    a <- tildea <- Var.a <- numeric(K)

    for (j in 1:K) {
      #Find marginals of table j
      mj <- apply(x[,,j], MARGIN=1, sum)
      nj <- apply(x[,,j], MARGIN=2, sum)

      #Solve for tilde(a)_j
      coef <- c(-mj[1]*nj[1] * or.hat.mh, nj[2]-mj[1]+or.hat.mh*(nj[1]+mj[1]),
                1-or.hat.mh)
      sols <- Re(polyroot(coef))
      #Take the root, which fulfills 0 < tilde(a)_j <= min(n1_j, m1_j)
      tildeaj <- sols[(0 < sols) &  (sols <= min(nj[1],mj[1]))]
      #Observed value
      aj <- x[1,1,j]

      #Determine other expected cell entries
      tildebj <- mj[1] - tildeaj
      tildecj <- nj[1] - tildeaj
      tildedj <- mj[2] - tildecj

      #Compute \hat{\Var}(a_j | \widehat{\OR}_MH)
      Var.aj <- (1/tildeaj + 1/tildebj + 1/tildecj + 1/tildedj)^(-1)

      #Compute contribution
      X2.HBD <- X2.HBD + as.numeric((aj - tildeaj)^2 / Var.aj)

      #Assign found value for later computations
      a[j] <- aj ;  tildea[j] <- tildeaj ; Var.a[j] <- Var.aj
    }

    #Compute Tarone corrected test
    X2.HBDT <-as.numeric( X2.HBD -  (sum(a) - sum(tildea))^2/sum(Var.aj) )

    #Compute p-value based on the Tarone corrected test
    p <- 1-pchisq(X2.HBDT, df=K-1)

    res <- list(X2.HBD=X2.HBD,X2.HBDT=X2.HBDT,p=p)
    class(res) <- "bdtest"
    return(res)
  }


  # Breslow-Day test for homogeneity
  bd_test_results <- breslowday.test(array_of_tables)

  # Annotations for conditional independence and homogeneity with key statistics
  cmh_annotation <- sprintf(
    "\n\n(B) The Cochran-Mantel-Haenszel test is %ssignificant (chi-sq: %.2f; df: %d; p-value: %.3f), suggesting %s.",
    ifelse(cmh_test$p.value < 0.05, "", "not "),
    cmh_test$statistic,
    cmh_test$parameter,
    cmh_test$p.value,
    ifelse(cmh_test$p.value < 0.05, "conditional dependence (the odds ratio in at least one of the partial tables is not equal to 1)", "conditional independence (the odds ratios in all of the partial tables are equal to 1)")
  )

  bd_annotation <- sprintf(
    "\n\n(C) The Breslow-Day-Tarone test is %ssignificant (chi-sq: %.2f; df: %d; p-value: %.3f), indicating %s of the odds ratios across strata.",
    ifelse(bd_test_results$p < 0.05, "", "not "),
    bd_test_results$X2.HBDT,
    dim(array_of_tables)[3] - 1,  # degrees of freedom for Breslow-Day test
    bd_test_results$p,
    ifelse(bd_test_results$p < 0.05, "heterogeneity", "homogeneity")
  )

  # Set an empty slot for interpretation
  interpretation <- ""

  # Add chi-squared test annotations to the interpretation string
  interpretation <- paste(interpretation, "\n\n(A) Chi-squared test results:\n", chi_sq_annotation)

  # Add the combineed annotations for CHM and BD test
  interpretation <- paste(interpretation, cmh_annotation, bd_annotation)

  # Extract common OR and CI from cmh_test
  common_or <- round(cmh_test$estimate,3)
  common_ci <- round(cmh_test$conf.int,3)

  # Flip OR and CI if flip.or is TRUE
  if (flip.or) {
    common_or <- 1 / common_or
    common_ci <- 1 / rev(common_ci)  # Invert and reverse the CI bounds
  }

  # Round the values for presentation
  common_or <- round(common_or, 3)
  common_ci <- round(common_ci, 3)

  # Extract Odds Ratio for the marginal table
  marginal_or <- round(marginal_or_and_ci["OR"],3)

  # Function to detect if there is a change in direction of ORs across strata
  detect_direction_change <- function(ors_and_cis) {
    # Get OR values with respect to the flip.or setting
    or_values <- sapply(ors_and_cis, function(orc) orc["OR"])

    # Check if any OR is less than 1 when not flipped, and greater than 1 when flipped
    any(or_values < 1) && any(or_values > 1)
  }

  # Interpretation based on homogeneity or heterogeneity of ORs
  if(bd_test_results$p < 0.05) {
    # Heterogeneity exists
    if (detect_direction_change(ors_and_cis)) {
      # Interaction effect due to change in direction
      interpretation <- paste0(interpretation, "\n\n(D) Significant heterogeneity of odds ratios across strata and a change in the direction of the association between '", variable.names[1], "' and '", variable.names[2], "' have been detected, indicating an interaction effect. The stratifying variable '", stratifying.variable.name, "' modifies the direction of the relationship between '", variable.names[1], "' and '", variable.names[2], "'. Note that, in this case, since not all the conditional odds ratios are in the same direction, the result of the CMH test is to be interpreted with caution.")
    } else {
      # Specification effect
      interpretation <- paste0(interpretation, "\n\n(D) Significant heterogeneity of odds ratios across strata and a change in the strength of the association between '", variable.names[1], "' and '", variable.names[2], "' have been detected, suggesting an interaction effect. The stratifying variable '", stratifying.variable.name, "' modifies the strength of the relationship between '", variable.names[1], "' and '", variable.names[2], "', but does not change the direction of the association.")
    }
  } else {
    # Homogeneity exists
    interpretation <- paste0(interpretation, "\n\n(D) Given the homogeneity of odds ratios across strata, '", stratifying.variable.name, "' does not significantly modify the association between '", variable.names[1], "' and '", variable.names[2], "'. This means that the conditional association between '", variable.names[1], "' and '", variable.names[2], "' is the same (in direction and magnitude) at each level of the stratifying variable '", stratifying.variable.name, "'. The association, which does not significantly differ across the levels of '", stratifying.variable.name, "', can be summarised using the Mantel-Haenszel estimate of a common odds ratio (", common_or, ").")
  }


  # Construct the results data frame
  results_df <- data.frame(
    Stratum = c(paste("Partial Table", seq_along(list.of.tables)), "Marginal Table", "MH common OR"),
    OR = c(sapply(ors_and_cis, `[[`, "OR"), marginal_or_and_ci["OR"], common_or),
    CI_lower = c(sapply(ors_and_cis, `[[`, "CI_lower"), marginal_or_and_ci["CI_lower"], common_ci[1]),
    CI_upper = c(sapply(ors_and_cis, `[[`, "CI_upper"), marginal_or_and_ci["CI_upper"], common_ci[2])
  )

  # Create the gt table
  gt_table <- gt(results_df) %>%
    tab_header(
      title = md("**Stratified Analysis Results**")
    ) %>%
    cols_label(
      Stratum = md("**Table**"),
      OR = md("**Odds Ratio**"),
      CI_lower = md("**95% CI Lower**"),
      CI_upper = md("**95% CI Upper**")
    ) %>%
    fmt_number(
      columns = c(OR, CI_lower, CI_upper),
      decimals = 3
    ) %>%
    tab_options(
      table.font.size = px(table.font.size),
      source_notes.font.size = px(source.note.font.size)
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_title(groups = "title")
    ) %>%
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_column_labels()
    ) %>%
    tab_source_note(
      source_note = md(paste(interpretation, sep = "\n"))
    )

  # Print the gt table
  print(gt_table)

  return(list(
    odds_ratios = results_df,
    partial_chi_sq_results = chi_sq_annotations,
    marginal_chi_sq_result = marginal_chi_sq_annotation,
    cmh_test_result = cmh_test,
    breslow_day_test_result = bd_test_results
  ))
}
