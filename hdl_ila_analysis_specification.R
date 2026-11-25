# Executable specification: longitudinal HDL-C and prevalent ILAs
# Review focus: covariate coding, functional forms, splines, and standardized CIs.
#
# USAGE
#   source("hdl_ila_analysis_specification.R")
#   input <- list(lcs = lcs, hdl = hdl)  # optionally: elig, posterior_fu, posterior_lcs
#   results <- run_hdl_ila(input, output_dir = "results")
# Or: Rscript hdl_ila_analysis_specification.R analysis_inputs.rds results
# The RDS must contain the input list above. No study data are embedded here.
# Only aggregate results/specifications are written; fitted models stay in memory.
# Optional trajectory-fitting functions are provided separately below; the main
# runner does not silently launch the computationally intensive LCMM search.
#
# PROVENANCE AND MODEL NUMBERING
# Based on the analyst's R2_code_specification.R, R2 analysis tables, and revised
# manuscript/supplement. Manuscript Model 1 = analyst Model 3; manuscript Model 2
# = analyst Model 4. The original less-adjusted analyst Models 1/2 are not exported.
# Lowest HDL-C is reported for manuscript Model 1 only.
# Historical TWM requires >=3 REMAINING measurements after removing the latest
# measurement, as specified for the revised sensitivity analysis.
# Original reported environment: R 4.0.5, lcmm 2.2.2. Numerical reproduction
# requires the original analysis inputs and trajectory starting fits.
# sessionInfo() records the executing environment.
#
# INPUT CONTRACT (canonical names; map source database fields before calling)
# lcs: one row/person, already screened for the study's non-lipid eligibility.
#   id: unique nonmissing participant identifier (character or numeric).
#   lcs_date, date_recent: Date or ISO YYYY-MM-DD; ila: 0/1.
#   age: years; sex: Men/Women (reference Men).
#   bmi_cat: 18.5-<25 (reference), <18.5, 25-<30, >=30; or supply numeric bmi.
#   smoking: Former (reference), Current; pack_years: continuous, linear.
#   drinking: None (reference), Light, Heavy; or numeric alcohol_g_week.
#     None = 0 g/week; Heavy >=210 g/week in men or >=112 in women;
#     Light = positive consumption below the sex-specific heavy threshold.
#   family_history, other_cancer, copd, low_income, statin, diabetes,
#     hypertension, ihd, stroke: 0/1; absence (0) is the reference.
#   family_history: family history of lung cancer; other_cancer: self-reported
#     history of other cancer. low_income: medical aid/lowest premium decile.
#   Comorbidities: >=2 separate claims with a principal/first secondary code
#     in the 3 years before LCS: COPD J42-J44 excluding J43.0; diabetes E10-E14;
#     hypertension I10-I13/I15; IHD I20-I25; stroke I60-I64.
#   statin: >=2 prescriptions during the 12 months before LCS.
#   region: Seoul metropolitan area (reference), Chungcheong, Honam-Jeju,
#     Yeongnam, Gangwon; or supply residence_code as a character string.
#   hdl_recent, ldl_recent, tg_recent: mg/dL from the same selected lipid panel.
#     date_recent identifies the closest qualifying panel on/before LCS, in
#     the four-calendar-year window specified in the analyst's source code.
#     The qualifying panel must have HDL-C, LDL-C, and triglycerides available.
#     Selection of that panel and claims/questionnaire preprocessing are inputs,
#     not inferred from the HDL-only longitudinal file.
#   Optional emphysema: CT-reported emphysema, 0/1/NA (NA = not assessed).
# hdl: id, date (Date/ISO), hdl (mg/dL), one adjudicated record/person/date.
#   Construction uses 2009-01-01 through LCS inclusive, HDL-C 10-150 inclusive,
#   >=3 valid observations, positive first-to-last observation duration.
# elig (optional): all participants reaching the ILA-assessment step, NOT just
#   the analytic cohort. Columns: id, assessed (0/1), screening_year, region
#   (or residence_code), age, sex, bmi_cat (or bmi), smoking, pack_years,
#   drinking (or alcohol_g_week), lung_rads, low_income, copd.
#   lung_rads is a nominal category; supply its original labels. Its reference
#   is set explicitly by the runner argument lung_rads_reference (default "1").
# posterior_fu/posterior_lcs (optional): id, class (1/2/3), prob1, prob2, prob3.
#   Labels must already correspond to the manuscript's trajectory Class 1/2/3.
#   Class 1 is the regression reference. Latent-class labels can permute on refit.
#
# DEPENDENCIES
# Base R + splines for primary models; sandwich for IPW/posterior weighting;
# lcmm for optional trajectory fitting. Install separately; no auto-installation.
# Reference documentation:
# https://stat.ethz.ch/R-manual/R-devel/library/splines/html/ns.html
# https://sandwich.r-forge.r-project.org/reference/vcovHC.html
# https://sandwich.r-forge.r-project.org/reference/vcovCL.html
# https://cecileproust-lima.github.io/lcmm/reference/lcmm.html
# https://cecileproust-lima.github.io/lcmm/reference/gridsearch.html

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop("Install the required R package: ", package, call. = FALSE)
}
require_columns <- function(d, columns, label) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop(label, " lacks: ", paste(missing, collapse = ", "), call. = FALSE)
}
check_id <- function(d, label) {
  if (anyNA(d$id) || any(!nzchar(as.character(d$id))) || anyDuplicated(d$id))
    stop(label, " must have unique, nonmissing participant IDs.")
}
check_numeric <- function(x, label, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || any(!is.na(x) & (!is.finite(x) | x < lower | x > upper)))
    stop(label, " must be numeric and within [", lower, ", ", upper, "].")
  x
}
binary <- function(x, label) {
  s <- as.character(x)
  if (any(!is.na(s) & !s %in% c("0", "1"))) stop(label, " must be coded 0/1/NA.")
  as.integer(s)
}
nominal <- function(x, levels, label) {
  s <- as.character(x)
  bad <- unique(s[!is.na(s) & !s %in% levels])
  if (length(bad)) stop(label, " has unmapped categories: ", paste(bad, collapse = ", "))
  factor(s, levels = levels, ordered = FALSE)
}
as_date <- function(x, label) {
  if (inherits(x, "Date")) return(x)
  if (!is.character(x)) stop(label, " must be Date or ISO date strings.")
  d <- as.Date(x, format = "%Y-%m-%d")
  if (any(!is.na(x) & is.na(d))) stop("Invalid dates in ", label)
  d
}
region_from_code <- function(code) {
  prefix <- substr(as.character(code), 1L, 2L)
  map <- list("Seoul metropolitan area" = c("11", "28", "41"),
              Chungcheong = c("30", "36", "43", "44"),
              "Honam-Jeju" = c("29", "45", "46", "49"),
              Yeongnam = c("26", "27", "31", "47", "48"), Gangwon = "42")
  out <- rep(NA_character_, length(prefix))
  for (area in names(map)) out[prefix %in% map[[area]]] <- area
  nominal(out, names(map), "region")
}
code_covariates <- function(data, selection = FALSE, lung_rads_reference = "1") {
  d <- as.data.frame(data, stringsAsFactors = FALSE)
  base <- c("id", "age", "sex", "smoking", "pack_years", "low_income", "copd")
  extra <- c("ila", "lcs_date", "date_recent", "hdl_recent", "ldl_recent", "tg_recent",
             "family_history", "other_cancer", "statin", "diabetes", "hypertension", "ihd", "stroke")
  require_columns(d, c(base, if (selection) c("assessed", "screening_year", "lung_rads") else extra), "data")
  check_id(d, "data"); d$id <- as.character(d$id)
  d$age <- check_numeric(d$age, "age", 0)
  d$pack_years <- check_numeric(d$pack_years, "pack_years", 0)
  d$sex <- nominal(d$sex, c("Men", "Women"), "sex")
  d$smoking <- nominal(d$smoking, c("Former", "Current"), "smoking")
  if (!"bmi_cat" %in% names(d)) {
    require_columns(d, "bmi", "data")
    d$bmi_cat <- cut(check_numeric(d$bmi, "bmi", 0), c(-Inf, 18.5, 25, 30, Inf),
                     right = FALSE, labels = c("<18.5", "18.5-<25", "25-<30", ">=30"))
  }
  d$bmi_cat <- nominal(d$bmi_cat, c("18.5-<25", "<18.5", "25-<30", ">=30"), "bmi_cat")
  if (!"drinking" %in% names(d)) {
    require_columns(d, "alcohol_g_week", "data")
    a <- check_numeric(d$alcohol_g_week, "alcohol_g_week", 0)
    d$drinking <- ifelse(a == 0, "None", ifelse(a >= ifelse(d$sex == "Men", 210, 112), "Heavy", "Light"))
  }
  d$drinking <- nominal(d$drinking, c("None", "Light", "Heavy"), "drinking")
  if (!"region" %in% names(d)) {
    require_columns(d, "residence_code", "data"); d$region <- region_from_code(d$residence_code)
  }
  d$region <- nominal(d$region, c("Seoul metropolitan area", "Chungcheong", "Honam-Jeju", "Yeongnam", "Gangwon"), "region")
  binary_names <- if (selection) c("assessed", "low_income", "copd") else
    c("ila", "family_history", "other_cancer", "copd", "low_income", "statin", "diabetes", "hypertension", "ihd", "stroke")
  for (v in binary_names) d[[v]] <- binary(d[[v]], v)
  if (selection) {
    yr <- as.character(d$screening_year)
    d$screening_year <- factor(yr, levels = sort(unique(yr[!is.na(yr)])))
    levels_lr <- sort(unique(as.character(d$lung_rads[!is.na(d$lung_rads)])))
    if (!lung_rads_reference %in% levels_lr) stop("Lung-RADS reference not present; supply its actual source label.")
    d$lung_rads <- nominal(d$lung_rads, c(lung_rads_reference, setdiff(levels_lr, lung_rads_reference)), "lung_rads")
    return(d)
  }
  for (v in c("lcs_date", "date_recent")) d[[v]] <- as_date(d[[v]], v)
  for (v in c("hdl_recent", "ldl_recent", "tg_recent")) d[[v]] <- check_numeric(d[[v]], v, 0)
  required_complete <- c("id", "sex", "ila", "lcs_date", "date_recent", "hdl_recent", "ldl_recent", "tg_recent")
  if (any(!complete.cases(d[required_complete]))) stop("Analytic input has missing outcome/sex/recent lipid-panel fields.")
  if (any(d$hdl_recent < 10 | d$hdl_recent > 150)) stop("Invalid recent HDL-C values.")
  start_recent <- as.Date(paste0(as.integer(format(d$lcs_date, "%Y")) - 4L, "-01-01"))
  if (any(d$date_recent < start_recent | d$date_recent > d$lcs_date)) stop("Recent lipid-panel date outside the specified window.")
  d$years_recent_to_lcs <- as.numeric(d$lcs_date - d$date_recent) / 365.25
  d$low_recent <- low_status(d$hdl_recent, d$sex)
  d$ldl_cat <- cut(d$ldl_recent, c(0, 100, 130, 160, 190, Inf), right = FALSE,
                   labels = c("<100", "100-<130", "130-<160", "160-<190", ">=190"))
  d$tg_cat <- cut(d$tg_recent, c(0, 150, 200, 500, Inf), right = FALSE,
                  labels = c("<150", "150-<200", "200-<500", ">=500"))
  if ("emphysema" %in% names(d)) d$emphysema <- binary(d$emphysema, "emphysema")
  d
}
low_status <- function(value, sex) {
  nominal(ifelse(value < ifelse(sex == "Men", 40, 50), "Low", "Normal"), c("Normal", "Low"), "HDL-C status")
}

# 1. Exposures: trapezoidal TWM, visit-based proportion, and nadir --------------
# TWM_i = sum_j [(HDL_ij + HDL_i,j+1)/2 * (t_i,j+1 - t_ij)] /
#         (t_i,last - t_i,first).
# Only the first-to-last observed span enters; no extension to the LCS date.
construct_exposures <- function(lcs, hdl) {
  d <- code_covariates(lcs)
  h <- as.data.frame(hdl, stringsAsFactors = FALSE)
  require_columns(h, c("id", "date", "hdl"), "hdl")
  h$id <- as.character(h$id); h$date <- as_date(h$date, "hdl$date")
  if (anyNA(h$id) || anyNA(h$date)) stop("Longitudinal HDL-C IDs/dates cannot be missing.")
  check_numeric(h$hdl, "hdl$hdl")
  j <- match(h$id, d$id)
  h$lcs_date <- d$lcs_date[j]; h$sex <- d$sex[j]
  valid <- !is.na(j) & !is.na(h$hdl) & h$hdl >= 10 & h$hdl <= 150 &
    h$date >= as.Date("2009-01-01") & h$date <= h$lcs_date
  h <- h[which(valid), , drop = FALSE]
  h <- h[order(h$id, h$date), , drop = FALSE]
  if (!nrow(h)) stop("No valid longitudinal HDL-C measurements.")
  if (anyDuplicated(h[c("id", "date")])) stop("Adjudicate duplicate participant/date HDL-C records before analysis.")
  grouped <- split(seq_len(nrow(h)), h$id)
  twm <- function(v, days) sum((head(v, -1L) + tail(v, -1L)) / 2 * diff(days)) / diff(range(days))
  summaries <- lapply(grouped, function(ix) {
    z <- h[ix, , drop = FALSE]; n <- nrow(z); days <- as.numeric(z$date)
    if (n < 3L || diff(range(days)) <= 0) return(NULL)
    minimum <- min(z$hdl); at_min <- max(which(z$hdl == minimum))
    use <- seq_len(n - 1L); historical_ok <- length(use) >= 3L && diff(range(days[use])) > 0
    data.frame(id = z$id[1L], n_meas = n, obs_years = diff(range(days)) / 365.25,
               twm = twm(z$hdl, days), prop_low = 100 * mean(z$hdl < ifelse(z$sex == "Men", 40, 50)),
               nadir = minimum, years_nadir_to_lcs = as.numeric(z$lcs_date[1L] - z$date[at_min]) / 365.25,
               n_meas_excl = if (historical_ok) n - 1L else NA_integer_,
               obs_years_excl = if (historical_ok) diff(range(days[use])) / 365.25 else NA_real_,
               twm_excl = if (historical_ok) twm(z$hdl[use], days[use]) else NA_real_, stringsAsFactors = FALSE)
  })
  summaries <- Filter(Negate(is.null), summaries)
  if (!length(summaries)) stop("No participant meets the >=3-measurement criterion.")
  s <- do.call(rbind, summaries); rownames(s) <- NULL
  keep <- d$id %in% s$id; removed <- sum(!keep); d <- d[keep, , drop = FALSE]
  # Replace any existing exposure fields rather than creating duplicate names.
  for (v in setdiff(names(s), "id")) d[[v]] <- s[[v]][match(d$id, s$id)]
  rownames(d) <- NULL
  d$low_twm <- low_status(d$twm, d$sex); d$low_twm_excl <- low_status(d$twm_excl, d$sex)
  joint_levels <- c("Normal / Normal", "Normal / Low", "Low / Normal", "Low / Low")
  joint <- function(twm, recent) nominal(ifelse(is.na(twm) | is.na(recent), NA_character_,
                                               paste(twm, recent, sep = " / ")), joint_levels, "joint category")
  d$joint <- joint(d$low_twm, d$low_recent); d$joint_excl <- joint(d$low_twm_excl, d$low_recent)
  h <- h[h$id %in% d$id, , drop = FALSE]
  list(data = d, hdl = h, excluded_for_longitudinal_eligibility = removed)
}

# 2. Final model definitions: all ordinary continuous covariates enter linearly.
# Categorical predictors use treatment coding, first declared level as reference.
# No interaction or nonlinear covariate effect is added to the primary models.
model_formula <- function(exposure = "low_twm", recent = FALSE, trajectory = FALSE,
                          historical = FALSE, emphysema = FALSE, nadir = FALSE) {
  covariates <- c("age", "sex", "bmi_cat", "smoking", "pack_years", "drinking",
                 "family_history", "other_cancer", if (emphysema) "emphysema" else "copd",
                 "diabetes", "hypertension", "ihd", "stroke", "region", "low_income", "statin",
                 "ldl_cat", "tg_cat", if (nadir) "years_nadir_to_lcs" else "years_recent_to_lcs")
  if (!trajectory) covariates <- c(covariates, if (historical) c("n_meas_excl", "obs_years_excl") else c("n_meas", "obs_years"))
  stats::reformulate(c(covariates, exposure, if (recent) "low_recent"), response = "ila")
}
fit_logistic <- function(formula, data, weight = NULL, variance = c("model", "HC0", "cluster"), cluster = "id") {
  variance <- match.arg(variance); d <- as.data.frame(data)
  mf <- stats::model.frame(formula, data = d, na.action = stats::na.pass)
  ok <- complete.cases(mf)
  if (!is.null(weight)) {
    check_numeric(d[[weight]], weight, 0)
    ok <- ok & !is.na(d[[weight]]) & d[[weight]] > 0
    if (variance == "model") stop("Specify robust covariance for weighted models.")
  }
  if (variance == "cluster") ok <- ok & !is.na(d[[cluster]])
  d <- droplevels(d[ok, , drop = FALSE]); rownames(d) <- NULL
  if (!nrow(d)) stop("No complete observations for: ", deparse(formula))
  predictors <- intersect(all.vars(stats::delete.response(stats::terms(formula))), names(d))
  contrast_list <- list()
  for (v in predictors) if (is.factor(d[[v]])) {
    if (nlevels(d[[v]]) < 2L) stop("Predictor has fewer than two observed levels: ", v)
    intended_reference <- levels(data[[v]])[1L]
    if (!identical(levels(d[[v]])[1L], intended_reference))
      stop("The prespecified reference category is absent after complete-case selection: ", v)
    contrast_list[[v]] <- stats::contr.treatment(levels(d[[v]]), base = 1L)
  }
  if (is.null(weight)) {
    fit <- stats::glm(formula, data = d, family = stats::binomial("logit"),
                      contrasts = contrast_list, na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE)
  } else {
    d$.analysis_weight <- d[[weight]]
    fit <- stats::glm(formula, data = d, family = stats::quasibinomial("logit"), weights = .analysis_weight,
                      contrasts = contrast_list, na.action = stats::na.fail, model = TRUE, x = TRUE, y = TRUE)
  }
  if (!isTRUE(fit$converged) || any(!is.finite(stats::coef(fit))))
    stop("Nonconverged or rank-deficient logistic fit; inspect the data/design matrix.")
  V <- stats::vcov(fit)
  if (variance != "model") {
    require_package("sandwich")
    V <- if (variance == "HC0") sandwich::vcovHC(fit, type = "HC0") else
      sandwich::vcovCL(fit, cluster = d[[cluster]], type = "HC0", cadjust = TRUE)
  }
  list(fit = fit, data = d, V = V, variance = variance, excluded = sum(!ok))
}

# 3. Marginal standardization and delta-method confidence intervals ----------
# For exposure value a, pbar_a = sum_i q_i expit(x_i(a)' beta), sum_i q_i = 1.
# Gradient g_a = sum_i q_i p_i(a)[1-p_i(a)] x_i(a).
# Var(pbar_a) = g_a' V g_a; Var(pbar_a-pbar_ref) = (g_a-g_ref)' V (g_a-g_ref).
# The shared-model covariance is retained in every difference CI.
# CIs are estimate +/- qnorm(.975)*SE on the probability/difference scale,
# multiplied by 100. They are not bootstrap/logit-transformed/clipped CIs.
# Standardization conditions on the empirical covariates and fitted exposures.
# Robust V for weighted models treats estimated selection weights/posterior
# probabilities as fixed; it does not propagate their first-stage uncertainty.
standardized_prevalence <- function(object, exposure, values, reference = values[1L],
                                    target = object$data, averaging_weights = NULL) {
  fit <- object$fit; beta <- stats::coef(fit); V <- object$V[names(beta), names(beta), drop = FALSE]
  d0 <- as.data.frame(target); n <- nrow(d0)
  if (!n || !exposure %in% names(d0)) stop("Invalid standardization target/exposure.")
  q <- if (is.null(averaging_weights)) rep(1 / n, n) else averaging_weights
  if (length(q) != n || any(!is.finite(q)) || any(q < 0) || sum(q) <= 0) stop("Invalid averaging weights.")
  q <- q / sum(q); tt <- stats::delete.response(stats::terms(fit))
  r <- match(reference, values); if (is.na(r)) stop("Reference must be among assigned values.")
  pg <- lapply(values, function(value) {
    d <- d0
    if (!is.null(fit$xlevels[[exposure]])) {
      if (!as.character(value) %in% fit$xlevels[[exposure]]) stop("Unfitted exposure level: ", value)
      d[[exposure]] <- factor(rep(as.character(value), n), levels = fit$xlevels[[exposure]])
    } else d[[exposure]] <- rep(as.numeric(value), n)
    mf <- stats::model.frame(tt, data = d, xlev = fit$xlevels, na.action = stats::na.fail)
    X <- stats::model.matrix(tt, mf, contrasts.arg = fit$contrasts)[, names(beta), drop = FALSE]
    p <- stats::plogis(drop(X %*% beta))
    list(p = sum(q * p), g = colSums(X * (q * p * (1 - p))))
  })
  p <- vapply(pg, function(x) x$p, numeric(1))
  G <- vapply(pg, function(x) x$g, numeric(length(beta)))
  S <- t(G) %*% V %*% G; se <- sqrt(pmax(diag(S), 0))
  GD <- sweep(G, 1L, G[, r], FUN = "-")
  se_d <- sqrt(pmax(colSums(GD * (V %*% GD)), 0)); z <- stats::qnorm(.975)
  data.frame(exposure = exposure, value = as.character(values), reference = as.character(reference),
             prevalence_percent = 100 * p, lower_percent = 100 * (p - z * se), upper_percent = 100 * (p + z * se),
             difference_pp = 100 * (p - p[r]), difference_lower_pp = 100 * (p - p[r] - z * se_d),
             difference_upper_pp = 100 * (p - p[r] + z * se_d), stringsAsFactors = FALSE)
}
format_p <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return(NA_character_)
    if (x < .001) return("<.001")
    if (x > .99) return(">.99")
    digits <- if (x < .01 || (x != .05 && round(x, 2) == .05)) 3L else 2L
    sub("^0", "", sprintf(paste0("%.", digits, "f"), x))
  }, character(1))
}
odds_ratios <- function(object) {
  b <- stats::coef(object$fit); se <- sqrt(diag(object$V)); z <- stats::qnorm(.975)
  p <- 2 * stats::pnorm(-abs(b / se))
  data.frame(term = names(b), log_OR = unname(b), SE = unname(se), OR = unname(exp(b)),
             lower = unname(exp(b - z * se)), upper = unname(exp(b + z * se)),
             P = unname(p), P_Radiology = format_p(p), row.names = NULL)
}
# Optional response-letter-only comparison; not a primary outcome/test.
wald_equal_hdl <- function(model2) {
  b <- stats::coef(model2$fit); V <- model2$V
  i <- "low_twmLow"; j <- "low_recentLow"
  if (!all(c(i, j) %in% names(b))) stop("Requires the mutually adjusted binary TWM/recent model.")
  contrast <- b[i] - b[j]; se <- sqrt(V[i, i] + V[j, j] - 2 * V[i, j])
  data.frame(log_OR_difference = unname(contrast), SE = unname(se),
             P = 2 * stats::pnorm(-abs(contrast / se)))
}

# 4. Primary analyses and prespecified sensitivity analyses ------------------
fit_primary <- function(d) {
  list(TWM_Model1 = fit_logistic(model_formula(), d),
       TWM_Model2 = fit_logistic(model_formula(recent = TRUE), d),
       Proportion_Model1 = fit_logistic(model_formula("I(prop_low / 10)"), d),
       Proportion_Model2 = fit_logistic(model_formula("I(prop_low / 10)", recent = TRUE), d),
       Lowest_Model1 = fit_logistic(model_formula("I(-nadir / 5)", nadir = TRUE), d))
}
# Scaling is performed before fitting: exact coefficients/SEs, not powers of
# rounded OR/CI values. Proportion: per 10 percentage points; nadir: per 5 mg/dL decrease.
fit_recent_continuous <- function(d, knots = c(33, 44, 54, 75)) {
  # Four published knots: full analytic-cohort percentiles .05/.35/.65/.95.
  # Set knots=NULL to recompute them (R quantile type 7) for a different dataset.
  empirical <- unname(stats::quantile(d$hdl_recent, c(.05, .35, .65, .95), type = 7, na.rm = TRUE))
  if (is.null(knots)) knots <- empirical
  if (length(knots) != 4L || any(!is.finite(knots)) || any(diff(knots) <= 0)) stop("Four distinct ordered knots are required.")
  # Natural and restricted cubic splines span the same 3-df function space here.
  # Interior knots: 35th/65th; boundaries: 5th/95th. Linear tails, no spline intercept.
  spline <- sprintf("splines::ns(hdl_recent, knots=c(%s), Boundary.knots=c(%s), intercept=FALSE)",
                    paste(format(knots[2:3], digits = 17), collapse = ","),
                    paste(format(knots[c(1, 4)], digits = 17), collapse = ","))
  f0 <- model_formula(); fl <- stats::update(f0, . ~ . + I(hdl_recent / 10))
  fs <- stats::as.formula(paste(paste(deparse(f0), collapse = " "), "+", spline))
  use <- complete.cases(stats::model.frame(fs, data = d, na.action = stats::na.pass))
  d <- d[use, , drop = FALSE]
  fits <- list(No_recent = fit_logistic(f0, d), Binary = fit_logistic(model_formula(recent = TRUE), d),
               Linear = fit_logistic(fl, d), RCS = fit_logistic(fs, d))
  lrt <- function(small, large) {
    stopifnot(identical(small$data$id, large$data$id))
    ll0 <- stats::logLik(small$fit); ll1 <- stats::logLik(large$fit)
    statistic <- 2 * as.numeric(ll1 - ll0); df <- attr(ll1, "df") - attr(ll0, "df")
    data.frame(chisq = statistic, df = df, P = stats::pchisq(statistic, df, lower.tail = FALSE))
  }
  tests <- rbind(cbind(test = "Recent HDL-C overall (RCS vs absent)", lrt(fits$No_recent, fits$RCS)),
                 cbind(test = "Recent HDL-C nonlinearity (RCS vs linear)", lrt(fits$Linear, fits$RCS)))
  list(models = fits, tests = tests, knots = knots, empirical_quantiles = empirical)
}
fit_historical <- function(d) {
  use <- !is.na(d$twm_excl); dh <- d[use, , drop = FALSE]
  list(Excluded_latest_Model2 = fit_logistic(model_formula("low_twm_excl", recent = TRUE, historical = TRUE), dh),
       All_measurements_same_participants_Model2 = fit_logistic(model_formula(recent = TRUE), dh),
       Joint_all_Model1 = fit_logistic(model_formula("joint"), d),
       Joint_excluded_latest_Model1 = fit_logistic(model_formula("joint_excl", historical = TRUE), dh))
}
fit_emphysema_sensitivity <- function(d) {
  if (!"emphysema" %in% names(d)) return(NULL)
  de <- d[!is.na(d$emphysema), , drop = FALSE]
  # Match the full covariate complete-case cohort between replacement/comparator.
  allvars <- union(all.vars(model_formula(recent = TRUE)), all.vars(model_formula(recent = TRUE, emphysema = TRUE)))
  de <- de[complete.cases(de[allvars]), , drop = FALSE]
  list(COPD_Model1 = fit_logistic(model_formula(), de), COPD_Model2 = fit_logistic(model_formula(recent = TRUE), de),
       Emphysema_Model1 = fit_logistic(model_formula(emphysema = TRUE), de),
       Emphysema_Model2 = fit_logistic(model_formula(recent = TRUE, emphysema = TRUE), de))
}
fit_selection_ipw <- function(d, elig, lung_rads_reference = "1") {
  e <- code_covariates(elig, selection = TRUE, lung_rads_reference = lung_rads_reference)
  # Screening year/region/Lung-RADS are nominal, not continuous ordinal scores.
  f <- assessed ~ screening_year + region + age + sex + bmi_cat + smoking +
    pack_years + drinking + lung_rads + low_income + copd
  selection <- fit_logistic(f, e)
  p <- stats::predict(selection$fit, type = "response")
  if (any(p <= 0 | p >= 1)) stop("Selection model has boundary probabilities; inspect separation/positivity.")
  ec <- selection$data; ec$sw <- mean(ec$assessed) / p
  assessed <- ec[ec$assessed == 1, , drop = FALSE]
  d$sw <- assessed$sw[match(d$id, assessed$id)]
  # Truncation percentiles are calculated among analytic participants with weights.
  limits <- unname(stats::quantile(d$sw, c(.01, .99), type = 7, na.rm = TRUE))
  d$sw_truncated <- pmin(pmax(d$sw, limits[1]), limits[2])
  fits <- list(Unweighted_Model1 = fit_logistic(model_formula(), d),
               Unweighted_Model2 = fit_logistic(model_formula(recent = TRUE), d),
               IPW_Model1 = fit_logistic(model_formula(), d, "sw", "HC0"),
               IPW_Model2 = fit_logistic(model_formula(recent = TRUE), d, "sw", "HC0"),
               IPW_truncated_Model1 = fit_logistic(model_formula(), d, "sw_truncated", "HC0"),
               IPW_truncated_Model2 = fit_logistic(model_formula(recent = TRUE), d, "sw_truncated", "HC0"))
  diagnostics <- function(w) {
    w <- w[!is.na(w)]
    data.frame(n = length(w), mean = mean(w), sd = stats::sd(w), min = min(w),
               p01 = unname(stats::quantile(w, .01)), median = stats::median(w),
               p99 = unname(stats::quantile(w, .99)), max = max(w), ESS = sum(w)^2 / sum(w^2))
  }
  list(models = fits, selection_model = selection,
       weights = rbind(cbind(weight = "stabilized", diagnostics(d$sw)),
                       cbind(weight = "truncated", diagnostics(d$sw_truncated))),
       selection_summary = data.frame(n_input = nrow(e), n_complete = nrow(ec),
                                      assessed_fraction = mean(ec$assessed), analytic_without_weight = sum(is.na(d$sw)),
                                      truncation_lower = limits[1], truncation_upper = limits[2]))
}

# 5. Trajectory construction and secondary association models ----------------
# Final selected specification: threshold-link LCMM of binary low HDL-C,
# 3 classes, 3-df natural cubic time spline, random intercept, common random-
# intercept variance (nwg=FALSE), intercept-only class membership model.
# Times are days, recorded at their actual visits; no time binning/interpolation.
# ns(df=3, intercept=FALSE) uses 2 interior knots at the pooled observed-time
# 1/3 and 2/3 quantiles; boundaries are the observed minimum/maximum.
# These TIME spline knots differ from the four RECENT-HDL concentration knots.
trajectory_time_data <- function(prepared, anchor = c("first", "LCS")) {
  anchor <- match.arg(anchor); h <- prepared$hdl
  ids <- sort(unique(h$id)); h$id_numeric <- match(h$id, ids)
  first <- ave(as.numeric(h$date), h$id, FUN = min)
  h$time_days <- if (anchor == "first") as.numeric(h$date) - first else as.numeric(h$date - h$lcs_date)
  h$low <- as.integer(low_status(h$hdl, h$sex) == "Low")
  list(data = h, id_map = data.frame(id_numeric = seq_along(ids), id = ids))
}
add_time_spline <- function(d) {
  basis <- splines::ns(d$time_days, df = 3, intercept = FALSE)
  if (length(unique(c(attr(basis, "knots"), attr(basis, "Boundary.knots")))) != 4L)
    stop("Repeated time knots: the requested 3-df natural spline is not identified.")
  d[c("s1", "s2", "s3")] <- as.data.frame(unclass(basis))
  list(data = d, knots = attr(basis, "knots"), boundary = attr(basis, "Boundary.knots"))
}
fit_final_trajectory <- function(prepared, anchor = c("first", "LCS"), starting_model = NULL,
                                 exploration_ids = NULL, seed = 123, maxiter = 500) {
  require_package("lcmm"); anchor <- match.arg(anchor)
  td <- trajectory_time_data(prepared, anchor); full <- add_time_spline(td$data)
  # To reproduce the original refit, pass the saved selected 3-class subsample
  # model as starting_model. If unavailable, the code below fits THAT selected
  # specification in a seeded 50% subsample before refitting the full cohort.
  # It does not repeat selection across 1-5 classes and spline df=2/3.
  # Saved original exploration IDs are preferable to recreating a random sample.
  if (is.null(starting_model)) {
    set.seed(seed)
    ids <- td$id_map$id
    if (is.null(exploration_ids)) exploration_ids <- sample(ids, floor(length(ids) / 2))
    if (!length(exploration_ids) || anyDuplicated(exploration_ids) || any(!exploration_ids %in% ids))
      stop("Invalid trajectory exploration IDs.")
    sub <- add_time_spline(td$data[td$data$id %in% exploration_ids, , drop = FALSE])$data
    one <- lcmm::lcmm(low ~ s1 + s2 + s3, random = ~1, subject = "id_numeric", link = "thresholds",
                       ng = 1, data = sub, maxiter = 500, convB = 1e-4, convL = 1e-4, convG = 1e-4)
    if (one$conv != 1L) stop("The initialization LCMM did not converge.")
    starting_model <- lcmm::gridsearch(rep = 10, maxiter = 10, minit = one,
      lcmm::lcmm(low ~ s1 + s2 + s3, mixture = ~s1 + s2 + s3, random = ~1,
                  subject = "id_numeric", link = "thresholds", ng = 3, nwg = FALSE,
                  classmb = ~1, data = sub, maxiter = 200, convB = 1e-4, convL = 1e-4, convG = 1e-4))
  }
  if (!inherits(starting_model, "lcmm") || starting_model$ng != 3L || starting_model$conv != 1L)
    stop("Supply a converged three-class starting model of the same specification.")
  fit <- lcmm::lcmm(low ~ s1 + s2 + s3, mixture = ~s1 + s2 + s3, random = ~1,
                     subject = "id_numeric", link = "thresholds", ng = 3, nwg = FALSE, classmb = ~1,
                     data = full$data, B = starting_model$best, maxiter = maxiter,
                     convB = 1e-4, convL = 1e-4, convG = 1e-4)
  if (fit$conv != 1L) stop("Full-cohort LCMM did not meet all convergence criteria; do not use its classes.")
  list(model = fit, anchor = anchor, interior_knots_days = full$knots,
       boundary_knots_days = full$boundary, id_map = td$id_map, exploration_ids = exploration_ids)
}
# Explicitly map a newly fitted model's raw class numbers to manuscript classes.
# class_order[j] = raw LCMM class corresponding to manuscript Class j.
posterior_table <- function(trajectory_fit, class_order) {
  if (!identical(sort(as.integer(class_order)), 1:3)) stop("class_order must be a permutation of 1:3.")
  pp <- trajectory_fit$model$pprob; map <- trajectory_fit$id_map
  id <- map$id[match(pp[[1L]], map$id_numeric)]
  if (anyNA(id)) stop("LCMM participant IDs cannot be mapped to the analytic cohort.")
  P <- pp[paste0("prob", class_order)]; names(P) <- paste0("prob", 1:3)
  data.frame(id = id, class = match(pp$class, class_order), P, check.names = FALSE)
}
fit_trajectory_associations <- function(d, posterior) {
  p <- as.data.frame(posterior); require_columns(p, c("id", "class", "prob1", "prob2", "prob3"), "posterior")
  p$id <- as.character(p$id); check_id(p, "posterior")
  if (!setequal(p$id, d$id)) stop("Posterior table must contain every analytic participant exactly once.")
  p <- p[match(d$id, p$id), , drop = FALSE]; P <- as.matrix(p[c("prob1", "prob2", "prob3")])
  if (!is.numeric(P) || any(!is.finite(P)) || any(P < 0 | P > 1) || any(abs(rowSums(P) - 1) > 1e-6))
    stop("Invalid posterior probabilities.")
  classes <- suppressWarnings(as.integer(as.character(p$class)))
  if (anyNA(classes) || any(!classes %in% 1:3)) stop("Invalid hard-assignment classes.")
  if (any(P[cbind(seq_len(nrow(P)), classes)] < apply(P, 1, max) - 1e-8)) stop("Hard assignment is not a maximum-posterior class.")
  hard <- d; hard$trajectory_class <- factor(classes, levels = 1:3)
  expanded <- do.call(rbind, lapply(1:3, function(k) {
    out <- d; out$trajectory_class <- factor(rep(k, nrow(d)), levels = 1:3); out$posterior_weight <- P[, k]; out
  }))
  f1 <- model_formula("trajectory_class", trajectory = TRUE)
  f2 <- model_formula("trajectory_class", trajectory = TRUE, recent = TRUE)
  # Zero-weight rows are excluded before the weighted fit AND clustered variance.
  # quasibinomial avoids noninteger-success warnings; robust inference uses V.
  fits <- list(Hard_Model1 = fit_logistic(f1, hard), Hard_Model2 = fit_logistic(f2, hard),
               Posterior_Model1 = fit_logistic(f1, expanded, "posterior_weight", "cluster"),
               Posterior_Model2 = fit_logistic(f2, expanded, "posterior_weight", "cluster"))
  prevalence <- list()
  for (nm in names(fits)) {
    obj <- fits[[nm]]
    # For assigned class contrasts, standardize over one row per participant.
    # This stays correct if some posterior probabilities are exactly zero.
    target <- obj$data[!duplicated(obj$data$id), , drop = FALSE]
    prevalence[[nm]] <- cbind(model = nm, standardized_prevalence(obj, "trajectory_class", as.character(1:3), target = target))
    if (grepl("Model2$", nm)) {
      # For recent-HDL contrasts, preserve each participant's class mixture.
      q <- if (grepl("Posterior", nm)) obj$data$posterior_weight else NULL
      prevalence[[paste0(nm, "_recent")]] <- cbind(model = nm, standardized_prevalence(obj, "low_recent", c("Normal", "Low"), averaging_weights = q))
    }
  }
  list(models = fits, prevalence = do.call(rbind, prevalence))
}

# 6. Intermeasurement intervals: pooled vs participant summaries --------------
interval_summaries <- function(prepared) {
  h <- prepared$hdl
  groups <- split(seq_len(nrow(h)), h$id)
  each <- lapply(groups, function(ix) {
    z <- h[ix, , drop = FALSE]; gap <- diff(as.numeric(sort(z$date))) / 365.25
    list(n = nrow(z), gap = gap, mean = mean(gap), median = stats::median(gap))
  })
  category <- vapply(each, function(x) if (x$n >= 7) ">=7" else as.character(x$n), character(1))
  rows <- lapply(c("3", "4", "5", "6", ">=7", "Total"), function(k) {
    chosen <- if (k == "Total") each else each[category == k]
    if (!length(chosen)) return(NULL)
    pooled <- unlist(lapply(chosen, function(x) x$gap), use.names = FALSE)
    means <- vapply(chosen, function(x) x$mean, numeric(1)); medians <- vapply(chosen, function(x) x$median, numeric(1))
    q <- function(x) unname(stats::quantile(x, c(.5, .25, .75), type = 7))
    values <- c(q(pooled), q(means), q(medians))
    names(values) <- as.vector(t(outer(c("pooled", "participant_mean", "participant_median"), c("median", "q25", "q75"), paste, sep = "_")))
    cbind(data.frame(measurements = k, participants = length(chosen), intervals = length(pooled)), as.data.frame(as.list(values)))
  })
  do.call(rbind, rows)
}

# 7. Run and export aggregate results ----------------------------------------
model_summary <- function(models) {
  do.call(rbind, lapply(names(models), function(nm) {
    o <- models[[nm]]
    data.frame(model = nm, formula = paste(deparse(stats::formula(o$fit)), collapse = " "),
               rows_used = nrow(o$data), participants = length(unique(o$data$id)),
               ILA_participants = sum(o$data$ila[!duplicated(o$data$id)]), rows_excluded = o$excluded,
               covariance = o$variance, stringsAsFactors = FALSE)
  }))
}
collect_or <- function(models) do.call(rbind, lapply(names(models), function(nm) cbind(model = nm, odds_ratios(models[[nm]]))))
run_hdl_ila <- function(input, output_dir = NULL, recent_knots = c(33, 44, 54, 75), lung_rads_reference = "1") {
  if (!is.list(input) || !all(c("lcs", "hdl") %in% names(input))) stop("Input must be a list with lcs and hdl.")
  prepared <- construct_exposures(input$lcs, input$hdl); d <- prepared$data
  primary <- fit_primary(d)
  q_prop <- unname(stats::quantile(d$prop_low, c(.25, .75), type = 7))
  q_nadir <- unname(stats::quantile(d$nadir, c(.75, .25), type = 7))
  prevalence <- lapply(names(primary), function(nm) {
    exposure <- if (grepl("TWM", nm)) "low_twm" else if (grepl("Proportion", nm)) "prop_low" else "nadir"
    values <- switch(exposure, low_twm = c("Normal", "Low"), prop_low = q_prop, nadir = q_nadir)
    cbind(model = nm, standardized_prevalence(primary[[nm]], exposure, values))
  })
  grid <- lapply(names(primary)[!grepl("TWM", names(primary))], function(nm) {
    exposure <- if (grepl("Proportion", nm)) "prop_low" else "nadir"
    values <- if (exposure == "prop_low") c(0, 25, 50, 75, 100) else
      unname(stats::quantile(d$nadir, c(.9, .75, .5, .25, .1), type = 7))
    cbind(model = nm, standardized_prevalence(primary[[nm]], exposure, values))
  })
  recent <- fit_recent_continuous(d, recent_knots)
  recent_prev <- do.call(rbind, lapply(c("Binary", "Linear", "RCS"), function(nm)
    cbind(model = paste0("Model2_", nm), standardized_prevalence(recent$models[[nm]], "low_twm", c("Normal", "Low")))))
  historical <- fit_historical(d)
  historical_prev <- rbind(
    cbind(model = "Excluded_latest_Model2", standardized_prevalence(historical$Excluded_latest_Model2, "low_twm_excl", c("Normal", "Low"))),
    cbind(model = "All_measurements_same_participants_Model2", standardized_prevalence(historical$All_measurements_same_participants_Model2, "low_twm", c("Normal", "Low"))))
  emphysema <- fit_emphysema_sensitivity(d)
  ipw <- if (!is.null(input$elig)) fit_selection_ipw(d, input$elig, lung_rads_reference) else NULL
  trajectories <- list()
  for (anchor in c("fu", "lcs")) {
    key <- paste0("posterior_", anchor)
    if (!is.null(input[[key]])) trajectories[[anchor]] <- fit_trajectory_associations(d, input[[key]])
  }
  results <- list(primary = primary, recent = recent, historical = historical, emphysema = emphysema,
                  ipw = ipw, trajectories = trajectories,
                  primary_prevalence = do.call(rbind, prevalence), primary_prevalence_grid = do.call(rbind, grid),
                  recent_prevalence = recent_prev, historical_prevalence = historical_prev,
                  intervals = interval_summaries(prepared),
                  cohort = data.frame(n_input = nrow(input$lcs), n_analytic = nrow(d), ila = sum(d$ila),
                                      excluded_longitudinal = prepared$excluded_for_longitudinal_eligibility,
                                      n_historical_eligible = sum(!is.na(d$twm_excl))))
  if (!is.null(output_dir)) {
    if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)))
      stop("Output directory is not empty; choose a new directory to avoid mixing different runs.")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    write_csv <- function(x, name) utils::write.csv(x, file.path(output_dir, paste0(name, ".csv")), row.names = FALSE, na = "")
    write_models <- function(models, prefix) {
      write_csv(collect_or(models), paste0(prefix, "_OR")); write_csv(model_summary(models), paste0(prefix, "_models"))
    }
    write_models(primary, "primary")
    for (nm in c("primary_prevalence", "primary_prevalence_grid", "recent_prevalence", "historical_prevalence", "intervals", "cohort"))
      write_csv(results[[nm]], nm)
    write_models(recent$models, "recent_continuous"); write_csv(recent$tests, "recent_HDL_likelihood_ratio_tests")
    write_csv(data.frame(percentile = c(5, 35, 65, 95), knot_used_mg_dL = recent$knots,
                         full_cohort_quantile_mg_dL = recent$empirical_quantiles), "recent_HDL_spline_knots")
    write_models(historical, "historical_and_joint")
    if (!is.null(emphysema)) write_models(emphysema, "emphysema")
    if (!is.null(ipw)) {
      write_models(ipw$models, "IPW"); write_csv(ipw$weights, "IPW_weight_diagnostics")
      write_csv(ipw$selection_summary, "IPW_selection_cohort"); write_csv(odds_ratios(ipw$selection_model), "IPW_selection_model_OR")
    }
    for (anchor in names(trajectories)) {
      write_models(trajectories[[anchor]]$models, paste0("trajectory_", anchor))
      write_csv(trajectories[[anchor]]$prevalence, paste0("trajectory_", anchor, "_prevalence"))
    }
    coding <- do.call(rbind, lapply(names(d), function(v) {
      if (!is.factor(d[[v]])) return(NULL)
      data.frame(variable = v, reference = levels(d[[v]])[1], levels = paste(levels(d[[v]]), collapse = "; "))
    }))
    write_csv(coding, "categorical_coding")
    writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
    writeLines(c("Model 1 = expanded adjustment (analyst Model 3).",
                 "Model 2 = Model 1 + binary recent HDL-C (analyst Model 4).",
                 "Lowest HDL-C: Model 1 only; OR per 5 mg/dL decrease.",
                 "Proportion: OR per 10-percentage-point increase.",
                 "All standardization CIs: two-sided 95% normal delta method; empirical covariates fixed.",
                 "Unweighted V: model-based; IPW: HC0; posterior weighting: participant-clustered HC0, cadjust=TRUE.",
                 "No bootstrap or clipping of standardized prevalence confidence intervals.",
                 paste0("IPW included: ", !is.null(ipw)), paste0("Emphysema included: ", !is.null(emphysema)),
                 paste0("Trajectory posterior inputs included: ", paste(names(trajectories), collapse = ", ")),
                 "LCMM fitting/model selection is not launched by run_hdl_ila().",
                 "Covariate preprocessing and trajectory fitting are documented in the R source."),
               file.path(output_dir, "analysis_specification.txt"))
  }
  invisible(results)
}

# Optional trajectory workflow (run deliberately on the secure analytic data):
# prepared <- construct_exposures(input$lcs, input$hdl)
# first_fit <- fit_final_trajectory(prepared, anchor="first", starting_model=saved_subsample_first)
# lcs_fit <- fit_final_trajectory(prepared, anchor="LCS", starting_model=saved_subsample_lcs)
# Inspect trajectories and map raw class numbers to the manuscript before:
# input$posterior_fu <- posterior_table(first_fit, class_order=c(1,2,3))
# input$posterior_lcs <- posterior_table(lcs_fit, class_order=c(1,2,3))
# run_hdl_ila(input, output_dir="results_with_trajectories")
# If saved starting fits are unavailable, omit starting_model to perform the
# explicitly documented subsample initialization of the selected specification.
# Export only aggregate model/curve summaries to a public repository.

# 8. Reproducible checks on synthetic data (no study records) ------------------
# Run: Rscript hdl_ila_analysis_specification.R --self-test
# Checks irregular-interval TWM, unit rescaling, complete-case alignment,
# and analytical delta gradients against central finite differences.
self_test <- function() {
  set.seed(417)
  n <- 2500L
  d <- data.frame(id = as.character(seq_len(n)), ila = stats::rbinom(n, 1, .15),
                  lcs_date = as.Date("2020-01-01"), date_recent = as.Date("2019-06-01"),
                  age = stats::runif(n, 54, 74), sex = sample(c("Men", "Women"), n, TRUE),
                  bmi = stats::runif(n, 16, 37), smoking = sample(c("Former", "Current"), n, TRUE),
                  pack_years = stats::runif(n, 20, 70), drinking = sample(c("None", "Light", "Heavy"), n, TRUE),
                  region = sample(c("Seoul metropolitan area", "Chungcheong", "Honam-Jeju", "Yeongnam", "Gangwon"), n, TRUE),
                  hdl_recent = stats::runif(n, 25, 80), ldl_recent = stats::runif(n, 50, 250), tg_recent = stats::runif(n, 50, 700))
  for (v in c("family_history", "other_cancer", "copd", "low_income", "statin", "diabetes", "hypertension", "ihd", "stroke"))
    d[[v]] <- stats::rbinom(n, 1, .25)
  h <- do.call(rbind, lapply(seq_len(n), function(i) {
    k <- sample(3:6, 1); days <- sort(sample(0:3000, k))
    data.frame(id = as.character(i), date = as.Date("2009-01-01") + days, hdl = stats::runif(k, 20, 80))
  }))
  # An independently calculable irregular-interval example: (40*365+60*730)/1095.
  h <- h[h$id != "1", ]; h <- rbind(h, data.frame(id = "1", date = as.Date("2009-01-01") + c(0,365,1095), hdl = c(30,50,70)))
  prep <- construct_exposures(d, h); a <- prep$data
  stopifnot(abs(a$twm[a$id == "1"] - 160/3) < 1e-10, is.na(a$twm_excl[a$id == "1"]))
  stopifnot(all(a$n_meas_excl[!is.na(a$n_meas_excl)] >= 3))
  # Nonconsecutive row names and omitted covariate values must not misalign standardization.
  a$age[3] <- NA_real_; rownames(a) <- paste0("participant_", seq_len(nrow(a)) * 7)
  obj <- fit_logistic(model_formula(), a)
  stopifnot(nrow(obj$data) == nrow(a)-1L, !"3" %in% obj$data$id)
  prev <- standardized_prevalence(obj, "low_twm", c("Normal", "Low"))
  nd <- obj$data; nd$low_twm <- factor("Low", levels=c("Normal","Low"))
  direct <- mean(stats::predict(obj$fit, newdata=nd, type="response"))
  stopifnot(abs(prev$prevalence_percent[2] / 100 - direct) < 1e-12,
            prev$difference_pp[1] == 0, prev$difference_lower_pp[1] == 0)
  tt <- stats::delete.response(stats::terms(obj$fit))
  X <- stats::model.matrix(tt, stats::model.frame(tt, nd, xlev=obj$fit$xlevels), contrasts.arg=obj$fit$contrasts)
  b <- stats::coef(obj$fit); X <- X[,names(b),drop=FALSE]
  p <- stats::plogis(drop(X %*% b)); analytical <- colMeans(X * (p*(1-p)))
  epsilon <- 1e-6
  numerical <- vapply(seq_along(b), function(j) {
    plus <- minus <- b; plus[j] <- plus[j]+epsilon; minus[j] <- minus[j]-epsilon
    (mean(stats::plogis(drop(X %*% plus))) - mean(stats::plogis(drop(X %*% minus))))/(2*epsilon)
  }, numeric(1))
  stopifnot(max(abs(analytical-numerical)) < 1e-6)
  raw <- fit_logistic(model_formula("nadir", nadir=TRUE), a)
  scaled <- fit_logistic(model_formula("I(-nadir / 5)", nadir=TRUE), a)
  term <- grep("nadir", names(stats::coef(scaled$fit)), value=TRUE)
  term <- term[!grepl("years", term)]
  stopifnot(abs(stats::coef(scaled$fit)[term] + 5*stats::coef(raw$fit)["nadir"]) < 1e-8)
  sens <- fit_recent_continuous(a)
  stopifnot(identical(sens$knots, c(33,44,54,75)), all(sens$tests$df == c(3,2)))
  message("Synthetic checks passed. Study-data numerical reproduction is a separate check.")
  invisible(TRUE)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (identical(args, "--self-test")) {
    self_test()
  } else if (length(args) == 2L) {
    run_hdl_ila(readRDS(args[1L]), output_dir = args[2L])
  } else {
    cat("Usage:\n  Rscript hdl_ila_analysis_specification.R analysis_inputs.rds results_dir\n",
        "  Rscript hdl_ila_analysis_specification.R --self-test\n",
        "Or source this script and call the documented functions.\n", sep = "")
  }
}
