## 08_analysis_console.R ------------------------------------------------------------------
## Everything from the unit root tests on, in one script.
## Console-only version: nothing is saved to the output folder.
## Main model:      ARDL(13,2) in levels  <->  ARDL(12,1) in differences  (P = 12)
## BIC robustness:  order chosen by BIC on the model in differences (P = BIC choice)
## All coefficient tests use Newey-West (HAC) standard errors; bounds tests use OLS.

library(urca); library(aTSA); library(ARDL); library(dynlm)
library(sandwich); library(lmtest); library(car); library(strucchange)

d <- readRDS("data/clean/oil_cpi_monthly.rds")
d$doil_post <- d$doil * d$post2000                   # oil change x post-2000 dummy
dt <- ts(d[, c("lcpi", "loil", "dcpi", "doil", "doil_pos", "doil_neg",
               "loil_pos", "loil_neg", "post2000", "doil_post")],
         start = c(1982, 1), frequency = 12)
nw <- function(m) NeweyWest(m, prewhite = FALSE)     # HAC covariance, used throughout



## 1. Unit root tests ----------------------------------------------------------------
# ADF (H0: unit root), KPSS (H0: stationary), Zivot-Andrews (H0: unit root, no break)
# Levels with trend and with drift; differences with drift. ADF lags chosen by BIC (max 12).
cat("\n########## 1. UNIT ROOT TESTS ##########\n")
spec <- data.frame(var = c("lcpi", "lcpi", "loil", "loil", "dcpi", "doil"),
                   det = c("trend", "drift", "trend", "drift", "drift", "drift"))
ur <- NULL
for (i in 1:nrow(spec)) {
  y     <- na.omit(d[[spec$var[i]]])
  dates <- d$date[!is.na(d[[spec$var[i]]])]
  trend <- spec$det[i] == "trend"

  adf  <- ur.df(y, type = spec$det[i], lags = 12, selectlags = "BIC")
  k    <- sum(grepl("z.diff.lag", rownames(adf@testreg$coefficients)))   # lags chosen
  kpss <- ur.kpss(y, type = if (trend) "tau" else "mu", lags = "short")
  za   <- ur.za(y, model = if (trend) "both" else "intercept", lag = max(k, 1))

  ur <- rbind(ur, data.frame(
    Series = spec$var[i], Deterministic = spec$det[i], Lags = k,
    ADF = round(adf@teststat[1], 2),  ADF_cv5 = adf@cval[1, "5pct"],
    KPSS = round(kpss@teststat, 3),   KPSS_cv5 = kpss@cval[1, "5pct"],
    ZA = round(za@teststat, 2),       ZA_cv5 = za@cval[2],
    ZA_break = format(dates[za@bpoint], "%Y-%m")))
}
ur$ADF_rejects  <- ur$ADF  < ur$ADF_cv5     # TRUE = reject unit root
ur$KPSS_rejects <- ur$KPSS > ur$KPSS_cv5    # TRUE = reject stationarity
ur$ZA_rejects   <- ur$ZA   < ur$ZA_cv5      # TRUE = reject unit root
print(ur, row.names = FALSE)


## 2. Engle-Granger --------------------------------------------------------------------
# H0: no cointegration. Compare tau with EG critical values (Wooldridge 2020, Ch. 18):
# constant 5% = -3.34 (10% = -3.04); constant + trend 5% = -3.78.
# Ignore aTSA's p-values: they use ordinary Dickey-Fuller critical values.
cat("\n########## 2. ENGLE-GRANGER ##########\n")
coint.test(d$lcpi, d$loil, nlag = 1)

# Robustness: tau for 0-12 lags in the residual ADF regression (constant case)
u <- residuals(lm(lcpi ~ loil, data = d))
eg_lags <- data.frame(lags = 0:12, tau = NA)
for (k in 0:12) eg_lags$tau[k + 1] <- round(ur.df(u, type = "none", lags = k)@teststat[1], 2)
cat("\nEG tau by lag length (5% cv = -3.34):\n"); print(eg_lags, row.names = FALSE)


## 3-6. ARDL models: main (P = 12), then BIC robustness (P = BIC choice) -----------------------
## BIC lag order, chosen on the model in differences (max 12 lags of each)
bic <- auto_ardl(dcpi ~ doil, data = dt, max_order = 12, selection = "BIC")$best_order
cat("\nBIC-selected order (differences, max 12):\n"); print(bic)
if (bic[2] != 1) cat("Note: BIC picks", bic[2], "oil lags; the models below keep oil at lags 0-1\n")

for (P in c(12, bic[1])) {

  cat("\n\n####################################################################")
  cat("\n#####  ARDL(", P + 1, ",2) in levels  /  ARDL(", P, ",1) in differences  #####")
  cat("\n####################################################################\n")


  ## 3. Bounds tests (case III: unrestricted intercept, no trend) -------------------
  # Decide on the I(0)/I(1) bounds, not the p-value
  cat("\n---------- 3. Bounds tests ----------\n")
  m_lev <- ardl(lcpi ~ loil, data = dt, order = c(P + 1, 2))
  print(bounds_f_test(m_lev, case = 3, alpha = 0.05))
  print(bounds_t_test(m_lev, case = 3, alpha = 0.05))

  ## 4. Full-sample pass-through ----------------------------------------------------
  cat("\n---------- 4. Full-sample pass-through ----------\n")
  # Levels (ARDL-ECM): read only the d(loil) rows; levels terms are not interpreted
  ecm <- uecm(m_lev)
  cat("\nLevels model (ARDL-ECM), HAC SEs:\n");  print(coeftest(ecm, vcov = nw(ecm)))

  # Differences (main)
  dm <- ardl(dcpi ~ doil, data = dt, order = c(P, 1))
  cat("\nDifferences model, HAC SEs:\n");        print(coeftest(dm, vcov = nw(dm)))

  # Multipliers with HAC (delta-method) SEs
  cat("\nImpact multiplier, levels model:\n")
  print(multipliers(m_lev, type = "sr", vcov_matrix = nw(m_lev)))
  cat("\nImpact multiplier, differences model:\n")
  print(multipliers(dm, type = "sr", vcov_matrix = nw(dm)))
  cat("\nTotal pass-through (theta0 + theta1)/(1 - sum phi), differences model:\n")
  print(multipliers(dm, type = "lr", vcov_matrix = nw(dm)))

  # Diagnostics
  cat("\nDiagnostics (p-values)\n")
  cat("Breusch-Godfrey(12): levels", round(bgtest(ecm, order = 12)$p.value, 3),
      "| differences", round(bgtest(dm, order = 12)$p.value, 3), "\n")
  cat("Breusch-Pagan:       levels", round(bptest(ecm)$p.value, 3),
      "| differences", round(bptest(dm)$p.value, 3), "\n")

  ## 5. Pre vs post 2000 ------------------------------------------------------------
  cat("\n---------- 5. Pre vs post 2000 ----------\n")
  # Only the oil-change terms are interacted with the dummy (no jump problem in 2000)
  int_lev <- dynlm(dcpi ~ L(lcpi, 1) + L(loil, 1) + L(dcpi, 1:P) + post2000 +
                     doil + L(doil, 1) + doil_post + L(doil_post, 1), data = dt)
  int_dif <- dynlm(dcpi ~ L(dcpi, 1:P) + post2000 +
                     doil + L(doil, 1) + doil_post + L(doil_post, 1), data = dt)
  cat("\nInteraction, levels model (doil_post rows = change after 2000):\n")
  print(coeftest(int_lev, vcov = nw(int_lev)))
  cat("\nInteraction, differences model:\n")
  print(coeftest(int_dif, vcov = nw(int_dif)))

  # Wald test, H0: both doil_post coefficients = 0 (no change after 2000)
  for (mod in list(list("levels", int_lev), list("differences", int_dif))) {
    b <- names(coef(mod[[2]])); idx <- grep("doil_post", b)
    H <- diag(length(b))[idx, , drop = FALSE]             # one row per term = 0
    cat("\nWald test,", mod[[1]], "model, terms:", b[idx], "\n")
    print(linearHypothesis(mod[[2]], hypothesis.matrix = H, rhs = rep(0, length(idx)),
                           vcov. = nw(mod[[2]])))
  }

  # 5b. Chow test at 2000m1 (aside): do ALL coefficients change?
  for (mod in list(list("levels", ecm), list("differences", dm))) {
    y <- model.response(model.frame(mod[[2]])); X <- model.matrix(mod[[2]])
    D <- as.numeric(tail(d$date, nrow(X)) >= as.Date("2000-01-01"))
    cat("\nChow test (classical F),", mod[[1]], "model:\n")
    print(sctest(y ~ X - 1, type = "Chow", point = sum(D == 0)))
    XD <- X * D; k <- ncol(X)
    mc <- lm(y ~ X + XD - 1)                              # every regressor interacted with D
    cat("Chow test (HAC-robust Wald, all interactions = 0),", mod[[1]], "model:\n")
    print(linearHypothesis(mc, hypothesis.matrix = cbind(matrix(0, k, k), diag(k)),
                           rhs = rep(0, k), vcov. = nw(mc)))
  }

  ## 6. Asymmetry -------------------------------------------------------------------
  cat("\n---------- 6. Asymmetry ----------\n")
  # 6a. Differences model: rises vs falls
  asy_dif <- dynlm(dcpi ~ L(dcpi, 1:P) + doil_pos + L(doil_pos, 1) +
                     doil_neg + L(doil_neg, 1), data = dt)
  cat("\nAsymmetric differences model, HAC SEs:\n"); print(coeftest(asy_dif, vcov = nw(asy_dif)))
  cat("\nShort-run symmetry Wald, differences model (H0: symmetric):\n")
  b <- names(coef(asy_dif))
  h <- grepl("doil_pos", b) - grepl("doil_neg", b)      # sum(rise terms) - sum(fall terms) = 0
  print(linearHypothesis(asy_dif, hypothesis.matrix = matrix(h, nrow = 1), rhs = 0,
                         vcov. = nw(asy_dif)))

  # 6b. NARDL in levels (partial sums of oil rises and falls)
  nardl <- ardl(lcpi ~ loil_pos + loil_neg, data = dt, order = c(P + 1, 2, 2))
  cat("\nNARDL bounds F test:\n"); print(bounds_f_test(nardl, case = 3, alpha = 0.05))
  nardl_ecm <- uecm(nardl)
  cat("\nNARDL (levels), HAC SEs:\n"); print(coeftest(nardl_ecm, vcov = nw(nardl_ecm)))

  b  <- names(coef(nardl_ecm))
  sr <- grepl("^d\\(", b)                               # short-run (differenced) terms
  h_sr <- (sr & grepl("loil_pos", b)) - (sr & grepl("loil_neg", b))
  h_lr <- (!sr & grepl("loil_pos", b)) - (!sr & grepl("loil_neg", b))
  cat("\nShort-run symmetry Wald, NARDL (H0: symmetric). Terms:", b[h_sr != 0], "\n")
  print(linearHypothesis(nardl_ecm, hypothesis.matrix = matrix(h_sr, nrow = 1), rhs = 0,
                         vcov. = nw(nardl_ecm)))
  cat("\nLong-run symmetry Wald, NARDL (only meaningful with cointegration). Terms:", b[h_lr != 0], "\n")
  print(linearHypothesis(nardl_ecm, hypothesis.matrix = matrix(h_lr, nrow = 1), rhs = 0,
                         vcov. = nw(nardl_ecm)))
}

