# EGARCH(1,1)-X volatility framework
#   * twelve cryptocurrency-stablecoin pairs (4 cryptos x 3 stablecoins)
#   * conditional variance: eGARCH(1,1), Student-t innovations
#   * external regressors entered contemporaneously in the variance eq.:
#       |Depeg|          absolute stablecoin depeg          -> delta1
#       Stress           binary market-stress dummy         -> delta2
#       |Depeg| x Stress stress interaction (key term)      -> delta3
#   * SEs from the observed information matrix (Hessian), valid for nu>2;
#     sandwich/robust SEs NOT used (would require nu>4, violated here)
#   * notation map:
#       alpha (size)      <- gamma1   (|z|-E|z| magnitude term)
#       gamma (asymmetry) <- alpha1   (z term, leverage)
#   * post-estimation: delta3 sign count, joint Wald (incl/excl BTC-DAI),
#     one-sample t-test and sign test on delta3; per-pair coefs to Excel


library(readxl)
library(rugarch)
library(dplyr)
library(openxlsx)

# 1. Load Data

btc  <- read_excel("btc_modified.xlsx")
eth  <- read_excel("eth_modified.xlsx")
xrp  <- read_excel("xrp_modified.xlsx")
bnb  <- read_excel("bnb_modified.xlsx")

usdt <- read_excel("usdt_modified.xlsx")
usdc <- read_excel("usdc_modified.xlsx")
dai  <- read_excel("dai_modified.xlsx")

stress <- read_excel("stress_periods_new.xlsx")

# 2. Build Dataset

data <- data.frame(
  btc = btc[[4]],
  eth = eth[[4]],
  xrp = xrp[[4]],
  bnb = bnb[[4]],
  usdt = abs(usdt[[3]]),
  usdc = abs(usdc[[3]]),
  dai = abs(dai[[3]]),
  S = stress$stress_daily
)

# Interaction terms
data$usdt_S <- data$usdt * data$S
data$usdc_S <- data$usdc * data$S
data$dai_S  <- data$dai  * data$S

data <- na.omit(data)

# 3. EGARCH FUNCTION (ONE COIN)

run_egarch_single <- function(ret, depeg, interaction, stress){
  
  X <- as.matrix(data[, c(depeg, stress, interaction)])
  
  spec <- ugarchspec(
    mean.model = list(
      armaOrder = c(0,0),
      include.mean = TRUE
    ),
    variance.model = list(
      model = "eGARCH",
      garchOrder = c(1,1),
      external.regressors = X
    ),
    distribution.model = "std"
  )
  
  fit <- ugarchfit(spec = spec, data = ret)
  return(fit)
}

# 4. RUN ALL MODELS

assets <- c("btc","eth","xrp","bnb")

coins <- list(
  USDT = c("usdt","S","usdt_S"),
  USDC = c("usdc","S","usdc_S"),
  DAI  = c("dai","S","dai_S")
)

results <- list()

for(asset in assets){
  for(coin in names(coins)){
    
    vars <- coins[[coin]]
    
    fit <- run_egarch_single(
      ret = data[[asset]],
      depeg = vars[1],
      interaction = vars[3],
      stress = vars[2]
    )
    
    results[[paste(asset, coin, sep = "_")]] <- fit
  }
}

# Wald test

ok <- sapply(results, function(f) f@fit$convergence == 0)
cat("converged:", sum(ok), "of", length(ok), "\n")
stopifnot(all(ok))

d3 <- t(sapply(results, function(f) f@fit$matcoef["vxreg3", c(1,2,4)]))
colnames(d3) <- c("est","se","p")
d3 <- as.data.frame(d3); d3$pair <- rownames(d3); rownames(d3) <- NULL
d3$t <- d3$est / d3$se
print(d3[, c("pair","est","se","t","p")], digits = 3)

cat("\nNegative delta3:", sum(d3$est < 0), "of", nrow(d3),
    "| positive:", d3$pair[d3$est > 0], "\n")

keep <- d3$pair != "btc_DAI"
W_all  <- sum(d3$t^2);  W_excl <- sum(d3$t[keep]^2)
cat(sprintf("Wald incl BTC-DAI: W=%.2f df=%d p=%.3g\n", W_all, nrow(d3), pchisq(W_all, nrow(d3), lower.tail=FALSE)))
cat(sprintf("Wald excl BTC-DAI: W=%.2f df=%d p=%.3g\n", W_excl, sum(keep), pchisq(W_excl, sum(keep), lower.tail=FALSE)))

tt <- t.test(d3$est[keep])
cat(sprintf("One-sample t (excl BTC-DAI): mean=%.4f t=%.2f p=%.3f\n", mean(d3$est[keep]), tt$statistic, tt$p.value))
print(binom.test(sum(d3$est[keep] < 0), sum(keep), 0.5))


# 5. EXTRACT RESULTS

extract_coefs <- function(fit){
  co <- fit@fit$matcoef          # Hessian SEs: valid for nu>2 (reviewer 3.1).
  # NOT robust.matcoef -- sandwich SEs need nu>4, violated here.
  pick <- function(rn) c(est = co[rn,1], p = co[rn,4])
  map <- list(
    "mu"               = pick("mu"),
    "omega"            = pick("omega"),
    "alpha (size)"     = pick("gamma1"),   # paper alpha = rugarch gamma1  (|z|-E|z| term)
    "beta"             = pick("beta1"),
    "gamma (asymmetry)"= pick("alpha1"),   # paper gamma = rugarch alpha1  (z term, leverage)
    "delta1 (|Depeg|)" = pick("vxreg1"),
    "delta2 (Stress)"  = pick("vxreg2"),
    "delta3 (Interact)"= pick("vxreg3"),
    "nu"               = pick("shape")
  )
  est <- sapply(map, `[`, "est"); p <- sapply(map, `[`, "p")
  stars <- ifelse(p<.01,"***", ifelse(p<.05,"**", ifelse(p<.1,"*","")))
  data.frame(Variable=names(map), Estimate=round(est,4), pvalue=round(p,4),
             Significance=stars, row.names=NULL)
}

# 6. SAVE TO EXCEL

wb <- createWorkbook()

for(name in names(results)){
  
  df <- extract_coefs(results[[name]])
  
  addWorksheet(wb, name)
  writeData(wb, name, df)
}

saveWorkbook(wb, "EGARCH_individual_models_28_06.xlsx", overwrite = TRUE)

# 7. PRINT CHECK

print("All individual models estimated:")
print(names(results))