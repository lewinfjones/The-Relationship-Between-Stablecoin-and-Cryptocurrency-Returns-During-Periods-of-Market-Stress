# =====================================================================
# Granger-causality / predictive-regression framework
#   * lag length = BIC (Schwarz, "SC(n)") per pair, applied to BOTH directions
#   * full sample, stress as contemporaneous dummy + interaction, HAC inference
#   * three F-tests:
#       Granger  : all gamma = all phi = 0   
#       baseline : all gamma = 0             
#       stress   : all phi   = 0            
#   * reports F-statistic AND p-value for each; Ljung-Box
# =====================================================================
library(readxl)
library(vars)
library(sandwich)
library(lmtest)
library(car)



BTC <- read_xlsx("btc_modified.xlsx")[[3]]
ETH <- read_xlsx("eth_modified.xlsx")[[3]]
BNB <- read_xlsx("bnb_modified.xlsx")[[3]]
XRP <- read_xlsx("xrp_modified.xlsx")[[3]]
USDT<- read_xlsx("usdt_modified.xlsx")[[3]]
USDC<- read_xlsx("usdc_modified.xlsx")[[3]]
DAI <- read_xlsx("dai_modified.xlsx")[[3]]
stress_all <- read_xlsx("stress_periods_new.xlsx")
St <- stress_all$stress_daily

crypto <- list(BTC=BTC, ETH=ETH, BNB=BNB, XRP=XRP)
stable <- list(USDT=USDT, USDC=USDC, DAI=DAI)

# ---- BIC lag per pair ----
bic_lag <- function(yc, ds, lag.max = 10) {
  d <- na.omit(cbind(yc, ds))
  as.integer(VARselect(d, lag.max = lag.max, type = "const")$selection["SC(n)"])
}

run_dir <- function(y, x, S, p, drop = NULL) {
  if (!is.null(drop)) { y[drop] <- NA; x[drop] <- NA; S[drop] <- NA }
  n <- length(y); L <- function(v,k) c(rep(NA,k), head(v, n-k))
  df <- data.frame(y = y, S = S)
  for (k in 1:p) df[[paste0("yl",k)]] <- L(y,k)
  for (k in 1:p) df[[paste0("xl",k)]] <- L(x,k)
  for (k in 1:p) df[[paste0("xS",k)]] <- L(x,k) * S
  df <- na.omit(df)
  xl <- paste0("xl",1:p); xS <- paste0("xS",1:p); yl <- paste0("yl",1:p)
  m <- lm(as.formula(paste("y ~", paste(c(yl, xl, "S", xS), collapse=" + "))), data = df)
  hac <- floor(4 * (nrow(df)/100)^(2/9))
  V <- NeweyWest(m, lag = hac, prewhite = FALSE, adjust = TRUE)
  lh <- function(terms){ h <- tryCatch(linearHypothesis(m, paste(terms,"= 0"), vcov.=V),
                                       error=function(e) NULL)
  if (is.null(h)) c(F=NA,p=NA) else c(F=h$F[2], p=h$`Pr(>F)`[2]) }
  z <- resid(m)
  list(p=p, all=lh(c(xl,xS)), gamma=lh(xl), phi=lh(xS),
       LB=Box.test(z,lag=10,type="Ljung-Box")$p.value,
       LBsq=Box.test(z^2,lag=10,type="Ljung-Box")$p.value)
}

st <- function(p) ifelse(is.na(p),"",ifelse(p<.01,"***",ifelse(p<.05,"**",ifelse(p<.10,"*",""))))
fp <- function(v) sprintf("%6.3f%-3s(%5.3f)", v["F"], st(v["p"]), v["p"])

collect <- function(drop=NULL, tag=""){
  cat("\n=========== Granger F-tests", tag, "=========== [F***(p)]\n")
  cat(sprintf("%-10s %2s | %-17s %-17s %-17s | %-17s %-17s %-17s\n",
              "pair","p","FWD Granger","FWD baseline g","FWD stress phi",
              "REV Granger","REV baseline g","REV stress phi"))
  rows <- list(); k <- 1
  for (sj in names(stable)) for (ci in names(crypto)){
    p <- bic_lag(crypto[[ci]], stable[[sj]])
    f <- run_dir(crypto[[ci]], stable[[sj]], St, p, drop)
    r <- run_dir(stable[[sj]], crypto[[ci]], St, p, drop)
    cat(sprintf("%-10s %2d | %-17s %-17s %-17s | %-17s %-17s %-17s\n",
                paste0(sj,"-",ci), p, fp(f$all),fp(f$gamma),fp(f$phi), fp(r$all),fp(r$gamma),fp(r$phi)))
    rows[[k]] <- data.frame(pair=paste0(sj,"-",ci), lag=p,
                            f_all_F=f$all["F"], f_all_p=f$all["p"], f_g_F=f$gamma["F"], f_g_p=f$gamma["p"], f_phi_F=f$phi["F"], f_phi_p=f$phi["p"],
                            r_all_F=r$all["F"], r_all_p=r$all["p"], r_g_F=r$gamma["F"], r_g_p=r$gamma["p"], r_phi_F=r$phi["F"], r_phi_p=r$phi["p"],
                            LBf=f$LB, LBf2=f$LBsq, LBr=r$LB, LBr2=r$LBsq); k <- k+1
  }
  do.call(rbind, rows)
}

full <- collect(tag="(FULL SAMPLE, BIC)")

cat("\nLjung-Box full sample -> share rejecting at 5%:  resid:",
    round(mean(full$LBf<0.05,na.rm=TRUE),2), " squared resid:", round(mean(full$LBf2<0.05,na.rm=TRUE),2), "\n")

write.csv(full, "granger_full_BIC_28_06.csv", row.names=FALSE)
cat("\nSaved granger_full_BIC_28_06.csv\n")