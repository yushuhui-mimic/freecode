suppressPackageStartupMessages({library(data.table)})
setwd(".")  # run from the repository root

prep <- function(f, sev) {
  d <- fread(f)
  d[, sae := as.integer(sae)]
  d[, male := fifelse(gender %in% c("M","Male"), 1L, 0L)]
  d[, hosp := as.integer(hospital_death)]
  d[]
}
m <- prep("off_mimic.csv"); e <- prep("off_eicu.csv")

cat("MIMIC:", nrow(m), "SAE", sum(m$sae), sprintf("(%.1f%%)", 100*mean(m$sae)),
    "| mortality", sprintf("%.1f%%", 100*mean(m$hosp, na.rm=TRUE)), "\n")
cat("eICU :", nrow(e), "SAE", sum(e$sae), sprintf("(%.1f%%)", 100*mean(e$sae)),
    "| mortality", sprintf("%.1f%%", 100*mean(e$hosp, na.rm=TRUE)), "\n")

feats <- c("hr_mean","hr_sd","hr_cv_pct","hr_instability","map_mean","map_sd","corr_hr_map")

cat("\n=== SAE vs 非 SAE (SMD) ===\n")
for (nm in c("MIMIC","eICU")) {
  d <- if (nm=="MIMIC") m else e
  for (v in feats) {
    a <- d[sae==1][[v]]; b <- d[sae==0][[v]]
    smd <- (mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE))/sqrt((sd(a,na.rm=TRUE)^2+sd(b,na.rm=TRUE)^2)/2)
    cat(sprintf("%-6s %-16s SAE %.3f  non %.3f  SMD %+.3f  p %.2e\n", nm, v,
        mean(a,na.rm=TRUE), mean(b,na.rm=TRUE), smd,
        tryCatch(wilcox.test(a,b)$p.value, error=function(x) NA)))
  }
}

cat("\n=== 校正后 OR (SAE ~ metric + age + sex + severity) ===\n")
for (nm in c("MIMIC","eICU")) {
  d <- if (nm=="MIMIC") m else e
  sev <- if (nm=="MIMIC") "sofa_first_day" else "apachescore"
  for (v in feats) {
    dd <- d[!is.na(get(v)) & !is.na(get(sev))]
    if (nrow(dd) < 300) next
    fit <- glm(as.formula(sprintf("sae ~ scale(%s) + age + male + %s", v, sev)), data=dd, family=binomial)
    s <- summary(fit)$coefficients; i <- grep("scale\\(", rownames(s))[1]
    cat(sprintf("%-6s %-16s OR %.2f (%.2f-%.2f) p %.2e n=%d\n", nm, v,
        exp(s[i,1]), exp(s[i,1]-1.96*s[i,2]), exp(s[i,1]+1.96*s[i,2]), s[i,4], nrow(dd)))
  }
}

cat("\n=== 校正后 OR (死亡 ~ metric + age + sex + severity + SAE) ===\n")
for (nm in c("MIMIC","eICU")) {
  d <- if (nm=="MIMIC") m else e
  sev <- if (nm=="MIMIC") "sofa_first_day" else "apachescore"
  for (v in c("hr_sd","hr_cv_pct","hr_instability","map_sd","corr_hr_map")) {
    dd <- d[!is.na(get(v)) & !is.na(hosp) & !is.na(get(sev))]
    if (nrow(dd) < 300) next
    fit <- glm(as.formula(sprintf("hosp ~ scale(%s) + age + male + %s + sae", v, sev)), data=dd, family=binomial)
    s <- summary(fit)$coefficients; i <- grep("scale\\(", rownames(s))[1]
    cat(sprintf("%-6s %-16s OR %.2f (%.2f-%.2f) p %.2e\n", nm, v,
        exp(s[i,1]), exp(s[i,1]-1.96*s[i,2]), exp(s[i,1]+1.96*s[i,2]), s[i,4]))
  }
}

cat("\n=== 三分位死亡率 ===\n")
for (nm in c("MIMIC","eICU")) {
  d <- if (nm=="MIMIC") m else e
  for (v in c("hr_sd","corr_hr_map")) {
    x <- d[!is.na(get(v)) & !is.na(hosp)]
    br <- unique(quantile(x[[v]], c(0,.33,.67,1), na.rm=TRUE))
    if (length(br) < 4) next
    x[, g := cut(get(v), br, include.lowest=TRUE, labels=c("Low","Mid","High"))]
    cat("--", nm, v, "\n"); print(x[, .(n=.N, mort=round(100*mean(hosp),1), sae=round(100*mean(sae),1)), by=g][order(g)])
  }
}
cat("\nDONE\n")
