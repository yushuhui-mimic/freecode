suppressPackageStartupMessages({library(data.table); library(ggplot2); library(splines)})
setwd(".")  # run from the repository root
outdir <- "hrv_mortality"; dir.create(outdir, showWarnings = FALSE)

prep <- function(f, sev) {
  d <- fread(f)
  d[, sae := as.integer(sae)]
  d[, male := fifelse(gender %in% c("M","Male"), 1L, 0L)]
  d[, hosp := as.integer(hospital_death)]
  d[]
}
m <- prep("off_mimic.csv"); e <- prep("off_eicu.csv")

auc <- function(obs, pred) {
  o <- order(pred, decreasing = TRUE)
  obs <- obs[o]; pred <- pred[o]
  n1 <- sum(obs == 1); n0 <- sum(obs == 0)
  if (n1 == 0 || n0 == 0) return(NA)
  ranks <- rank(pred, ties.method = "average")
  (sum(ranks[obs == 1]) - n1*(n1+1)/2) / (n1*n0)
}

res <- list(); curves <- list(); aucs <- list()

for (nm in c("MIMIC","eICU")) {
  d <- if (nm == "MIMIC") m else e
  sev <- if (nm == "MIMIC") "sofa_first_day" else "apachescore"
  dd <- d[!is.na(hr_sd) & !is.na(hosp) & !is.na(get(sev)) & !is.na(age) & !is.na(male)]
  dd[, hr_sd_z := scale(hr_sd)]

  ## 1) per-SD OR, adjusted
  f1 <- glm(as.formula(sprintf("hosp ~ hr_sd_z + age + male + %s + sae", sev)), data = dd, family = binomial)
  s <- summary(f1)$coefficients
  res[[length(res)+1]] <- data.table(cohort = nm, model = "adjusted per SD",
    OR = exp(s["hr_sd_z",1]), lci = exp(s["hr_sd_z",1]-1.96*s["hr_sd_z",2]),
    uci = exp(s["hr_sd_z",1]+1.96*s["hr_sd_z",2]), p = s["hr_sd_z",4], n = nrow(dd))

  ## 2) tertiles
  dd[, tert := cut(hr_sd, quantile(hr_sd, c(0,.33,.67,1)), include.lowest = TRUE, labels = c("T1","T2","T3"))]
  f2 <- glm(as.formula(sprintf("hosp ~ tert + age + male + %s + sae", sev)), data = dd, family = binomial)
  s2 <- summary(f2)$coefficients
  for (lv in c("T2","T3")) {
    if (paste0("tert", lv) %in% rownames(s2)) {
      res[[length(res)+1]] <- data.table(cohort = nm, model = paste0(lv, " vs T1"),
        OR = exp(s2[paste0("tert",lv),1]), lci = exp(s2[paste0("tert",lv),1]-1.96*s2[paste0("tert",lv),2]),
        uci = exp(s2[paste0("tert",lv),1]+1.96*s2[paste0("tert",lv),2]), p = s2[paste0("tert",lv),4], n = nrow(dd))
    }
  }

  ## 3) RCS (natural spline, 3 df) predicted risk curve
  f3 <- glm(as.formula(sprintf("hosp ~ ns(hr_sd, 3) + age + male + %s + sae", sev)), data = dd, family = binomial)
  grid <- data.table(hr_sd = seq(quantile(dd$hr_sd, .01), quantile(dd$hr_sd, .99), length.out = 100))
  grid[, age := median(dd$age)]; grid[, male := 0]; grid[, sae := 0]
  grid[[sev]] <- median(dd[[sev]], na.rm = TRUE)
  pr <- predict(f3, newdata = grid, type = "link", se.fit = TRUE)
  grid[, fit := plogis(pr$fit)]
  grid[, lo := plogis(pr$fit - 1.96*pr$se.fit)]
  grid[, hi := plogis(pr$fit + 1.96*pr$se.fit)]
  grid[, cohort := nm]
  curves[[length(curves)+1]] <- grid

  ## 4) AUC: base vs base+hr_sd
  fb <- glm(as.formula(sprintf("hosp ~ age + male + %s + sae", sev)), data = dd, family = binomial)
  a0 <- auc(dd$hosp, fitted(fb)); a1 <- auc(dd$hosp, fitted(f1))
  aucs[[length(aucs)+1]] <- data.table(cohort = nm, auc_base = round(a0,4),
    auc_plus_hrsd = round(a1,4), delta = round(a1-a0,4))
}

R <- rbindlist(res); C <- rbindlist(curves); A <- rbindlist(aucs)
print(R); print(A)
fwrite(R, file.path(outdir, "hrv_OR_table.csv"))
fwrite(A, file.path(outdir, "hrv_AUC_table.csv"))
fwrite(C, file.path(outdir, "hrv_RCS_curves.csv"))

p <- ggplot(C, aes(x = hr_sd, y = fit*100, colour = cohort, fill = cohort)) +
  geom_ribbon(aes(ymin = lo*100, ymax = hi*100), alpha = .15, colour = NA) +
  geom_line(size = 1) +
  labs(x = "Heart-rate SD in first 48 h (bpm)", y = "Predicted hospital mortality (%)",
       title = "Dose-response: heart-rate variability and mortality",
       subtitle = "Adjusted for age, sex, severity score and SAE; 1st-99th percentile") +
  theme_bw(base_size = 11)
ggsave(file.path(outdir, "Fig_RCS_hrsd.png"), p, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(outdir, "Fig_RCS_hrsd.pdf"), p, width = 7, height = 4.5)

cat("\nDONE\n")
