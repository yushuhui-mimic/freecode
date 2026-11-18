suppressPackageStartupMessages({library(data.table); library(poLCA)})
setwd(".")  # run from the repository root

LABS <- c("creatinine_max","platelets_min","wbc_max","hemoglobin_min","lactate_max","bilirubin_max")
base <- "./output"

d <- fread(file.path(base, "lab_cluster_assignment_mice.csv"))
k2 <- fread(file.path(base, "k2_assignment.csv"))
for (c in c("patient_id","stay_id")) { d[[c]] <- as.character(d[[c]]); k2[[c]] <- as.character(k2[[c]]) }
d <- merge(d, k2[, .(source, patient_id, stay_id, k2)], by = c("source","patient_id","stay_id"))
cat("merged:", nrow(d), "\n")

# 每个指标按库内三分位切成 3 类（LCA 需要分类指示变量）
for (v in LABS) {
  d[, (paste0(v, "_c")) := {
    br <- quantile(get(v), c(0, 1/3, 2/3, 1), na.rm = TRUE)
    br <- unique(br)
    if (length(br) < 4) as.integer(cut(get(v), 2, labels = FALSE)) + 1L
    else as.integer(cut(get(v), br, include.lowest = TRUE, labels = FALSE))
  }, by = source]
}
ind <- paste0(LABS, "_c")
f <- as.formula(paste("cbind(", paste(ind, collapse = ","), ") ~ 1"))

res <- list()
for (k in 2:4) {
  set.seed(42)
  m <- poLCA(f, data = as.data.frame(d[, ..ind]), nclass = k, maxiter = 300,
             nrep = 3, verbose = FALSE, na.rm = FALSE)
  # entropy
  pp <- m$posterior
  ent <- -sum(pp * log(pp + 1e-12)) / (nrow(pp) * log(k))
  ent <- 1 - ent
  d[, lca := m$predclass]
  tab <- d[, .(n = .N, mort = round(100*mean(hosp_mort),1)), by = lca][order(mort)]
  res[[length(res)+1]] <- data.table(k = k, BIC = round(m$bic,1), AIC = round(m$aic,1),
                                     entropy = round(ent,3), npar = m$npar)
  cat(sprintf("\n=== k=%d  BIC=%.0f  AIC=%.0f  entropy=%.3f  npar=%d\n", k, m$bic, m$aic, ent, m$npar))
  print(tab)
  if (k == 2) {
    # 与 k-means 的一致性
    tb <- table(d$k2, d$lca)
    print(tb)
    p0 <- sum(diag(tb)) / sum(tb)
    pe <- sum(rowSums(tb) * colSums(tb)) / sum(tb)^2
    cat(sprintf("与 k-means k=2 的一致率 = %.3f, Kappa = %.3f\n", p0, (p0-pe)/(1-pe)))
    # 每库 LCA 高危类死亡率
    print(d[, .(n=.N, mort=round(100*mean(hosp_mort),1)), by=.(source, lca)][order(source, lca)])
  }
}
fwrite(rbindlist(res), file.path(base, "LCA_fit_indices.csv"))
cat("\nDONE\n")
