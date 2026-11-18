## ============================================================
## BMJ 风格图件（Arial / Okabe-Ito / 无网格 / 300 dpi）
## 输出：PDF 矢量 + PNG 300dpi
## ============================================================
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(scales); library(grid)
})

BASE <- "./output"
OUT  <- file.path(BASE, "BMJ_figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 全局样式 ----
windowsFonts(Arial = windowsFont("Arial"))
OK <- c(blue = "#0072B2", orange = "#E69F00", green = "#009E73",
        vermillion = "#D55E00", purple = "#CC79A7", sky = "#56B4E9",
        yellow = "#F0E442", black = "#000000")

theme_bmj <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line        = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.35, colour = "black"),
      axis.text        = element_text(colour = "black", size = base_size),
      axis.title       = element_text(colour = "black", size = base_size),
      plot.title       = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle    = element_text(size = base_size, hjust = 0),
      plot.caption     = element_text(size = base_size - 1, hjust = 0),
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size - 1),
      legend.key       = element_blank(),
      panel.grid       = element_blank(),
      strip.background = element_blank(),
      plot.margin      = margin(6, 8, 6, 6)
    )
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, family = "Arial")
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300, bg = "white")
}

## ============================================================
## 数据
## ============================================================
prof <- fread(file.path(BASE, "cluster_profiles_mice.csv"))
ors  <- fread(file.path(BASE, "per_source_OR_k2.csv"))
mort <- fread(file.path(BASE, "mortality_by_cluster_source_mice.csv"))
sens <- fread(file.path(BASE, "sensitivity_summary.csv"))
s56  <- fread(file.path(BASE, "sensitivity_5vs6_labs.csv"))

## ============================================================
## Figure 1  Study flow（STROBE 风格）
## ============================================================
flow <- data.table(
  x = c(3, 3, 3, 3, 1.6, 4.4, 1.6, 4.4),
  y = c(9.2, 7.4, 5.6, 3.6, 2.0, 2.0, 0.9, 0.9),
  w = c(4.6, 4.6, 4.6, 4.6, 2.4, 2.4, 2.4, 2.4),
  h = c(0.95, 1.15, 1.15, 1.0, 0.95, 0.95, 0.75, 0.75),
  lab = c(
    "Critical care databases screened\nMIMIC-IV, eICU-CRD, SICdb, AmsterdamUMCdb,\nNWICU, HiRID, Zigong (n = 136,257 episodes)",
    "Excluded (n = 57,109)\nProxy sepsis definitions only\n(NWICU, HiRID)",
    "Excluded (n = 2,702)\nMissing hospital mortality\nor all six laboratory variables",
    "Sepsis episodes analysed\nn = 76,446 (five databases)",
    "Phenotype 1\nn = 52,310 (68.4%)",
    "Phenotype 2\nn = 24,136 (31.6%)",
    "Hospital mortality 9.1%",
    "Hospital mortality 22.0%"
  )
)
arrows <- data.table(
  x  = c(3, 3, 3, 3, 1.6, 4.4),
  y  = c(8.72, 6.82, 5.02, 3.10, 1.52, 1.52),
  xe = c(3, 3, 3, 3, 1.6, 4.4),
  ye = c(7.90, 6.05, 4.50, 2.45, 1.05, 1.05)
)

p1 <- ggplot() +
  geom_rect(data = flow, aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = "white", colour = "black", linewidth = 0.35) +
  geom_text(data = flow, aes(x = x, y = y, label = lab), size = 2.5, family = "Arial",
            lineheight = 0.95) +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xe, yend = ye),
               arrow = arrow(length = unit(1.6, "mm"), type = "closed"),
               linewidth = 0.35) +
  coord_cartesian(xlim = c(0, 6), ylim = c(0.3, 9.9), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(4, 4, 4, 4))
save_fig(p1, "Figure1_flow", 6.8, 5.4)

## ============================================================
## Figure 2  Phenotype profiles (heatmap)
## ============================================================
nice <- c(creatinine_max="Creatinine", platelets_min="Platelets", wbc_max="White cells",
          hemoglobin_min="Haemoglobin", lactate_max="Lactate", bilirubin_max="Bilirubin")
hm <- melt(prof, id.vars = "cluster", variable.name = "lab", value.name = "z")
hm[, lab := factor(nice[lab], levels = unname(nice))]
hm[, cluster_lab := factor(paste0("Phenotype ", cluster), levels = c("Phenotype 1","Phenotype 2"))]

p2 <- ggplot(hm, aes(x = lab, y = cluster_lab, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", z)),
            size = 2.6, family = "Arial",
            colour = ifelse(abs(hm$z) > 0.55, "white", "black")) +
  scale_fill_gradient2(low = OK[["blue"]], mid = "white", high = OK[["vermillion"]],
                       midpoint = 0, name = "z-score") +
  labs(x = NULL, y = NULL) +
  theme_bmj(9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        axis.line = element_blank(), axis.ticks = element_blank(),
        legend.position = "right")
save_fig(p2, "Figure2_phenotype_heatmap", 5.2, 2.4)

## ============================================================
## Figure 3  Forest（逐库 + 合并）
## ============================================================
f3 <- ors[cluster == 2, .(source, OR, lci, uci)]
setorder(f3, OR)
f3[, src := factor(source, levels = source)]

## DerSimonian–Laird
dl <- function(logOR, se) {
  w <- 1/se^2; fe <- sum(w*logOR)/sum(w)
  Q <- sum(w*(logOR-fe)^2); df <- length(logOR)-1
  C <- sum(w) - sum(w^2)/sum(w)
  tau2 <- max(0, (Q-df)/C)
  wr <- 1/(se^2+tau2)
  est <- sum(wr*logOR)/sum(wr); s <- sqrt(1/sum(wr))
  c(est=est, se=s, I2=max(0,(Q-df)/Q*100))
}
m3 <- dl(log(f3$OR), (log(f3$uci)-log(f3$lci))/(2*1.96))
pooled <- data.table(src = "Pooled (random effects)", OR = exp(m3["est"]),
                     lci = exp(m3["est"]-1.96*m3["se"]), uci = exp(m3["est"]+1.96*m3["se"]))
plotdf <- rbind(f3[, .(src, OR, lci, uci)], pooled, fill = TRUE)
plotdf[, src := factor(src, levels = c("Pooled (random effects)", as.character(f3$src)))]
plotdf[, is_pool := src == "Pooled (random effects)"]

p3 <- ggplot(plotdf, aes(x = OR, y = src)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0, linewidth = 0.4,
                 colour = OK[["blue"]]) +
  geom_point(aes(shape = is_pool), size = 2.1, colour = OK[["blue"]], fill = OK[["blue"]]) +
  geom_text(aes(x = uci, label = sprintf("%.2f (%.2f–%.2f)", OR, lci, uci)),
            hjust = -0.10, size = 2.3, family = "Arial") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 23), guide = "none") +
  scale_x_log10(breaks = c(1, 2, 3, 5, 7), labels = c("1", "2", "3", "5", "7")) +
  coord_cartesian(xlim = c(0.9, 26)) +
  labs(x = "Odds ratio for hospital mortality (phenotype 2 vs 1)", y = NULL) +
  theme_bmj(9)
save_fig(p3, "Figure3_forest_OR", 6.2, 2.9)

## ============================================================
## Figure 4  Mortality by phenotype and database
## ============================================================
m4 <- copy(mort)
setnames(m4, c("cluster","mean"), c("phenotype","mortality"), skip_absent = TRUE)
if (!"mortality" %in% names(m4)) m4[, mortality := mean]
m4[, mortality := mortality * 100]
m4[, phenotype := factor(paste0("Phenotype ", phenotype), levels = c("Phenotype 1","Phenotype 2"))]
ord <- m4[phenotype == "Phenotype 2"][order(mortality), source]
m4[, source := factor(source, levels = ord)]

p4 <- ggplot(m4, aes(x = source, y = mortality, fill = phenotype)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62,
           colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f", mortality)),
            position = position_dodge(width = 0.7), vjust = -0.45, size = 2.3, family = "Arial") +
  scale_fill_manual(values = c("Phenotype 1" = OK[["sky"]], "Phenotype 2" = OK[["vermillion"]])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Hospital mortality (%)", fill = NULL) +
  theme_bmj(9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = c(0.82, 0.85))
save_fig(p4, "Figure4_mortality_bars", 5.6, 3.0)

## ============================================================
## Figure 5  Sensitivity analyses
## ============================================================
lab_map <- c("主分析 k=2" = "Main analysis (k = 2)",
             "S1 排除 ziGongCC" = "Excluding Zigong",
             "S2 仅主分析集" = "Primary analysis set only",
             "S3 缩尾 0.5%/99.5%" = "Winsorising 0.5th/99.5th",
             "S4 k=3" = "Three clusters",
             "S5 k=4" = "Four clusters")
s5 <- copy(sens)
s5[, scenario_en := lab_map[scenario]]
s5[is.na(scenario_en), scenario_en := scenario]
s5[, scenario_en := factor(scenario_en, levels = rev(s5$scenario_en))]

p5 <- ggplot(s5, aes(x = pooled_OR, y = scenario_en)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0, linewidth = 0.4, colour = OK[["green"]]) +
  geom_point(size = 2.0, shape = 15, colour = OK[["green"]]) +
  geom_text(aes(x = uci, label = sprintf("%.2f (%.2f–%.2f)", pooled_OR, lci, uci)),
            hjust = -0.10, size = 2.3, family = "Arial") +
  scale_x_log10(breaks = c(1, 1.5, 2, 3, 5), labels = c("1", "1.5", "2", "3", "5")) +
  coord_cartesian(xlim = c(0.9, 14)) +
  labs(x = "Pooled odds ratio for hospital mortality", y = NULL) +
  theme_bmj(9)
save_fig(p5, "Figure5_sensitivity", 6.4, 2.6)

## ============================================================
## Figure 6  Six vs five variables
## ============================================================
s6 <- dcast(s56[cluster == 2], source ~ model, value.var = "OR")
setnames(s6, c("6 labs", "5 labs"), c("six", "five"))
setorder(s6, six)
s6[, source := factor(source, levels = source)]
long6 <- melt(s6, id.vars = "source", variable.name = "set", value.name = "OR")
long6[, set := factor(fifelse(set == "six", "Six variables (with lactate)",
                              "Five variables (lactate excluded)"),
                      levels = c("Six variables (with lactate)", "Five variables (lactate excluded)"))]

p6 <- ggplot(long6, aes(x = source, y = OR, fill = set)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62,
           colour = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.1f", OR)), position = position_dodge(width = 0.7),
            vjust = -0.45, size = 2.3, family = "Arial") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  scale_fill_manual(values = c("Six variables (with lactate)" = OK[["blue"]],
                               "Five variables (lactate excluded)" = OK[["orange"]])) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Odds ratio for hospital mortality", fill = NULL) +
  theme_bmj(9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = c(0.75, 0.85))
save_fig(p6, "Figure6_six_vs_five", 5.8, 3.0)

cat("Figures written to", OUT, "\n")
cat("Pooled OR:", round(exp(m3["est"]), 2), "(",
    round(exp(m3["est"]-1.96*m3["se"]), 2), "-",
    round(exp(m3["est"]+1.96*m3["se"]), 2), ") I2 =", round(m3["I2"], 1), "%\n")
