#############################################
## Neutral Community Model (NCM) with inset pie chart
## - Fits Sloan NCM (nlsLM), reports R² and Nm
## - Classifies taxa: Above / Neutral / Below (Wilson 95% CI)
## - Exports stats and draws scatter + prediction + inset pie
#############################################

suppressPackageStartupMessages({
  library(Hmisc)
  library(minpack.lm)
  library(stats4)
  library(grid)
})

in_file <- "HFE-180.txt"   # taxa (rows) x samples (cols) abundance table
out_dir <- "ncm_out"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

spp <- read.csv(in_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
spp <- t(spp)

N   <- mean(rowSums(spp))
p.m <- colMeans(spp)
p.m <- p.m[p.m != 0]
p   <- p.m / N
spp.bi <- 1 * (spp > 0)
freq <- colMeans(spp.bi)
freq <- freq[freq != 0]

C <- merge(p, freq, by = 0)
C <- C[order(C[, 2]), ]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))), ]
p    <- C.0[, 2]
freq <- C.0[, 3]
names(p)    <- C.0[, 1]
names(freq) <- C.0[, 1]
taxa_keep <- names(p)

d <- 1 / N
m.fit <- nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE),
               start = list(m = 0.1))
m.ci <- suppressMessages(confint(m.fit, "m", level = 0.95))

freq.pred <- pbeta(d, N * coef(m.fit) * p, N * coef(m.fit) * (1 - p), lower.tail = FALSE)
pred.ci   <- binconf(freq.pred * nrow(spp), nrow(spp), alpha = 0.05,
                     method = "wilson", return.df = TRUE)
Rsqr <- 1 - sum((freq - freq.pred)^2) / sum((freq - mean(freq))^2)

write.csv(p,          file.path(out_dir, "p.csv"))
write.csv(freq,       file.path(out_dir, "freq.csv"))
write.csv(freq.pred,  file.path(out_dir, "freq.pred.csv"))
write.csv(pred.ci,    file.path(out_dir, "freq.pred.ci_wilson.csv"), row.names = FALSE)

plot_dat <- data.frame(taxon = names(p), p = as.numeric(p), freq = as.numeric(freq),
                       freq.pred = as.numeric(freq.pred),
                       Lower = pred.ci$Lower, Upper = pred.ci$Upper,
                       row.names = NULL)

col_vec <- rep("black", nrow(plot_dat))
col_vec[plot_dat$freq <= plot_dat$Lower] <- "#A52A2A"
col_vec[plot_dat$freq >= plot_dat$Upper] <- "#2c4976"

above   <- sum(plot_dat$freq >= plot_dat$Upper)
below   <- sum(plot_dat$freq <= plot_dat$Lower)
neutral <- sum(plot_dat$freq >  plot_dat$Lower & plot_dat$freq <  plot_dat$Upper)
total   <- nrow(plot_dat)

prop_above   <- above   / total
prop_below   <- below   / total
prop_neutral <- neutral / total

cat(sprintf("Proportion above:  %.3f\n", prop_above))
cat(sprintf("Proportion below:  %.3f\n", prop_below))
cat(sprintf("Proportion neutral: %.3f\n", prop_neutral))
cat(sprintf("R²: %.3f\n", Rsqr))
cat(sprintf("Nm: %.3f\n", as.numeric(coef(m.fit)) * N))

png(file.path(out_dir, "NCM_scatter_with_pie.png"), width = 1800, height = 1500, res = 180)
grid.newpage()
pushViewport(viewport(h = 0.6, w = 0.6))
pushViewport(dataViewport(xData = range(log10(plot_dat$p)), yData = c(0, 1.02), extension = c(0.02, 0)))

grid.rect()
grid.points(log10(plot_dat$p), plot_dat$freq, pch = 20,
            gp = gpar(col = col_vec, cex = 0.7, alpha = 0.3))
grid.yaxis(); grid.xaxis()
grid.lines(log10(plot_dat$p), plot_dat$freq.pred, gp = gpar(col = "blue", lwd = 2), default = "native")
grid.lines(log10(plot_dat$p), plot_dat$Lower,     gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")
grid.lines(log10(plot_dat$p), plot_dat$Upper,     gp = gpar(col = "blue", lwd = 2, lty = 2), default = "native")

grid.text(y = unit(0, "npc") - unit(2.5, "lines"),
          label = "Mean Relative Abundance (log10)", gp = gpar(fontface = 2))
grid.text(x = unit(0, "npc") - unit(3, "lines"),
          label = "Frequency of Occurrence", gp = gpar(fontface = 2), rot = 90)

x <- unit(4 / 5, "npc"); y <- unit(1 / 5, "npc")
grid.text(paste0("R² = ", sprintf("%.3f", Rsqr),
                 "\nNm = ", sprintf("%.0f", as.numeric(coef(m.fit)) * N)),
          x = x, y = y, just = c("centre", "bottom"))

proportions <- c(prop_above, prop_neutral, prop_below)
labels     <- c("Above", "Neutral", "Below")
colors     <- c("#2c4976", "grey", "#A52A2A")

pie_vp <- viewport(x = 0.85, y = 0.85, width = 0.2, height = 0.2, just = c("center", "center"))
pushViewport(pie_vp)

cumsum_props  <- cumsum(proportions)
start_angles  <- c(0, head(cumsum_props, -1)) * 2 * pi
end_angles    <- cumsum_props * 2 * pi

draw_pie_slice <- function(x, y, r, start, end, fill, col = "white") {
  n <- 100
  ang <- seq(start, end, length.out = n)
  x_coords <- c(x, x + r * cos(ang), x)
  y_coords <- c(y, y + r * sin(ang), y)
  grid.polygon(x = x_coords, y = y_coords, gp = gpar(fill = fill, col = col))
}
for (i in seq_along(proportions)) {
  draw_pie_slice(0.5, 0.5, 0.5, start_angles[i], end_angles[i], colors[i], "white")
}

mid_angles <- (start_angles + end_angles) / 2
for (i in seq_along(proportions)) {
  lx <- 0.5 + 0.7 * cos(mid_angles[i])
  ly <- 0.5 + 0.7 * sin(mid_angles[i])
  grid.text(sprintf("%.1f%%", proportions[i] * 100), x = lx, y = ly, gp = gpar(fontsize = 8))
}
popViewport()

legend_x <- unit(0.85, "npc"); legend_y <- unit(0.65, "npc")
for (i in seq_along(labels)) {
  grid.rect(x = legend_x - unit(0.5, "lines"),
            y = legend_y - unit((i - 1) * 1.2, "lines"),
            width = unit(0.5, "lines"), height = unit(0.5, "lines"),
            gp = gpar(fill = colors[i], col = "white"))
  grid.text(labels[i],
            x = legend_x + unit(1, "lines"),
            y = legend_y - unit((i - 1) * 1.2, "lines"),
            gp = gpar(fontsize = 10))
}
dev.off()

set.seed(1)
B <- 1000
boot_Nm <- boot_R2 <- rep(NA_real_, B)
n_samp <- nrow(spp)
m_start <- as.numeric(coef(m.fit)[["m"]])

fit_ncm_once <- function(mat, taxa_ids, m_start) {
  X   <- as.matrix(mat[, taxa_ids, drop = FALSE])
  Nb  <- mean(rowSums(X))
  pm  <- colMeans(X)
  p_b <- pm / Nb
  freq_b <- colMeans(1 * (X > 0))

  eps <- 1e-9
  p_b[p_b <= 0]       <- eps; p_b[p_b >= 1]       <- 1 - eps
  freq_b[freq_b <= 0] <- eps; freq_b[freq_b >= 1] <- 1 - eps

  d_b <- 1 / Nb
  fit_try <- try(
    nlsLM(freq_b ~ pbeta(d_b, Nb * m * p_b, Nb * m * (1 - p_b), lower.tail = FALSE),
          start = list(m = m_start), lower = c(m = 1e-8), upper = c(m = 1),
          control = nls.lm.control(maxiter = 500)),
    silent = TRUE
  )
  if (inherits(fit_try, "try-error")) return(c(Nm = NA_real_, R2 = NA_real_))
  mhat <- as.numeric(coef(fit_try))
  pred <- pbeta(d_b, Nb * mhat * p_b, Nb * mhat * (1 - p_b), lower.tail = FALSE)
  R2   <- 1 - sum((freq_b - pred)^2) / sum((freq_b - mean(freq_b))^2)
  Nm   <- Nb * mhat
  c(Nm = Nm, R2 = R2)
}

for (b in seq_len(B)) {
  idx_b <- sample.int(n_samp, n_samp, replace = TRUE)
  spp_b <- spp[idx_b, taxa_keep, drop = FALSE]
  tmp <- fit_ncm_once(spp_b, taxa_keep, m_start)
  boot_Nm[b] <- tmp["Nm"]; boot_R2[b] <- tmp["R2"]
}

Nm_CI <- quantile(boot_Nm, c(0.025, 0.975), na.rm = TRUE)
R2_CI <- quantile(boot_R2, c(0.025, 0.975), na.rm = TRUE)

cat("\n=== Bootstrap 95% CIs ===\n")
cat("Nm (point): ", round(as.numeric(coef(m.fit)) * N), "\n", sep = "")
cat("Nm 95% CI: [", round(Nm_CI[1]), ", ", round(Nm_CI[2]), "]\n", sep = "")
cat("R² (point): ", round(Rsqr, 3), "\n", sep = "")
cat("R² 95% CI: [", round(R2_CI[1], 3), ", ", round(R2_CI[2], 3), "]\n", sep = "")
