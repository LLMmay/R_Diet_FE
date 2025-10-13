#############################################
## Neutral Community Model (NCM) with inset pie chart
## - Fits Sloan NCM (nlsLM), reports R² and Nm
## - Classifies taxa: Above / Neutral / Below (Wilson 95% CI)
## - Exports stats and draws scatter + prediction + inset pie
#############################################
library(Hmisc)
library(minpack.lm)
library(stats4)
library(grid)
spp<-read.csv('HFE-180.txt',head=T,stringsAsFactors=F,row.names=1,sep = "\t")
spp <- t(spp)
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
taxa_keep <- names(p)
d = 1/N
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
m.fit
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr
write.csv(p, file="p.csv")
write.csv(freq, file="freq.csv")
write.csv(freq.pred, file="freq.pred.csv")
bacnlsALL <- data.frame(p, freq, freq.pred, pred.ci[,2:3])
inter.col <- rep('black', nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower] <- '#A52A2A'
inter.col[bacnlsALL$freq >= bacnlsALL$Upper] <- '#2c4976'
above <- sum(bacnlsALL$freq >= bacnlsALL$Upper)
below <- sum(bacnlsALL$freq <= bacnlsALL$Lower)
neutral <- sum(bacnlsALL$freq > bacnlsALL$Lower & bacnlsALL$freq < bacnlsALL$Upper)
total <- nrow(bacnlsALL)
prop_above <- above / total
prop_below <- below / total
prop_neutral <- neutral / total
cat("Proportion of OTUs above the NCM prediction: ", round(prop_above, 3), "\n")
cat("Proportion of OTUs below the NCM prediction: ", round(prop_below, 3), "\n")
cat("Proportion of OTUs within the NCM prediction (neutral): ", round(prop_neutral, 3), "\n")
grid.newpage()
pushViewport(viewport(h=0.6, w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02), extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq, pch=20, gp=gpar(col=inter.col, cex=0.7, alpha=0.3))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p), bacnlsALL$freq.pred, gp=gpar(col='blue', lwd=2), default='native')
grid.lines(log10(bacnlsALL$p), bacnlsALL$Lower, gp=gpar(col='blue', lwd=2, lty=2), default='native')
grid.lines(log10(bacnlsALL$p), bacnlsALL$Upper, gp=gpar(col='blue', lwd=2, lty=2), default='native')
grid.text(y=unit(0,'npc')-unit(2.5,'lines'), label='Mean Relative Abundance (log10)', gp=gpar(fontface=2))
grid.text(x=unit(0,'npc')-unit(3,'lines'), label='Frequency of Occurrence', gp=gpar(fontface=2), rot=90)
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=", round(Rsqr, 3), "\n", "Nm=", round(coef(m.fit)*N)), 
            x=x[j], y=y[i], just=just)
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)
proportions <- c(prop_above, prop_neutral, prop_below)
labels <- c("Above", "Neutral", "Below")
colors <- c("#2c4976", "grey", "#A52A2A")
pie_vp <- viewport(x=0.85, y=0.85, width=0.2, height=0.2, just=c("center", "center"))
pushViewport(pie_vp)
cumsum_props <- cumsum(proportions)
start_angles <- c(0, cumsum_props[-length(cumsum_props)] * 2 * pi)
end_angles <- cumsum_props * 2 * pi
draw_pie_slice <- function(x, y, r, start, end, fill, col="white") {
  n <- 100
  angles <- seq(start, end, length.out=n)
  x_coords <- c(x, x + r * cos(angles), x)
  y_coords <- c(y, y + r * sin(angles), y)
  grid.polygon(x=x_coords, y=y_coords, gp=gpar(fill=fill, col=col))
}
for (i in 1:length(proportions)) {
  draw_pie_slice(x=0.5, y=0.5, r=0.5, start=start_angles[i], end=end_angles[i], 
                 fill=colors[i], col="white")
}
mid_angles <- (start_angles + end_angles) / 2
for (i in 1:length(proportions)) {
  label_x <- 0.5 + 0.7 * cos(mid_angles[i])
  label_y <- 0.5 + 0.7 * sin(mid_angles[i])
  grid.text(sprintf("%.1f%%", proportions[i] * 100), 
            x=label_x, y=label_y, 
            gp=gpar(fontsize=8))
}
popViewport()
legend_x <- unit(0.85, "npc")
legend_y <- unit(0.65, "npc")
for (i in 1:length(labels)) {
  grid.rect(x=legend_x - unit(0.5, "lines"), 
            y=legend_y - unit((i-1)*1.2, "lines"), 
            width=unit(0.5, "lines"), height=unit(0.5, "lines"),
            gp=gpar(fill=colors[i], col="white"))
  grid.text(labels[i], 
            x=legend_x + unit(1, "lines"), 
            y=legend_y - unit((i-1)*1.2, "lines"), 
            gp=gpar(fontsize=10))
}
set.seed(1)
B <- 1000
boot_Nm  <- rep(NA_real_, B)
boot_R2  <- rep(NA_real_, B)
fit_ncm_once <- function(mat, taxa_ids, m_start) {
  X   <- as.matrix(mat[, taxa_ids, drop = FALSE])
  Nb  <- mean(rowSums(X))
  pm  <- colMeans(X)
  p_b <- pm / Nb
  freq_b <- colMeans(1 * (X > 0))
  eps <- 1e-9
  p_b[p_b <= 0]       <- eps
  p_b[p_b >= 1]       <- 1 - eps
  freq_b[freq_b <= 0] <- eps
  freq_b[freq_b >= 1] <- 1 - eps
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
set.seed(1)
B <- 1000
boot_Nm <- boot_R2 <- rep(NA_real_, B)
n_samp <- nrow(spp)
m_start <- as.numeric(coef(m.fit)[["m"]])
for (b in seq_len(B)) {
  idx_b <- sample.int(n_samp, n_samp, replace = TRUE)
  spp_b <- spp[idx_b, taxa_keep, drop = FALSE]
  tmp <- fit_ncm_once(spp_b, taxa_keep, m_start)
  boot_Nm[b] <- tmp["Nm"]; boot_R2[b] <- tmp["R2"]
}
Nm_CI <- quantile(boot_Nm, c(0.025, 0.975), na.rm = TRUE)
R2_CI <- quantile(boot_R2, c(0.025
