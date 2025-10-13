#############################################
## PCoA with top/right density plots
## Distance: Aitchison (CLR + Euclidean)
## PERMANOVA with block permutation (Animal_ID)
## 95% CI via jackknife by Animal_ID
## betadisper and ANOSIM diagnostics
#############################################
library(tidyverse)
library(vegan)
library(ape)
library(gridExtra)
library(ggpubr)
library(permute)
set.seed(1)
nperm <- 9999
esp_clr <- read.csv("Genus-G.csv", header = TRUE, row.names = 1, check.names = FALSE)
Xclr <- t(as.matrix(esp_clr))
meta_raw <- read.csv("Group_G.csv", header = TRUE, row.names = 1, check.names = FALSE)

common <- intersect(rownames(Xclr), rownames(meta_raw))
stopifnot(length(common) > 0)
Xclr <- Xclr[common, , drop = FALSE]
meta <- droplevels(meta_raw[common, , drop = FALSE])

first_present <- function(x, choices) {
  for (nm in choices) if (nm %in% colnames(x)) return(nm)
  NA_character_
}
id_col   <- first_present(meta, c("Animal_ID","AnimalID","animal_id","Animal"))
diet_col <- first_present(meta, c("Diet","diet","Group","group"))
if (is.na(id_col) || is.na(diet_col)) stop("meta 必须包含 Animal_ID 与 Diet 两列（容错大小写已尝试）。")

meta[[id_col]]   <- factor(meta[[id_col]])
meta[[diet_col]] <- factor(meta[[diet_col]])

D_ait <- dist(Xclr, method = "euclidean")

pcoa_plot <- function(D, meta, diet_col, title = "PCoA (Aitchison; Diet)") {
  pc <- ape::pcoa(D)
  ax1 <- round(100 * pc$values$Relative_eig[1], 2)
  ax2 <- round(100 * pc$values$Relative_eig[2], 2)
  df  <- as.data.frame(pc$vectors[, 1:2])
  colnames(df) <- c("PC1","PC2")
  df$SampleID <- rownames(df)
  df <- cbind(df, meta[rownames(df), , drop = FALSE])
  pal <- c("#4a6d8c","#946769")
  base <- ggplot(df, aes(PC1, PC2, colour = .data[[diet_col]])) +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(level = 0.80, alpha = 0.35) +
    stat_ellipse(level = 0.95, linetype = "dashed", alpha = 0.35) +
    labs(x = paste0("PCoA1 (", ax1, "%)"),
         y = paste0("PCoA2 (", ax2, "%)"),
         title = title) +
    scale_color_manual(values = pal[seq_len(nlevels(df[[diet_col]]))]) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  top <- ggplot(df, aes(x = PC1, fill = .data[[diet_col]])) +
    geom_density(alpha = 0.5, colour = "black") +
    scale_fill_manual(values = pal[seq_len(nlevels(df[[diet_col]]))]) +
    theme_void() + theme(legend.position = "none")
  right <- ggplot(df, aes(x = PC2, fill = .data[[diet_col]])) +
    geom_density(alpha = 0.5, colour = "black") +
    scale_fill_manual(values = pal[seq_len(nlevels(df[[diet_col]]))]) +
    coord_flip() + theme_void() + theme(legend.position = "none")
  empty <- ggplot() + theme_void()
  grid.arrange(top, empty, base, right, ncol = 2, nrow = 2,
               widths = c(4, 1), heights = c(1, 4))
  invisible(df)
}
pcoa_plot(D_ait, meta, diet_col)

make_block_perm <- function(id_vec, nperm = 999) {
  id_vec <- factor(id_vec)
  how(
    nperm = nperm,
    plots  = Plots(strata = id_vec, type = "free"),
    within = Within(type = "none")
  )
}

ctrl_block <- make_block_perm(meta[[id_col]], nperm = nperm)
fit_diet <- adonis2(as.formula(paste0("D_ait ~ ", diet_col)),
                    data = meta,
                    permutations = ctrl_block,
                    by = "margin")
print(fit_diet)
if (!dir.exists("beta_tests_out")) dir.create("beta_tests_out")
write.table(as.data.frame(fit_diet),
            "beta_tests_out/PERMANOVA_Aitchison_Diet.tsv",
            sep = "\t", quote = FALSE)

stopifnot(is.factor(meta[[diet_col]]), is.factor(meta[[id_col]]))
stopifnot(all(rownames(meta) == attr(D_ait, "Labels")))

fit_full <- adonis2(D_ait ~ get(diet_col), data = meta, permutations = 0)
R2_full  <- fit_full$R2[1]
cat("Diet R² point estimate:", round(R2_full, 4), "\n")

boot_r2_diet_safe <- function(D, meta, diet_col, id_col, B = 1000, max_tries = 1000) {
  id_df <- unique(data.frame(ID = meta[[id_col]], Diet = meta[[diet_col]]))
  id_df$ID   <- factor(as.character(id_df$ID))
  id_df$Diet <- factor(as.character(id_df$Diet), levels = levels(meta[[diet_col]]))
  diet_lvls   <- levels(id_df$Diet)
  ids_by_diet <- lapply(diet_lvls, function(dl) id_df$ID[id_df$Diet == dl])
  names(ids_by_diet) <- diet_lvls
  n_by_diet <- vapply(ids_by_diet, length, integer(1))
  if (any(n_by_diet == 0)) stop("某个 Diet 组没有 Animal_ID。")
  r2s <- numeric(B)
  Dm  <- as.matrix(D)
  for (b in seq_len(B)) {
    tries <- 0
    repeat {
      tries <- tries + 1
      boot_ids <- unlist(mapply(function(ids, k) sample(ids, size = k, replace = TRUE),
                                ids_by_diet, n_by_diet, SIMPLIFY = TRUE), use.names = FALSE)
      keep <- meta[[id_col]] %in% boot_ids
      metasub <- droplevels(meta[keep, , drop = FALSE])
      if (nlevels(metasub[[diet_col]]) < 2) {
        if (tries >= max_tries) { r2s[b] <- NA_real_; break }
      } else {
        idx  <- which(keep)
        Dsub <- as.dist(Dm[idx, idx])
        fitb <- adonis2(Dsub ~ get(diet_col), data = metasub, permutations = 0)
        r2s[b] <- fitb$R2[1]
        break
      }
    }
  }
  quantile(r2s, c(0.025, 0.975), na.rm = TRUE)
}
ci_diet <- boot_r2_diet_safe(D_ait, meta, diet_col = diet_col, id_col = id_col, B = 999)
cat("Diet R² 95% CI (Aitchison):", paste(round(ci_diet, 4), collapse = " - "), "\n")

bd_ait <- betadisper(D_ait, meta[[diet_col]])
print(permutest(bd_ait, permutations = nperm))

anos <- anosim(D_ait, meta[[diet_col]], permutations = nperm)
cat(sprintf("ANOSIM (Aitchison; Diet): R = %.3f, P = %.4g\n",
            anos$statistic, anos$signif))
