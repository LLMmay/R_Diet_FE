library(phyloseq)
library(tidyverse)
library(ggpubr)
library(vegan)
library(Maaslin2)
library(ggplot2)
library(ggprism)
library(ggrepel)
library(compositions) 
library(stringr)
setwd("C:/Users/Lin limei/Desktop/CAZy")
outdir <- "C:/Users/Lin limei/Desktop/MAG_clr/CAZy_FE_G"
input_data <- read.csv("CAZy_G.csv", header = TRUE, row.names = 1, check.names = FALSE)
input_metadata <- read.csv("Group_G.csv",  header = TRUE, row.names = 1, check.names = FALSE)

fix_ids <- function(x) {
  x <- iconv(x, from = "", to = "UTF-8", sub = "")
  x <- gsub("\\s+", "_", x)
  x <- gsub("[^[:alnum:]_.-]", "_", x)
  make.names(x, unique = TRUE, allow_ = TRUE)
}

colnames(input_data)  <- fix_ids(colnames(input_data))
rownames(input_metadata) <- fix_ids(rownames(input_metadata))
rownames(input_data) <- fix_ids(rownames(input_data))
rownames(input_data) <- fix_ids(rownames(input_data))

maas <- Maaslin2(
  input_data       = input_data,
  input_metadata   = input_metadata,
  output           = outdir,
  min_abundance    = 0.0,
  min_prevalence   = 0,
  min_variance     = 0.0,
  normalization    = "CLR",
  transform        = "NONE",
  analysis_method  = "LM",
  max_significance = 0.05,
  fixed_effects    = c("Diet","Time"),
  random_effects   = c("Animal"),
  reference        = c("Diet,F"),
  correction       = "BH",
  standardize      = FALSE,
  cores            = 1,
  plot_heatmap     = FALSE,
  plot_scatter     = FALSE,
  save_scatter     = FALSE,
  save_models      = FALSE
)

data <- maas$results
head(data)
