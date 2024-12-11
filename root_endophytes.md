---
title: "Root endophyte analysis"
output: 
  html_document:
    keep_md: true
    cache: false
---

Built with R version  4.3.1.



## Setup

### Libraries


```r
# library(car)
library(cowplot)
library(data.table)
library(DESeq2)
library(DHARMa)
library(ggpubr)
library(grid)
library(gridExtra)
library(iNEXT)
library(kableExtra)
library(knitr)
library(lmPerm)
library(MASS)
library(pscl)
# library(rcompanion)
library(seqinr)
library(tidyverse)
library(vegan)
library(viridis)

# devtools::install_github("eastmallingresearch/Metabarcoding_pipeline/scripts")
library(metafuncs)
```

### Functions and constants


```r
ALPHA <-        0.1   # DESeq2 alpha value
OTUFILTER <-    0.01  # Remove OTUs with proportion of total reads below value
READFILTER <-   0.05  # Remove samples with read sum < sample_median_reads*READFILTER 
PAIREDONLY <-   FALSE # Will remove the pair of samples which fail the readfilter - probably only useful for DESeq separated by type NOTE removes pairs before DESeq object is created   
TAXCONF <-      0.80  # Sets the taxonomy confidence level to get "rank" in taxonomy files
TOPOTU <-       10    # Number of Top OTUs for summary information
DIFFOTU <-      200   # Number of Top OTUs for correlation analysis
FUNASVFILTER <- 0 # 0.05  # Filter fungal ASVs with proportion of total reads below value
BACASVFILTER <- 0 # 0.5   # Filter bacterial ASVs with proportion of total reads below value

# graphics
DEVICE <- "png"
DPI <-    1200
WIDTH <-  9
HEIGHT <- 9

# Model design
FACTORS <-            c("Site", "Storage", "Scion")
DESIGN <-             y ~ Site + Storage + Scion
FULL_DESIGN <-        y ~ Site * Storage * Scion
design_with_canker <- y ~ Site * Storage * Scion * Cankers
canker_design <-      "Cankers ~ Site * Storage * Scion"

# Control
ASV_MODELS <- FALSE # Toggle for ASV model cells
```


```r
# colour blind palette
cbPalette <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73", 
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

source("functions/metabarcoding.R")
source("functions/loadme.R")
source("functions/rarefaction.R")
```

## Load data

Bacterial and fungal ASV (ZOTU) tables, sample metadata, and taxonomy files are
loaded into named lists using the `loadData` function from Greg's `metafuncs` 
package.

Site names are encoded as follows according to previous work:

- Avalon (East Sutton) -> 1
- Scripps (Brenchley) -> 2
- WWF (Pluckley) -> 3


```r
metadata <- "sample_metadata.txt"

# Load data
ubiome_FUN <- loadData(
  "data/FUN.zotu_table.txt", metadata,
  "data/zFUN.sintax.taxa", RHB="FUN"
)

ubiome_BAC <- loadData(
  "data/BAC.zotu_table.txt", metadata,
  "data/zBAC.sintax.taxa", RHB="BAC"
)

# Change sites Avalon -> 1, Scripps -> 2, and WWF -> 3.
# Storage from planting date.
# No storage for December plantings, yes for March and April (4 months).
mutate_factors <- function(data){
  data <- data %>%
    rename(location = site, Scion = cultivar) %>%
    mutate(
      Site = case_when(
        location == "Avalon" ~ 1,
        location == "Scripps" ~ 2,
        location == "WWF" ~ 3
      ) %>% as.factor(),
      Storage = case_when(
        planting_date %in% c("march", "april") ~ "yes",
        planting_date %in% c("dec") ~ "no"
      ) %>% as.factor(),
      Scion = as.factor(Scion)
    )
  return(data)
}

ubiome_FUN$colData <- mutate_factors(ubiome_FUN$colData)
ubiome_BAC$colData <- mutate_factors(ubiome_BAC$colData)

# In taxData and countData replace 'OTU' with 'ASV'
rownames(ubiome_FUN$taxData) <- gsub("OTU", "ASV", rownames(ubiome_FUN$taxData))
rownames(ubiome_BAC$taxData) <- gsub("OTU", "ASV", rownames(ubiome_BAC$taxData))

rownames(ubiome_FUN$countData) <- gsub("OTU", "ASV", rownames(ubiome_FUN$countData))
rownames(ubiome_BAC$countData) <- gsub("OTU", "ASV", rownames(ubiome_BAC$countData))
```

### Global removals


```r
# Sample "A2-7" removed due to missampling.
ubiome_BAC$colData <- ubiome_BAC$colData[!rownames(ubiome_BAC$colData) %in% "HMA27", ]
ubiome_BAC$countData <- ubiome_BAC$countData[, !colnames(ubiome_BAC$countData) %in% "HMA27"]
ubiome_FUN$colData <- ubiome_FUN$colData[!rownames(ubiome_FUN$colData) %in% "HMA27", ]
ubiome_FUN$countData <- ubiome_FUN$countData[, !colnames(ubiome_FUN$countData) %in% "HMA27"]
```

## Filter samples and ASVs

### Filtering taxa

Plantae taxa are filtered from fungal `taxData`.
Chloroplast and Eukaryote  taxa are filtered from bacterial `taxData`.
Corresponding ASVs are removed from `countData`.


```r
# Filter Plant, Chloroplast, and Eukaryote ASVs

# Fungi: Plantae ASVs
cat("Fungi:", length(grep("Plantae", ubiome_FUN$taxData$kingdom)), "Plantae ASVs\n")
```

```
# Fungi: 0 Plantae ASVs
```

```r
# Bacteria: Chloroplast (Streptophyta) and Eukaryote ASVs
cat(
  "Bacteria:", length(grep("Streptophyta", ubiome_BAC$taxData$genus)), "Chloroplast ASVs;", 
  length(grep("Eukaryota", ubiome_BAC$taxData$kingdom)), "Eukaryote ASVs\n"
)
```

```
# Bacteria: 37 Chloroplast ASVs; 188 Eukaryote ASVs
```

```r
# Filter Chloroplast and Eukaryote
filt <- rownames(
  ubiome_BAC$taxData[
    grepl("Streptophyta", ubiome_BAC$taxData$genus) & 
    as.numeric(ubiome_BAC$taxData$g_conf) >= TAXCONF,
  ]
)

filt <- c(filt, rownames(ubiome_BAC$taxData[grep("Eukaryota", ubiome_BAC$taxData$kingdom), ]))

cat("Bacteria: removing", length(filt), "ASVs")
```

```
# Bacteria: removing 198 ASVs
```

```r
ubiome_BAC$taxData <- ubiome_BAC$taxData[!rownames(ubiome_BAC$taxData) %in% filt, ]
ubiome_BAC$countData <- ubiome_BAC$countData[!rownames(ubiome_BAC$countData) %in% filt, ]
```

### Filtering samples

Plot rarefaction curves.

Remove samples with read count below 5 % of median.


```r
invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs = list(envir = globalenv())))
rare_bac <- gfunc(countData, colData, "Bacteria")
# rare_bac <- gfunc(as.data.frame(counts(dds)), as.data.frame(colData(dds)), "Bacteria ZOTU")
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs = list(envir = globalenv())))
rare_fun <- gfunc(countData, colData, "Fungi")
# rare_fun <- gfunc(as.data.frame(counts(dds)), as.data.frame(colData(dds)), "Fungi ZOTU")

rarefaction_plots <- grid.arrange(
  rare_bac, rare_fun,
  left = textGrob(label = expression("log"[10] * " aligned sequences"), rot = 90),
  bottom = "ASV count", nrow = 2
)

ggsave(filename = "rarefaction_plots.png", plot = rarefaction_plots, path = "figures/")

rarefaction_plots

# Fungi
med <- median(colSums(ubiome_FUN$countData))
filt <- !colSums(ubiome_FUN$countData) > med * READFILTER
cat("Fungi: ",sum(filt),"sample(s) removed\n")

# Bacteria
med <- median(colSums(ubiome_BAC$countData))
filt <- !colSums(ubiome_BAC$countData) > med * READFILTER
cat("Bacteria: ",sum(filt),"sample(s) removed\n")
```

### Filter ASVs

#### ASV read count

Number of ASVs which account for 50 %, 80 %, and 99 % of total reads.


```r
asv_propotions <- function(countData, proportion){
  i <- sum(countData)
  y <- rowSums(countData)
  y <- y[order(y, decreasing = T)]
  asvs <- length(y[(cumsum(y) / i <= proportion)])
  return(asvs)
}

proportions <- c(0.5, 0.9, 0.99, 1)

top_asvs <- data.table(
  "proportion" = proportions,
  "Fungi" = lapply(proportions, function(x) asv_propotions(ubiome_FUN$countData, x)),
  "Bacteria" = lapply(proportions, function(x) asv_propotions(ubiome_BAC$countData, x))
)

top_asvs %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> proportion </th>
   <th style="text-align:left;"> Fungi </th>
   <th style="text-align:left;"> Bacteria </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 0.50 </td>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:left;"> 169 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.90 </td>
   <td style="text-align:left;"> 171 </td>
   <td style="text-align:left;"> 2186 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 0.99 </td>
   <td style="text-align:left;"> 995 </td>
   <td style="text-align:left;"> 5883 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 1.00 </td>
   <td style="text-align:left;"> 2401 </td>
   <td style="text-align:left;"> 7265 </td>
  </tr>
</tbody>
</table>

#### Filter ASVs

Remove ASVs with read count below 1 % of total reads.


```r
# Fungi
keep <- filter_otus(ubiome_FUN$countData, OTUFILTER)
cat(
  "Fungi:\n", 
  " - total ASVs:", nrow(ubiome_FUN$countData), "\n",
  " - removing", nrow(ubiome_FUN$countData) - length(keep), "ASVs\n",
  " - remaining ASVs:", length(keep), "\n"
)
```

```
# Fungi:
#   - total ASVs: 2401 
#   - removing 1406 ASVs
#   - remaining ASVs: 995
```

```r
ubiome_FUN$taxData <- ubiome_FUN$taxData[rownames(ubiome_FUN$taxData) %in% keep,]
ubiome_FUN$countData <- ubiome_FUN$countData[rownames(ubiome_FUN$countData) %in% keep,]

# Bacteria
keep <-  filter_otus(ubiome_BAC$countData, OTUFILTER)
cat(
  "Bacteria:\n",
  " - total ASVs:", nrow(ubiome_BAC$countData), "\n",
  " - removing", nrow(ubiome_BAC$countData) - length(keep), "ASVs\n",
  " - remaining ASVs:", length(keep), "\n"
)
```

```
# Bacteria:
#   - total ASVs: 7265 
#   - removing 1382 ASVs
#   - remaining ASVs: 5883
```

```r
ubiome_BAC$taxData <- ubiome_BAC$taxData[rownames(ubiome_BAC$taxData) %in% keep,]
ubiome_BAC$countData <- ubiome_BAC$countData[rownames(ubiome_BAC$countData) %in% keep,]
```

## Absolute abundance normalisation

ASV normalisation is performed using qPCR theoretical copy number data.
Copy number is calculated per mg of root sample from the qPCR data.

### Prepare qPCR abundance data


```r
abundance <- fread("mean_abundance.csv")

# Add sample ID to abundance data
abundance$id <- paste0("HM", gsub("-", "", abundance$Sample))
# abundance$id <- abundance$Sample
abundance$copy_number <- abundance$MeanAdjustedTCN_mg
abundance$log_copy_number <- log10(abundance$copy_number)

# Add bacterial (16S) and fungal (ITS) abundance to ubiome BAC and FUN named lists
ubiome_FUN$abundance <- abundance[abundance$Target == "ITS"] %>%
  column_to_rownames(var = "id")
ubiome_BAC$abundance <- abundance[abundance$Target == "16S"] %>%
  column_to_rownames(var = "id")

# Merge copy number from abundance with colData
ubiome_FUN$colData <- merge(
  ubiome_FUN$colData, 
  ubiome_FUN$abundance[, c("Target", "copy_number", "log_copy_number")], 
  by = 0
) %>% column_to_rownames(var = "Row.names")
ubiome_BAC$colData <- merge(
  ubiome_BAC$colData, 
  ubiome_BAC$abundance[, c("Target", "copy_number", "log_copy_number")], 
  by = 0
) %>% column_to_rownames(var = "Row.names")
```

#### Remove outliers


```r
# Detect outliers with std > threshold from the median
detect_outliers <- function(x, val, threshold, na.rm = TRUE) {
  med_x <- median(x[[val]], na.rm = na.rm)
  sd_x <- sd(x[[val]], na.rm = na.rm)
  outliers <- x[x[[val]] > (med_x + threshold * sd_x) | x[[val]] < (med_x - threshold * sd_x), ]
  return(outliers)
}

outliers_FUN <- detect_outliers(ubiome_FUN$abundance, "MeanAdjustedTCN_mg", 3)
outliers_BAC <- detect_outliers(ubiome_BAC$abundance, "MeanAdjustedTCN_mg", 3)

# Remove samples with copy number > 3 std from the median
outliers <- rownames(outliers_FUN)
ubiome_FUN$abundance <- ubiome_FUN$abundance[!rownames(ubiome_FUN$abundance) %in% outliers, ]
ubiome_FUN$countData <- ubiome_FUN$countData[, !colnames(ubiome_FUN$countData) %in% outliers]
ubiome_FUN$colData <- ubiome_FUN$colData[!rownames(ubiome_FUN$colData) %in% outliers, ]

cat("Fungi: removing", length(outliers), "outlier(s)\n")
```

```
# Fungi: removing 1 outlier(s)
```

Sample A1-3 is removed from the fungal data due to abnormally high copy number.

## Canker count data

Canker count data for sampled trees only.


```r
# Canker count data for sampled trees only

canker_data <- fread("canker_data.csv", select = c(1:5, 7:34))

# Remove spaces from column names and convert to lowercase
colnames(canker_data) <- tolower(gsub(" ", "_", colnames(canker_data)))

# Codify site names, add storage and total canker count for timepoint 4
canker_data <- mutate(
  canker_data,
  Site = case_when(
    site == "Avalon" ~ 1,
    site == "Scripps" ~ 2,
    site == "WWF" ~ 3
  ) %>% as.factor(),
  Storage = case_when(
    planting_date %in% c("March", "April") ~ "yes",
    planting_date %in% c("Dec") ~ "no"
  ),
  Scion = as.factor(cultivar),
  total_cankers = a4 + b4 + c4 + d4 + e4
)

# Identify samples with missing values
missing <- unique(canker_data[!complete.cases(canker_data), code])

# Also remove sample A2-7 due to missampling
missing <- c(missing, "HMA27")

# Remove missing samples from canker data
canker_data <- canker_data[!canker_data$code %in% missing, ]

# Verify that there are two trees for each sample
canker_data %>% group_by(code) %>% summarise(n = n()) %>% filter(n != 2)
```

```
# # A tibble: 0 × 2
# # ℹ 2 variables: code <chr>, n <int>
```

```r
# Sum of total cankers for each pair of trees with matching code
cankers <- canker_data %>% 
  group_by(code) %>% 
  summarise(
    Site = first(Site),
    Storage = first(Storage),
    Scion = first(Scion),
    Cankers = sum(total_cankers)
  ) %>% 
  column_to_rownames("code")

# Add total canker count to colData for both FUN and BAC
ubiome_FUN$colData <- merge(
  ubiome_FUN$colData, 
  cankers["Cankers"], 
  by = 0,
  all.x = T
) %>% column_to_rownames("Row.names")

ubiome_BAC$colData <- merge(
  ubiome_BAC$colData, 
  cankers["Cankers"], 
  by = 0,
  all.x = T
) %>% column_to_rownames("Row.names")
```

Summary stats


```r
# png("figures/hist.png", width = 800, height = 600)
# hist(cankers$Cankers, breaks = 20, main = "Total canker count", xlab = "Total canker count")
# dev.off()

cankers_hist <- ggdensity(
  cankers, x = "Cankers", fill = "Site", facet.by = "Site", ncol = 1,
  add = "mean", rug = T, palette = cbPalette,
  title = "Total canker count", xlab = "Total canker count"
)

cankers_hist
```

![](root_endophytes_files/figure-html/canker summary-1.png)<!-- -->

```r
ggsave(filename = "cankers_hist.png", plot = cankers_hist, path = "figures/")

cankers_box <- ggboxplot(
  cankers, x = "Site", y = "Cankers", palette = cbPalette,
  color = "Scion", add = "jitter", legend = "top", 
  title = "Total canker count", xlab = "Site", ylab = "Total canker count"
)

cankers_box
```

![](root_endophytes_files/figure-html/canker summary-2.png)<!-- -->

```r
ggsave(filename = "cankers_box.png", plot = cankers_box, path = "figures/")

cankers_bar <- ggbarplot(
  cankers, x = "Site", y = "Cankers", fill = "Scion", 
  palette = cbPalette, add = "mean_se", position = position_dodge(0.8),
  title = "Total canker count", xlab = "Site", ylab = "Total canker count"
)

cankers_bar
```

![](root_endophytes_files/figure-html/canker summary-3.png)<!-- -->

```r
ggsave(filename = "cankers_bar.png", plot = cankers_bar, path = "figures/")
```

GLM


```r
# Effect of Site, Scion, and Storage on canker count

# Formula
formula <- update(FULL_DESIGN, Cankers ~ .)
# formula <- Cankers ~ Site + Storage + Scion + site:Storage + site:Scion + Storage:Scion

# Log-linear model
canker_lm <- lm(update(FULL_DESIGN, log(Cankers + 1) ~ .), data = cankers)

par(mfrow = c(2, 2))
plot(canker_lm)
```

![](root_endophytes_files/figure-html/canker GLM-1.png)<!-- -->

```r
# Residual checking
res <- resid(canker_lm, type = "pearson")

# Poisson model
canker_poisson <- glm(formula, data = cankers, family = "poisson")

poisson_plot <- plot(simulateResiduals(canker_poisson), title = "Poisson model")
```

![](root_endophytes_files/figure-html/canker GLM-2.png)<!-- -->

```r
# Model overdispersed

# Negative binomial model
canker_negbin <- glm.nb(formula, data = cankers)

sim <- simulateResiduals(canker_negbin)

plot(sim, title = "Negative binomial model")
```

![](root_endophytes_files/figure-html/canker GLM-3.png)<!-- -->

```r
# canker_model_plots <- ggarrange(lm_plot, poisson_plot, negbin_plot, ncol = 3)

# ggsave(filename = "canker_model_plots.png", plot = canker_model_plots, path = "figures/")

# png("figures/canker_residuals.png", width = 800, height = 600)
# plot(sim)
# dev.off()

testZeroInflation(sim)
```

```
# 
# 	DHARMa zero-inflation test via comparison to expected zeros with
# 	simulation under H0 = fitted model
# 
# data:  simulationOutput
# ratioObsSim = 0.77399, p-value = 0.76
# alternative hypothesis: two.sided
```

```r
nagelkerke(canker_negbin)
```

```
# Error in nagelkerke(canker_negbin): could not find function "nagelkerke"
```

```r
# Model good fit

canker_anova <- anova(canker_negbin, test = "Chisq") %>% data.frame()
total_deviance <- sum(canker_anova$Deviance, na.rm = T) + tail(canker_anova$Resid..Dev, 1)
canker_anova$Perc..Dev <- canker_anova$Deviance / total_deviance * 100

canker_anova
```

```
#                    Df     Deviance Resid..Df Resid..Dev     Pr..Chi.
# NULL               NA           NA        73  314.98196           NA
# Site                2 118.75744595        71  196.22452 1.629852e-26
# Storage             1   0.02066259        70  196.20386 8.857019e-01
# Scion               6  24.29375368        64  171.91010 4.611042e-04
# Site:Storage        2  28.86115598        62  143.04895 5.406044e-07
# Site:Scion         12  31.65405909        50  111.39489 1.564383e-03
# Storage:Scion       6   7.40666778        44  103.98822 2.848693e-01
# Site:Storage:Scion 11  26.62949240        33   77.35873 5.225415e-03
#                       Perc..Dev
# NULL                         NA
# Site               37.702935303
# Storage             0.006559927
# Scion               7.712744377
# Site:Storage        9.162796386
# Site:Scion         10.049483066
# Storage:Scion       2.351457744
# Site:Storage:Scion  8.454291190
```

![](root_endophytes_files/figure-html/canker GLM-4.png)<!-- -->

## Create DESeq objects


```r
# Make sure countData and colData still match, if they do, create DESeq objects, if not throw error
if(identical(colnames(ubiome_FUN$countData), rownames(ubiome_FUN$colData))) {
  # Create DESeq object
  ubiome_FUN$dds <- ubiom_to_des(ubiome_FUN)
  print("FUN DESeq object created")
} else {
  stop("FUN countData and colData do not match")
}
```

```
# [1] "FUN DESeq object created"
```

```r
if(identical(colnames(ubiome_BAC$countData), rownames(ubiome_BAC$colData))) {
  # Create DESeq object
  ubiome_BAC$dds <- ubiom_to_des(ubiome_BAC)
  print("BAC DESeq object created")
} else {
  stop("BAC countData and colData do not match")
}
```

```
# [1] "BAC DESeq object created"
```

## Abundance normalisation

Absolute abundance normalisation using DESeq2 size factors.

Values are centred around the mean of the copy number.


```r
# Normalise count data using DESeq2 size factors

ubiome_FUN$dds$sizeFactor <- ubiome_FUN$dds$copy_number / mean(ubiome_FUN$dds$copy_number)
ubiome_BAC$dds$sizeFactor <- ubiome_BAC$dds$copy_number / mean(ubiome_BAC$dds$copy_number)
```


```r
# Save environment
save.image("data_loaded.RData")
```

<!-- #=============================================================================== -->
# **Fungi**
<!-- #=============================================================================== -->


```r
# Unpack fungi data
invisible(mapply(assign, names(ubiome_FUN), ubiome_FUN, MoreArgs = list(envir = globalenv())))
```

## ASV and sample summary

### Read and sample summary


```r
cat(
  "Raw reads", "\n\n",
  "Total raw reads:\t\t", sum(countData), "\n",
  "Mean raw reads per sample:\t", mean(colSums(countData)), "\n",
  "Median raw reads per sample:\t", median(colSums(countData)), "\n",
  "Max raw reads per sample:\t", max(colSums(countData)), "\n",
  "Min raw reads per sample:\t", min(colSums(countData)), "\n\n"
)
```

```
# Raw reads 
# 
#  Total raw reads:		 7293776 
#  Mean raw reads per sample:	 90046.62 
#  Median raw reads per sample:	 93435 
#  Max raw reads per sample:	 113518 
#  Min raw reads per sample:	 38472
```

```r
#colSums(countData)

nct <- counts(dds, normalize = T)
cat("Normalised reads", "\n\n",
  "Total normalised reads:\t\t", sum(nct), "\n",
  "Mean normalised reads per sample:\t", mean(colSums(nct)), "\n",
  "Median normalised reads per sample:\t", median(colSums(nct)), "\n",
  "Min normalised reads per sample:\t", min(colSums(nct)), "\n",
  "Max normalised reads per sample:\t", max(colSums(nct)), "\n\n"
)
```

```
# Normalised reads 
# 
#  Total normalised reads:		 12468857 
#  Mean normalised reads per sample:	 153936.5 
#  Median normalised reads per sample:	 98624.28 
#  Min normalised reads per sample:	 28901.7 
#  Max normalised reads per sample:	 881441.3
```

```r
#round(colSums(counts(dds,normalize = T)),0)
```

### ASV summary 


```r
cat(
  "Total ASVs:\t\t", nrow(taxData),"\n\n",
  "Raw reads per ASV summary", "\n\n",
  "Mean raw reads per ASV:\t", mean(rowSums(countData)),"\n",
  "Median raw per ASV:\t\t", median(rowSums(countData)),"\n",
  "ASV raw Min reads:\t\t", min(rowSums(countData)),"\n",
  "ASV raw Max reads:\t\t", max(rowSums(countData)),"\n\n"
)
```

```
# Total ASVs:		 995 
# 
#  Raw reads per ASV summary 
# 
#  Mean raw reads per ASV:	 7330.428 
#  Median raw per ASV:		 588 
#  ASV raw Min reads:		 115 
#  ASV raw Max reads:		 714327
```

```r
cat(
  "Normalised reads per ASV summary","\n\n",
  "Mean normalised reads per ASV:\t\t", mean(rowSums(nct)),"\n",
  "Median normalised reads per ASV:\t", median(rowSums(nct)),"\n",
  "ASV normalised Min reads:\t\t", min(rowSums(nct)),"\n",
  "ASV normalised Max reads:\t\t", max(rowSums(nct)),"\n\n"
)
```

```
# Normalised reads per ASV summary 
# 
#  Mean normalised reads per ASV:		 12531.51 
#  Median normalised reads per ASV:	 1025.725 
#  ASV normalised Min reads:		 101.2814 
#  ASV normalised Max reads:		 1509459
```

```r
y <- rowSums(nct)
y <- y[order(y, decreasing = T)]
# proportion
xy <- y / sum(y)

cat("Top " ,TOPOTU, "ASVs:\n")
```

```
# Top  10 ASVs:
```

```r
data.frame(
  counts = y[1:TOPOTU], 
  proportion = xy[1:TOPOTU], 
  rank = taxData[names(y)[1:TOPOTU],]$rank
) %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> counts </th>
   <th style="text-align:right;"> proportion </th>
   <th style="text-align:left;"> rank </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ASV2 </td>
   <td style="text-align:right;"> 1509458.8 </td>
   <td style="text-align:right;"> 0.1210583 </td>
   <td style="text-align:left;"> Ascomycota(p) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV1 </td>
   <td style="text-align:right;"> 1490469.5 </td>
   <td style="text-align:right;"> 0.1195354 </td>
   <td style="text-align:left;"> Dactylonectria macrodidyma(s) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV5 </td>
   <td style="text-align:right;"> 1068164.1 </td>
   <td style="text-align:right;"> 0.0856666 </td>
   <td style="text-align:left;"> Leotiomycetes(c) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV4 </td>
   <td style="text-align:right;"> 1059908.0 </td>
   <td style="text-align:right;"> 0.0850044 </td>
   <td style="text-align:left;"> Ascomycota(p) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV3 </td>
   <td style="text-align:right;"> 480660.1 </td>
   <td style="text-align:right;"> 0.0385488 </td>
   <td style="text-align:left;"> Ilyonectria destructans(s) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV7 </td>
   <td style="text-align:right;"> 290896.6 </td>
   <td style="text-align:right;"> 0.0233299 </td>
   <td style="text-align:left;"> Fusarium(g) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV6 </td>
   <td style="text-align:right;"> 227927.9 </td>
   <td style="text-align:right;"> 0.0182798 </td>
   <td style="text-align:left;"> Ilyonectria robusta(s) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV9 </td>
   <td style="text-align:right;"> 201690.6 </td>
   <td style="text-align:right;"> 0.0161755 </td>
   <td style="text-align:left;"> Ascomycota(p) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV8 </td>
   <td style="text-align:right;"> 191083.5 </td>
   <td style="text-align:right;"> 0.0153249 </td>
   <td style="text-align:left;"> Fusarium(g) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV11 </td>
   <td style="text-align:right;"> 131684.2 </td>
   <td style="text-align:right;"> 0.0105610 </td>
   <td style="text-align:left;"> Truncatella angustata(s) </td>
  </tr>
</tbody>
</table>

## Taxonomy Summary

### Taxonomy identifiable

Proportion of ASVs which can be assigned (with the given confidence) at each taxonomic rank.


```r
# Proportion of ASVs which can be assigned (with the given confidence) at each taxonomic rank

tx <- copy(taxData)
setDT(tx)
cols <- names(tx)[9:15]

tx[, (cols) := lapply(.SD, as.factor), .SDcols = cols]

data.table(
  rank = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
  "0.8" = round(unlist(lapply(cols, function(col) sum(as.number(tx[[col]]) >= 0.8) / nrow(tx))), 2),
  "0.65" = round(unlist(lapply(cols, function(col) sum(as.number(tx[[col]]) >= 0.65) / nrow(tx))), 2),
  "0.5" = round(unlist(lapply(cols, function(col) sum(as.number(tx[[col]]) >= 0.5) / nrow(tx))), 2)
) %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> rank </th>
   <th style="text-align:right;"> 0.8 </th>
   <th style="text-align:right;"> 0.65 </th>
   <th style="text-align:right;"> 0.5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> kingdom </td>
   <td style="text-align:right;"> 1.00 </td>
   <td style="text-align:right;"> 1.00 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> phylum </td>
   <td style="text-align:right;"> 0.84 </td>
   <td style="text-align:right;"> 0.87 </td>
   <td style="text-align:right;"> 0.90 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> class </td>
   <td style="text-align:right;"> 0.70 </td>
   <td style="text-align:right;"> 0.74 </td>
   <td style="text-align:right;"> 0.78 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> order </td>
   <td style="text-align:right;"> 0.54 </td>
   <td style="text-align:right;"> 0.60 </td>
   <td style="text-align:right;"> 0.64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> family </td>
   <td style="text-align:right;"> 0.42 </td>
   <td style="text-align:right;"> 0.45 </td>
   <td style="text-align:right;"> 0.49 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> genus </td>
   <td style="text-align:right;"> 0.38 </td>
   <td style="text-align:right;"> 0.42 </td>
   <td style="text-align:right;"> 0.47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> species </td>
   <td style="text-align:right;"> 0.24 </td>
   <td style="text-align:right;"> 0.30 </td>
   <td style="text-align:right;"> 0.35 </td>
  </tr>
</tbody>
</table>

% of reads which can be assigned to each taxonomic ranks


```r
tx <-taxData[rownames(dds),]
nc <- counts(dds, normalize = T)
ac <- sum(nc)

data.table(
  rank = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
  "0.8" = round(unlist(lapply(cols, function(col)(sum(nc[which(as.numeric(tx[[col]]) >= 0.8),]) / ac * 100))), 2),
  "0.65" = round(unlist(lapply(cols, function(col)(sum(nc[which(as.numeric(tx[[col]]) >= 0.65),]) / ac * 100))), 2),
  "0.5" = round(unlist(lapply(cols, function(col)(sum(nc[which(as.numeric(tx[[col]]) >= 0.5),]) / ac * 100))), 2)
) %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> rank </th>
   <th style="text-align:right;"> 0.8 </th>
   <th style="text-align:right;"> 0.65 </th>
   <th style="text-align:right;"> 0.5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> kingdom </td>
   <td style="text-align:right;"> 100.00 </td>
   <td style="text-align:right;"> 100.00 </td>
   <td style="text-align:right;"> 100.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> phylum </td>
   <td style="text-align:right;"> 84.14 </td>
   <td style="text-align:right;"> 96.59 </td>
   <td style="text-align:right;"> 96.83 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> class </td>
   <td style="text-align:right;"> 60.12 </td>
   <td style="text-align:right;"> 70.92 </td>
   <td style="text-align:right;"> 71.47 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> order </td>
   <td style="text-align:right;"> 53.49 </td>
   <td style="text-align:right;"> 58.87 </td>
   <td style="text-align:right;"> 68.76 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> family </td>
   <td style="text-align:right;"> 44.97 </td>
   <td style="text-align:right;"> 46.80 </td>
   <td style="text-align:right;"> 50.25 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> genus </td>
   <td style="text-align:right;"> 46.06 </td>
   <td style="text-align:right;"> 48.01 </td>
   <td style="text-align:right;"> 50.72 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> species </td>
   <td style="text-align:right;"> 30.44 </td>
   <td style="text-align:right;"> 36.70 </td>
   <td style="text-align:right;"> 41.62 </td>
  </tr>
</tbody>
</table>

### Taxonomy plots

Plots of proportion of normalised reads assigned to members of phylum and class.


```r
dat <- list(as.data.frame(counts(dds, normalize = T)), taxData, as.data.frame(colData(dds)))

design <- c("Site", "Storage")

# md1 <- getSummedTaxa(dat, conf = TAXCONF, design = design, cutoff = 0.1)
md1 <- getSummedTaxa(dat, conf = TAXCONF, design = design, taxon = "phylum", cutoff = 0.1)

md1[, Site := factor(Site, levels = c(1, 2, 3))]
md1[, Storage := factor(Storage, levels = c("no", "yes"))]
md1[, taxon := factor(taxon, levels = unique(taxon[order(value, decreasing = T)]))]

removals <- md1[, .(value = mean(value)), by = "taxon"][value < 0.5, taxon]
md1 <- md1[!taxon %in% removals, ]

fun_phylum_plot <- plotfun1(md1, x = "taxon", fill = "Site") +
  facet_wrap(~ Storage)

ggsave("figures/fun_phylum.png", fun_phylum_plot, width = 25, height = 15, units = "cm")

fun_phylum_plot
```

![](root_endophytes_files/figure-html/FUN taxonomy plots-1.png)<!-- -->

```r
md2 <- getSummedTaxa(dat, conf = TAXCONF, design = design, taxon = "class", cutoff = 0.1)

md2[, Site := factor(Site, levels = c(1, 2, 3))]
md2[, Storage := factor(Storage, levels = c("no", "yes"))]
md2[, taxon := factor(taxon, levels = unique(taxon[order(value, decreasing = T)]))]

removals <- md2[, .(value = mean(value)), by = "taxon"][value < 0.5, taxon]
md2 <- md2[!taxon %in% removals, ]

fun_class_plot <- plotfun1(md2, x = "taxon", fill = "Site") +
  facet_wrap(~ Storage)

ggsave("figures/fun_class.png", fun_class_plot, width = 25, height = 15, units = "cm")

fun_class_plot
```

![](root_endophytes_files/figure-html/FUN taxonomy plots-2.png)<!-- -->

## Community size

Plot copy number for each sample grouped by site, Scion, and Storage.
Test the effect of site, Scion, and Storage on copy number using ANOVA.


```r
# abundance_plot <- ggplot(
#   data = as.data.frame(colData(dds)), 
#   aes(x = site, y = log_copy_number, colour = Scion, shape = Storage)
# ) + geom_jitter() + 
#   scale_colour_manual(values = cbPalette)

fun_abundance_box <- ggboxplot(
  data = as.data.frame(colData(dds)), x = "Site", y = "log_copy_number", 
  color = "Scion", add = "jitter", legend = "top", 
  title = "Fungal abundance", xlab = "Site", ylab = "log10 copy number"
)

ggsave(
  filename = "fun_abundance.png", plot = fun_abundance_box, path = "figures/", 
  height = 20, width = 20, units = "cm"
)

fun_abundance_box
```

![](root_endophytes_files/figure-html/FUN size-1.png)<!-- -->

```r
fun_abundance_bar <- ggbarplot(
  data = as.data.frame(colData(dds)), x = "Storage", y = "log_copy_number", 
  fill = "Site", add = "mean_se", 
  palette = cbPalette, position = position_dodge(0.8),
  title = "(a) Fungal abundance", xlab = "Storage ", ylab = "Mean copy number (log10)"
) + guides(fill = guide_legend(title = "Site"))

ggsave(
  filename = "fun_abundance_bar.png", plot = fun_abundance_bar, path = "figures/", 
  height = 20, width = 20, units = "cm"
)

fun_abundance_bar
```

![](root_endophytes_files/figure-html/FUN size-2.png)<!-- -->

```r
# Formula for ANOVA
formula <- update(FULL_DESIGN, log_copy_number ~ .)

abundance_anova <- aov(formula, data = as.data.frame(colData(dds)))

# Normality check
par(mfrow = c(2, 2))
plot(abundance_anova)
```

![](root_endophytes_files/figure-html/FUN size-3.png)<!-- -->

```r
png("figures/fun_abundance_norm.png", width = 800, height = 600)
par(mfrow = c(2, 2))
plot(abundance_anova)
dev.off()
```

```
# png 
#   2
```

```r
# Results
summary(abundance_anova)
```

```
#                    Df Sum Sq Mean Sq F value Pr(>F)  
# Site                2  0.861  0.4306   4.663 0.0153 *
# Storage             1  0.501  0.5012   5.427 0.0251 *
# Scion               6  0.477  0.0795   0.860 0.5324  
# Site:Storage        2  0.683  0.3415   3.698 0.0338 *
# Site:Scion         12  1.203  0.1003   1.086 0.3981  
# Storage:Scion       6  0.458  0.0763   0.827 0.5564  
# Site:Storage:Scion 12  0.918  0.0765   0.828 0.6214  
# Residuals          39  3.602  0.0924                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
abundance_results <- abundance_anova %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(abundance_results$Sum.Sq)
abundance_results$Perc.Var <- abundance_results$Sum.Sq / total_variance * 100

abundance_results %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> Sum.Sq </th>
   <th style="text-align:right;"> Mean.Sq </th>
   <th style="text-align:right;"> F.value </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.8612429 </td>
   <td style="text-align:right;"> 0.4306214 </td>
   <td style="text-align:right;"> 4.6626824 </td>
   <td style="text-align:right;"> 0.0152878 </td>
   <td style="text-align:right;"> 9.895484 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.5012156 </td>
   <td style="text-align:right;"> 0.5012156 </td>
   <td style="text-align:right;"> 5.4270616 </td>
   <td style="text-align:right;"> 0.0250969 </td>
   <td style="text-align:right;"> 5.758852 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.4767020 </td>
   <td style="text-align:right;"> 0.0794503 </td>
   <td style="text-align:right;"> 0.8602723 </td>
   <td style="text-align:right;"> 0.5323908 </td>
   <td style="text-align:right;"> 5.477197 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.6830842 </td>
   <td style="text-align:right;"> 0.3415421 </td>
   <td style="text-align:right;"> 3.6981494 </td>
   <td style="text-align:right;"> 0.0338326 </td>
   <td style="text-align:right;"> 7.848482 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 1.2031665 </td>
   <td style="text-align:right;"> 0.1002639 </td>
   <td style="text-align:right;"> 1.0856371 </td>
   <td style="text-align:right;"> 0.3981291 </td>
   <td style="text-align:right;"> 13.824108 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.4580790 </td>
   <td style="text-align:right;"> 0.0763465 </td>
   <td style="text-align:right;"> 0.8266645 </td>
   <td style="text-align:right;"> 0.5564216 </td>
   <td style="text-align:right;"> 5.263223 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.9180631 </td>
   <td style="text-align:right;"> 0.0765053 </td>
   <td style="text-align:right;"> 0.8283835 </td>
   <td style="text-align:right;"> 0.6213710 </td>
   <td style="text-align:right;"> 10.548335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 3.6018401 </td>
   <td style="text-align:right;"> 0.0923549 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 41.384319 </td>
  </tr>
</tbody>
</table>

### Communtiy size with canker count


```r
cat("Model formula: ", deparse(update(design_with_canker, log_copy_number ~ .)), "\n")
```

```
# Model formula:  log_copy_number ~ Site + Storage + Scion + Cankers + Site:Storage +      Site:Scion + Storage:Scion + Site:Cankers + Storage:Cankers +      Scion:Cankers + Site:Storage:Scion + Site:Storage:Cankers +      Site:Scion:Cankers + Storage:Scion:Cankers + Site:Storage:Scion:Cankers
```

```r
abundance_canker_anova <- aov(update(design_with_canker, log_copy_number ~ .), data = as.data.frame(colData(dds)))
summary(abundance_canker_anova)
```

```
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2 1.0706  0.5353   4.031  0.199
# Storage                     1 0.4432  0.4432   3.338  0.209
# Scion                       6 0.5775  0.0963   0.725  0.679
# Cankers                     1 0.0231  0.0231   0.174  0.717
# Site:Storage                2 0.7982  0.3991   3.006  0.250
# Site:Scion                 12 0.7337  0.0611   0.460  0.843
# Storage:Scion               6 0.4909  0.0818   0.616  0.727
# Site:Cankers                2 0.1441  0.0720   0.543  0.648
# Storage:Cankers             1 0.0010  0.0010   0.007  0.940
# Scion:Cankers               6 0.2230  0.0372   0.280  0.905
# Site:Storage:Scion         11 1.7497  0.1591   1.198  0.540
# Site:Storage:Cankers        2 0.1581  0.0791   0.595  0.627
# Site:Scion:Cankers         11 0.8236  0.0749   0.564  0.785
# Storage:Scion:Cankers       5 0.1443  0.0289   0.217  0.926
# Site:Storage:Scion:Cankers  2 0.0944  0.0472   0.356  0.738
# Residuals                   2 0.2656  0.1328               
# 8 observations deleted due to missingness
```

```r
abundance_canker_results <- abundance_canker_anova %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(abundance_canker_results$Sum.Sq)
abundance_canker_results$Perc.Var <- abundance_canker_results$Sum.Sq / total_variance * 100
abundance_canker_results %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> Sum.Sq </th>
   <th style="text-align:right;"> Mean.Sq </th>
   <th style="text-align:right;"> F.value </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1.0705544 </td>
   <td style="text-align:right;"> 0.5352772 </td>
   <td style="text-align:right;"> 4.0312124 </td>
   <td style="text-align:right;"> 0.1987592 </td>
   <td style="text-align:right;"> 13.8300271 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.4431962 </td>
   <td style="text-align:right;"> 0.4431962 </td>
   <td style="text-align:right;"> 3.3377436 </td>
   <td style="text-align:right;"> 0.2092346 </td>
   <td style="text-align:right;"> 5.7254592 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.5775253 </td>
   <td style="text-align:right;"> 0.0962542 </td>
   <td style="text-align:right;"> 0.7248976 </td>
   <td style="text-align:right;"> 0.6785683 </td>
   <td style="text-align:right;"> 7.4607980 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0231057 </td>
   <td style="text-align:right;"> 0.0231057 </td>
   <td style="text-align:right;"> 0.1740104 </td>
   <td style="text-align:right;"> 0.7170845 </td>
   <td style="text-align:right;"> 0.2984919 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.7981764 </td>
   <td style="text-align:right;"> 0.3990882 </td>
   <td style="text-align:right;"> 3.0055628 </td>
   <td style="text-align:right;"> 0.2496528 </td>
   <td style="text-align:right;"> 10.3112936 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.7336815 </td>
   <td style="text-align:right;"> 0.0611401 </td>
   <td style="text-align:right;"> 0.4604508 </td>
   <td style="text-align:right;"> 0.8433227 </td>
   <td style="text-align:right;"> 9.4781128 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.4908622 </td>
   <td style="text-align:right;"> 0.0818104 </td>
   <td style="text-align:right;"> 0.6161199 </td>
   <td style="text-align:right;"> 0.7267407 </td>
   <td style="text-align:right;"> 6.3412351 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.1440854 </td>
   <td style="text-align:right;"> 0.0720427 </td>
   <td style="text-align:right;"> 0.5425589 </td>
   <td style="text-align:right;"> 0.6482735 </td>
   <td style="text-align:right;"> 1.8613764 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0009707 </td>
   <td style="text-align:right;"> 0.0009707 </td>
   <td style="text-align:right;"> 0.0073101 </td>
   <td style="text-align:right;"> 0.9396531 </td>
   <td style="text-align:right;"> 0.0125395 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.2229625 </td>
   <td style="text-align:right;"> 0.0371604 </td>
   <td style="text-align:right;"> 0.2798579 </td>
   <td style="text-align:right;"> 0.9049341 </td>
   <td style="text-align:right;"> 2.8803557 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 1.7497057 </td>
   <td style="text-align:right;"> 0.1590642 </td>
   <td style="text-align:right;"> 1.1979240 </td>
   <td style="text-align:right;"> 0.5403023 </td>
   <td style="text-align:right;"> 22.6036881 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.1581039 </td>
   <td style="text-align:right;"> 0.0790519 </td>
   <td style="text-align:right;"> 0.5953460 </td>
   <td style="text-align:right;"> 0.6268233 </td>
   <td style="text-align:right;"> 2.0424753 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.8235730 </td>
   <td style="text-align:right;"> 0.0748703 </td>
   <td style="text-align:right;"> 0.5638536 </td>
   <td style="text-align:right;"> 0.7850178 </td>
   <td style="text-align:right;"> 10.6393816 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.1442801 </td>
   <td style="text-align:right;"> 0.0288560 </td>
   <td style="text-align:right;"> 0.2173168 </td>
   <td style="text-align:right;"> 0.9264703 </td>
   <td style="text-align:right;"> 1.8638917 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0944484 </td>
   <td style="text-align:right;"> 0.0472242 </td>
   <td style="text-align:right;"> 0.3556489 </td>
   <td style="text-align:right;"> 0.7376541 </td>
   <td style="text-align:right;"> 1.2201376 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.2655664 </td>
   <td style="text-align:right;"> 0.1327832 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 3.4307364 </td>
  </tr>
</tbody>
</table>

## Alpha diversity analysis

### Alpha diversity plot


```r
# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

fun_alpha_plot <- plot_alpha(
  counts(dds, normalize = F), colData(dds),
  design = "Scion", colour = "Site",
  measures = c("Shannon", "Simpson"),
  type = "bar"
) + scale_colour_manual(values = cbPalette) + 
  theme(axis.title.x = element_blank()) +
  ggtitle("Fungal α-diversity")

ggsave(
  filename = "fun_alpha.png", plot = fun_alpha_plot, path = "figures/", 
  height = 20, width = 40, units = "cm"
)
```

```
# Error in `geom_errorbar()` at metafuncs/R/plot_phylo.R:243:11:
# ! Problem while computing aesthetics.
# ℹ Error occurred in the 2nd layer.
# Caused by error:
# ! object 'se' not found
```

```r
fun_alpha_plot
```

```
# Error in `geom_errorbar()` at metafuncs/R/plot_phylo.R:243:11:
# ! Problem while computing aesthetics.
# ℹ Error occurred in the 2nd layer.
# Caused by error:
# ! object 'se' not found
```

### Permutation based anova on α-diversity index ranks


```r
# get the diversity index data
all_alpha_ord <- plot_alpha(
  counts(dds, normalize = F),
  colData(dds),
  returnData = T
)

# join diversity indices and metadata
all_alpha_ord <- all_alpha_ord[
  as.data.table(colData(dds), keep.rownames = "Samples"), 
  on = "Samples"
]

fun_alpha <- all_alpha_ord

formula <- FULL_DESIGN # x ~ Site * Storage * Scion + Site / Site.block

# Chao1

setkey(all_alpha_ord, S.chao1)
all_alpha_ord[, measure := as.numeric(as.factor(S.chao1))]
result <- aovp(update(formula, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS "
```

```r
summary(result)
```

```
# Component 1 :
#                    Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# Site                2  11554.5    5777.2 5000   <2e-16 ***
# Storage             1   2056.4    2056.4 1827   0.0520 .  
# Site:Storage        2    812.4     406.2  461   0.3731    
# Scion               6    875.7     145.9   80   0.9375    
# Site:Scion         12   2817.7     234.8  566   0.8975    
# Storage:Scion       6   2046.5     341.1  455   0.7560    
# Site:Storage:Scion 12   2735.3     227.9 2678   0.9459    
# Residuals          39  21381.5     548.2                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 11554.4684 </td>
   <td style="text-align:right;"> 5777.2342 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 26.094102 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2056.4245 </td>
   <td style="text-align:right;"> 2056.4245 </td>
   <td style="text-align:right;"> 1827 </td>
   <td style="text-align:right;"> 0.0519978 </td>
   <td style="text-align:right;"> 4.644139 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 812.3544 </td>
   <td style="text-align:right;"> 406.1772 </td>
   <td style="text-align:right;"> 461 </td>
   <td style="text-align:right;"> 0.3731020 </td>
   <td style="text-align:right;"> 1.834585 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 875.6763 </td>
   <td style="text-align:right;"> 145.9461 </td>
   <td style="text-align:right;"> 80 </td>
   <td style="text-align:right;"> 0.9375000 </td>
   <td style="text-align:right;"> 1.977589 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2817.7201 </td>
   <td style="text-align:right;"> 234.8100 </td>
   <td style="text-align:right;"> 566 </td>
   <td style="text-align:right;"> 0.8975265 </td>
   <td style="text-align:right;"> 6.363415 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2046.5389 </td>
   <td style="text-align:right;"> 341.0898 </td>
   <td style="text-align:right;"> 455 </td>
   <td style="text-align:right;"> 0.7560440 </td>
   <td style="text-align:right;"> 4.621813 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2735.3173 </td>
   <td style="text-align:right;"> 227.9431 </td>
   <td style="text-align:right;"> 2678 </td>
   <td style="text-align:right;"> 0.9458551 </td>
   <td style="text-align:right;"> 6.177320 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 21381.5000 </td>
   <td style="text-align:right;"> 548.2436 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 48.287037 </td>
  </tr>
</tbody>
</table>

```r
# Shannon

setkey(all_alpha_ord, shannon)
all_alpha_ord[, measure := as.numeric(as.factor(shannon))]
result <- aovp(update(formula, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS "
```

```r
summary(result)
```

```
# Component 1 :
#                    Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# Site                2  12291.8    6145.9 5000  < 2e-16 ***
# Storage             1   1077.6    1077.6  230  0.30435    
# Site:Storage        2   2320.4    1160.2 1380  0.07101 .  
# Scion               6    570.9      95.1  107  1.00000    
# Site:Scion         12   3082.1     256.8  385  0.90130    
# Storage:Scion       6   2730.6     455.1  396  0.69192    
# Site:Storage:Scion 12   5311.1     442.6 1397  0.50966    
# Residuals          39  16895.5     433.2                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 12291.7504 </td>
   <td style="text-align:right;"> 6145.87520 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 27.759147 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1077.5749 </td>
   <td style="text-align:right;"> 1077.57489 </td>
   <td style="text-align:right;"> 230 </td>
   <td style="text-align:right;"> 0.3043478 </td>
   <td style="text-align:right;"> 2.433548 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2320.4220 </td>
   <td style="text-align:right;"> 1160.21098 </td>
   <td style="text-align:right;"> 1380 </td>
   <td style="text-align:right;"> 0.0710145 </td>
   <td style="text-align:right;"> 5.240339 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 570.8675 </td>
   <td style="text-align:right;"> 95.14458 </td>
   <td style="text-align:right;"> 107 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 1.289222 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3082.1462 </td>
   <td style="text-align:right;"> 256.84552 </td>
   <td style="text-align:right;"> 385 </td>
   <td style="text-align:right;"> 0.9012987 </td>
   <td style="text-align:right;"> 6.960583 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2730.6068 </td>
   <td style="text-align:right;"> 455.10114 </td>
   <td style="text-align:right;"> 396 </td>
   <td style="text-align:right;"> 0.6919192 </td>
   <td style="text-align:right;"> 6.166682 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 5311.1323 </td>
   <td style="text-align:right;"> 442.59435 </td>
   <td style="text-align:right;"> 1397 </td>
   <td style="text-align:right;"> 0.5096636 </td>
   <td style="text-align:right;"> 11.994427 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 16895.5000 </td>
   <td style="text-align:right;"> 433.21795 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 38.156052 </td>
  </tr>
</tbody>
</table>

```r
# Simpson

setkey(all_alpha_ord, simpson)
all_alpha_ord[, measure := as.numeric(as.factor(simpson))]
result <- aovp(update(formula, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS "
```

```r
summary(result)
```

```
# Component 1 :
#                    Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# Site                2  12937.5    6468.8 5000  < 2e-16 ***
# Storage             1    764.2     764.2  339  0.23009    
# Site:Storage        2   2484.2    1242.1 2525  0.08079 .  
# Scion               6   1188.1     198.0  834  0.87650    
# Site:Scion         12   2027.8     169.0  637  0.98273    
# Storage:Scion       6   2529.6     421.6  486  0.67695    
# Site:Storage:Scion 12   5334.6     444.6  537  0.55493    
# Residuals          39  17014.0     436.3                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 12937.5044 </td>
   <td style="text-align:right;"> 6468.7522 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 29.217490 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 764.2007 </td>
   <td style="text-align:right;"> 764.2007 </td>
   <td style="text-align:right;"> 339 </td>
   <td style="text-align:right;"> 0.2300885 </td>
   <td style="text-align:right;"> 1.725837 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2484.2454 </td>
   <td style="text-align:right;"> 1242.1227 </td>
   <td style="text-align:right;"> 2525 </td>
   <td style="text-align:right;"> 0.0807921 </td>
   <td style="text-align:right;"> 5.610310 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1188.0542 </td>
   <td style="text-align:right;"> 198.0090 </td>
   <td style="text-align:right;"> 834 </td>
   <td style="text-align:right;"> 0.8764988 </td>
   <td style="text-align:right;"> 2.683049 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2027.7819 </td>
   <td style="text-align:right;"> 168.9818 </td>
   <td style="text-align:right;"> 637 </td>
   <td style="text-align:right;"> 0.9827316 </td>
   <td style="text-align:right;"> 4.579453 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2529.5940 </td>
   <td style="text-align:right;"> 421.5990 </td>
   <td style="text-align:right;"> 486 </td>
   <td style="text-align:right;"> 0.6769547 </td>
   <td style="text-align:right;"> 5.712724 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 5334.6194 </td>
   <td style="text-align:right;"> 444.5516 </td>
   <td style="text-align:right;"> 537 </td>
   <td style="text-align:right;"> 0.5549348 </td>
   <td style="text-align:right;"> 12.047469 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 17014.0000 </td>
   <td style="text-align:right;"> 436.2564 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 38.423668 </td>
  </tr>
</tbody>
</table>

#### Permutation based anova on α-diversity index ranks iwth canker lesion counts


```r
cat("Model formula: ", deparse(design_with_canker), "\n")
```

```
# Model formula:  y ~ Site * Storage * Scion * Cankers
```

```r
# Chao1

setkey(all_alpha_ord, S.chao1)
all_alpha_ord[, measure := as.numeric(as.factor(S.chao1))]
result <- aovp(update(design_with_canker, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS : numeric variables centered"
```

```r
summary(result)
```

```
# Component 1 :
#                            Df R Sum Sq R Mean Sq Iter Pr(Prob)  
# Site                        2   9641.2    4820.6 1933  0.08484 .
# Storage                     1   2245.6    2245.6  850  0.10588  
# Site:Storage                2   1732.7     866.4  510  0.38039  
# Scion                       6    866.6     144.4   77  0.87013  
# Site:Scion                 12   2746.9     228.9  233  0.73391  
# Storage:Scion               6   2160.7     360.1  330  0.65152  
# Site:Storage:Scion         11   2414.5     219.5  227  0.74449  
# Cankers                     1    407.0     407.0   51  0.86275  
# Site:Cankers                2     38.4      19.2  136  0.89706  
# Storage:Cankers             1   1732.2    1732.2  919  0.09902 .
# Site:Storage:Cankers        2   1183.3     591.7  505  0.47525  
# Scion:Cankers               6   2236.3     372.7  480  0.72083  
# Site:Scion:Cankers         11   6761.8     614.7  845  0.51953  
# Storage:Scion:Cankers       5   3278.4     655.7  853  0.55334  
# Site:Storage:Scion:Cankers  2    117.2      58.6   51  0.98039  
# Residuals                   2   1341.0     670.5                
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 9641.20503 </td>
   <td style="text-align:right;"> 4820.60252 </td>
   <td style="text-align:right;"> 1933 </td>
   <td style="text-align:right;"> 0.0848422 </td>
   <td style="text-align:right;"> 24.7821976 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2245.63863 </td>
   <td style="text-align:right;"> 2245.63863 </td>
   <td style="text-align:right;"> 850 </td>
   <td style="text-align:right;"> 0.1058824 </td>
   <td style="text-align:right;"> 5.7722930 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1732.74465 </td>
   <td style="text-align:right;"> 866.37232 </td>
   <td style="text-align:right;"> 510 </td>
   <td style="text-align:right;"> 0.3803922 </td>
   <td style="text-align:right;"> 4.4539267 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 866.59233 </td>
   <td style="text-align:right;"> 144.43205 </td>
   <td style="text-align:right;"> 77 </td>
   <td style="text-align:right;"> 0.8701299 </td>
   <td style="text-align:right;"> 2.2275288 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2746.94513 </td>
   <td style="text-align:right;"> 228.91209 </td>
   <td style="text-align:right;"> 233 </td>
   <td style="text-align:right;"> 0.7339056 </td>
   <td style="text-align:right;"> 7.0608743 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2160.66774 </td>
   <td style="text-align:right;"> 360.11129 </td>
   <td style="text-align:right;"> 330 </td>
   <td style="text-align:right;"> 0.6515152 </td>
   <td style="text-align:right;"> 5.5538799 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 2414.45992 </td>
   <td style="text-align:right;"> 219.49636 </td>
   <td style="text-align:right;"> 227 </td>
   <td style="text-align:right;"> 0.7444934 </td>
   <td style="text-align:right;"> 6.2062390 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 407.01518 </td>
   <td style="text-align:right;"> 407.01518 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 0.8627451 </td>
   <td style="text-align:right;"> 1.0462106 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 38.38784 </td>
   <td style="text-align:right;"> 19.19392 </td>
   <td style="text-align:right;"> 136 </td>
   <td style="text-align:right;"> 0.8970588 </td>
   <td style="text-align:right;"> 0.0986739 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1732.18898 </td>
   <td style="text-align:right;"> 1732.18898 </td>
   <td style="text-align:right;"> 919 </td>
   <td style="text-align:right;"> 0.0990207 </td>
   <td style="text-align:right;"> 4.4524983 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1183.30973 </td>
   <td style="text-align:right;"> 591.65487 </td>
   <td style="text-align:right;"> 505 </td>
   <td style="text-align:right;"> 0.4752475 </td>
   <td style="text-align:right;"> 3.0416338 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2236.26630 </td>
   <td style="text-align:right;"> 372.71105 </td>
   <td style="text-align:right;"> 480 </td>
   <td style="text-align:right;"> 0.7208333 </td>
   <td style="text-align:right;"> 5.7482019 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 6761.80566 </td>
   <td style="text-align:right;"> 614.70961 </td>
   <td style="text-align:right;"> 845 </td>
   <td style="text-align:right;"> 0.5195266 </td>
   <td style="text-align:right;"> 17.3808568 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3278.37298 </td>
   <td style="text-align:right;"> 655.67460 </td>
   <td style="text-align:right;"> 853 </td>
   <td style="text-align:right;"> 0.5533411 </td>
   <td style="text-align:right;"> 8.4268809 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 117.15333 </td>
   <td style="text-align:right;"> 58.57666 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 0.9803922 </td>
   <td style="text-align:right;"> 0.3011363 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1341.00000 </td>
   <td style="text-align:right;"> 670.50000 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 3.4469682 </td>
  </tr>
</tbody>
</table>

```r
# Shannon

setkey(all_alpha_ord, shannon)
all_alpha_ord[, measure := as.numeric(as.factor(shannon))]
result <- aovp(update(design_with_canker, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS : numeric variables centered"
```

```r
summary(result)
```

```
# Component 1 :
#                            Df R Sum Sq R Mean Sq Iter Pr(Prob)   
# Site                        2  14993.4    7496.7 5000  0.00520 **
# Storage                     1   1402.0    1402.0 3879  0.02526 * 
# Site:Storage                2   1038.5     519.3  261  0.33333   
# Scion                       6    529.2      88.2  235  0.68511   
# Site:Scion                 12   2999.7     250.0 1431  0.24598   
# Storage:Scion               6   2819.4     469.9 3034  0.12129   
# Site:Storage:Scion         11   2787.0     253.4 1068  0.25843   
# Cankers                     1     10.7      10.7   51  0.98039   
# Site:Cankers                2    514.8     257.4  928  0.15409   
# Storage:Cankers             1    718.3     718.3  311  0.24437   
# Site:Storage:Cankers        2    149.3      74.7  196  0.43878   
# Scion:Cankers               6   3676.2     612.7 1423  0.13985   
# Site:Scion:Cankers         11   4331.3     393.8 3041  0.14864   
# Storage:Scion:Cankers       5   3439.8     688.0 5000  0.07620 . 
# Site:Storage:Scion:Cankers  2   1409.9     705.0 1374  0.10189   
# Residuals                   2    102.5      51.2                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 14993.40799 </td>
   <td style="text-align:right;"> 7496.70399 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0052000 </td>
   <td style="text-align:right;"> 36.6387950 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1402.02677 </td>
   <td style="text-align:right;"> 1402.02677 </td>
   <td style="text-align:right;"> 3879 </td>
   <td style="text-align:right;"> 0.0252642 </td>
   <td style="text-align:right;"> 3.4260771 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1038.54266 </td>
   <td style="text-align:right;"> 519.27133 </td>
   <td style="text-align:right;"> 261 </td>
   <td style="text-align:right;"> 0.3333333 </td>
   <td style="text-align:right;"> 2.5378454 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 529.21633 </td>
   <td style="text-align:right;"> 88.20272 </td>
   <td style="text-align:right;"> 235 </td>
   <td style="text-align:right;"> 0.6851064 </td>
   <td style="text-align:right;"> 1.2932249 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2999.65362 </td>
   <td style="text-align:right;"> 249.97114 </td>
   <td style="text-align:right;"> 1431 </td>
   <td style="text-align:right;"> 0.2459818 </td>
   <td style="text-align:right;"> 7.3301343 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2819.38267 </td>
   <td style="text-align:right;"> 469.89711 </td>
   <td style="text-align:right;"> 3034 </td>
   <td style="text-align:right;"> 0.1212920 </td>
   <td style="text-align:right;"> 6.8896133 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 2786.98914 </td>
   <td style="text-align:right;"> 253.36265 </td>
   <td style="text-align:right;"> 1068 </td>
   <td style="text-align:right;"> 0.2584270 </td>
   <td style="text-align:right;"> 6.8104546 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 10.72218 </td>
   <td style="text-align:right;"> 10.72218 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 0.9803922 </td>
   <td style="text-align:right;"> 0.0262014 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 514.75261 </td>
   <td style="text-align:right;"> 257.37630 </td>
   <td style="text-align:right;"> 928 </td>
   <td style="text-align:right;"> 0.1540948 </td>
   <td style="text-align:right;"> 1.2578805 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 718.34815 </td>
   <td style="text-align:right;"> 718.34815 </td>
   <td style="text-align:right;"> 311 </td>
   <td style="text-align:right;"> 0.2443730 </td>
   <td style="text-align:right;"> 1.7553988 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 149.31807 </td>
   <td style="text-align:right;"> 74.65903 </td>
   <td style="text-align:right;"> 196 </td>
   <td style="text-align:right;"> 0.4387755 </td>
   <td style="text-align:right;"> 0.3648826 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 3676.23412 </td>
   <td style="text-align:right;"> 612.70569 </td>
   <td style="text-align:right;"> 1423 </td>
   <td style="text-align:right;"> 0.1398454 </td>
   <td style="text-align:right;"> 8.9834672 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 4331.34949 </td>
   <td style="text-align:right;"> 393.75904 </td>
   <td style="text-align:right;"> 3041 </td>
   <td style="text-align:right;"> 0.1486353 </td>
   <td style="text-align:right;"> 10.5843465 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 3439.82927 </td>
   <td style="text-align:right;"> 687.96585 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0762000 </td>
   <td style="text-align:right;"> 8.4057740 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1409.94613 </td>
   <td style="text-align:right;"> 704.97306 </td>
   <td style="text-align:right;"> 1374 </td>
   <td style="text-align:right;"> 0.1018923 </td>
   <td style="text-align:right;"> 3.4454293 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 102.50000 </td>
   <td style="text-align:right;"> 51.25000 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.2504752 </td>
  </tr>
</tbody>
</table>

```r
# Simpson

setkey(all_alpha_ord, simpson)
all_alpha_ord[, measure := as.numeric(as.factor(simpson))]
result <- aovp(update(design_with_canker, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS : numeric variables centered"
```

```r
summary(result)
```

```
# Component 1 :
#                            Df R Sum Sq R Mean Sq Iter Pr(Prob)  
# Site                        2  14576.4    7288.2 3836  0.03989 *
# Storage                     1    681.9     681.9  491  0.17108  
# Site:Storage                2   1252.8     626.4  296  0.26689  
# Scion                       6   1038.9     173.2  281  0.67616  
# Site:Scion                 12   2119.4     176.6  280  0.83929  
# Storage:Scion               6   2788.4     464.7  770  0.37922  
# Site:Storage:Scion         11   3985.5     362.3  796  0.47739  
# Cankers                     1      2.5       2.5   51  0.82353  
# Site:Cankers                2   1134.2     567.1  478  0.25314  
# Storage:Cankers             1    126.6     126.6   83  0.55422  
# Site:Storage:Cankers        2    921.6     460.8  212  0.37736  
# Scion:Cankers               6   5011.5     835.3  551  0.29038  
# Site:Scion:Cankers         11   3823.1     347.6  771  0.50195  
# Storage:Scion:Cankers       5   2957.8     591.6  843  0.28588  
# Site:Storage:Scion:Cankers  2   1607.3     803.6  873  0.19817  
# Residuals                   2    305.0     152.5                
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 14576.428600 </td>
   <td style="text-align:right;"> 7288.214300 </td>
   <td style="text-align:right;"> 3836 </td>
   <td style="text-align:right;"> 0.0398853 </td>
   <td style="text-align:right;"> 34.4328799 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 681.922331 </td>
   <td style="text-align:right;"> 681.922331 </td>
   <td style="text-align:right;"> 491 </td>
   <td style="text-align:right;"> 0.1710794 </td>
   <td style="text-align:right;"> 1.6108575 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1252.847015 </td>
   <td style="text-align:right;"> 626.423508 </td>
   <td style="text-align:right;"> 296 </td>
   <td style="text-align:right;"> 0.2668919 </td>
   <td style="text-align:right;"> 2.9595131 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1038.901626 </td>
   <td style="text-align:right;"> 173.150271 </td>
   <td style="text-align:right;"> 281 </td>
   <td style="text-align:right;"> 0.6761566 </td>
   <td style="text-align:right;"> 2.4541248 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2119.430716 </td>
   <td style="text-align:right;"> 176.619226 </td>
   <td style="text-align:right;"> 280 </td>
   <td style="text-align:right;"> 0.8392857 </td>
   <td style="text-align:right;"> 5.0065833 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2788.363440 </td>
   <td style="text-align:right;"> 464.727240 </td>
   <td style="text-align:right;"> 770 </td>
   <td style="text-align:right;"> 0.3792208 </td>
   <td style="text-align:right;"> 6.5867563 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 3985.482984 </td>
   <td style="text-align:right;"> 362.316635 </td>
   <td style="text-align:right;"> 796 </td>
   <td style="text-align:right;"> 0.4773869 </td>
   <td style="text-align:right;"> 9.4146283 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2.453294 </td>
   <td style="text-align:right;"> 2.453294 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 0.8235294 </td>
   <td style="text-align:right;"> 0.0057952 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1134.173594 </td>
   <td style="text-align:right;"> 567.086797 </td>
   <td style="text-align:right;"> 478 </td>
   <td style="text-align:right;"> 0.2531381 </td>
   <td style="text-align:right;"> 2.6791791 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 126.554905 </td>
   <td style="text-align:right;"> 126.554905 </td>
   <td style="text-align:right;"> 83 </td>
   <td style="text-align:right;"> 0.5542169 </td>
   <td style="text-align:right;"> 0.2989518 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 921.611262 </td>
   <td style="text-align:right;"> 460.805631 </td>
   <td style="text-align:right;"> 212 </td>
   <td style="text-align:right;"> 0.3773585 </td>
   <td style="text-align:right;"> 2.1770580 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 5011.520578 </td>
   <td style="text-align:right;"> 835.253430 </td>
   <td style="text-align:right;"> 551 </td>
   <td style="text-align:right;"> 0.2903811 </td>
   <td style="text-align:right;"> 11.8383653 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 3823.119329 </td>
   <td style="text-align:right;"> 347.556303 </td>
   <td style="text-align:right;"> 771 </td>
   <td style="text-align:right;"> 0.5019455 </td>
   <td style="text-align:right;"> 9.0310880 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 2957.787348 </td>
   <td style="text-align:right;"> 591.557470 </td>
   <td style="text-align:right;"> 843 </td>
   <td style="text-align:right;"> 0.2858837 </td>
   <td style="text-align:right;"> 6.9869746 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1607.279689 </td>
   <td style="text-align:right;"> 803.639845 </td>
   <td style="text-align:right;"> 873 </td>
   <td style="text-align:right;"> 0.1981672 </td>
   <td style="text-align:right;"> 3.7967646 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 305.000000 </td>
   <td style="text-align:right;"> 152.500000 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.7204802 </td>
  </tr>
</tbody>
</table>

## Beta diversity PCA/NMDS

### PCA


```r
# Number of PCs to include
n_pcs <- 10

# Perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# To get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
fun_pca <- t(data.frame(t(mypca$x) * mypca$percentVar))

formula = FULL_DESIGN
```

#### Percent variation in first 10 PCs 


```r
# Cumulative percentage of variance explained
pca_cum_var <- data.frame(
  cumulative = cumsum(mypca$percentVar * 100),
  no = 1:length(mypca$percentVar)
)

# Plot cumulative percentage of variance explained
fun_cum_pca <- ggline(
  pca_cum_var, x = "no", y = "cumulative", plot_type = "l",
  xlab = "Number of PCs", ylab = "Cumulative % variance explained",
  title = "Fungi: cumulative % variance explained by PCs"
)
ggsave(filename = "fun_cum_pca.png", plot = fun_cum_pca, path = "figures/",)
fun_cum_pca
```

![](root_endophytes_files/figure-html/FUN PCA var-1.png)<!-- -->

```r
pca_var <- data.frame(
  PC = paste0("PC", 1:n_pcs),
  perc_var = round(mypca$percentVar[1:n_pcs] * 100, 1)
)

pca_var %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:right;"> perc_var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:right;"> 27.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:right;"> 21.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:right;"> 8.4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:right;"> 4.2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:right;"> 3.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:right;"> 2.6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:right;"> 1.9 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:right;"> 1.8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:right;"> 1.7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:right;"> 1.6 </td>
  </tr>
</tbody>
</table>

#### ANOVA of first 10 PCs


```r
pca_summary <- apply(
  mypca$x[, 1:n_pcs], 2, 
  function(x){
    summary(aov(update(formula, x ~ .), data = as.data.frame(cbind(x, colData(dds)))))
  }
)

pca_summary
```

```
# $PC1
#                    Df Sum Sq Mean Sq F value  Pr(>F)    
# Site                2  71073   35537 238.166 < 2e-16 ***
# Storage             1     31      31   0.208 0.65102    
# Scion               6    871     145   0.973 0.45635    
# Site:Storage        2   2232    1116   7.480 0.00178 ** 
# Site:Scion         12   1299     108   0.726 0.71791    
# Storage:Scion       6    696     116   0.777 0.59290    
# Site:Storage:Scion 12   1781     148   0.994 0.47159    
# Residuals          39   5819     149                    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC2
#                    Df Sum Sq Mean Sq F value   Pr(>F)    
# Site                2  41657   20829  62.213 7.35e-13 ***
# Storage             1    411     411   1.228  0.27459    
# Scion               6    125      21   0.062  0.99892    
# Site:Storage        2   5403    2701   8.069  0.00117 ** 
# Site:Scion         12   2269     189   0.565  0.85649    
# Storage:Scion       6    416      69   0.207  0.97246    
# Site:Storage:Scion 12   2080     173   0.518  0.89032    
# Residuals          39  13057     335                     
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC3
#                    Df Sum Sq Mean Sq F value Pr(>F)  
# Site                2   2091  1045.4   3.075 0.0575 .
# Storage             1    250   250.4   0.737 0.3960  
# Scion               6    976   162.7   0.479 0.8202  
# Site:Storage        2   1598   798.8   2.350 0.1087  
# Site:Scion         12   3182   265.2   0.780 0.6668  
# Storage:Scion       6   1399   233.2   0.686 0.6619  
# Site:Storage:Scion 12   3176   264.7   0.779 0.6682  
# Residuals          39  13257   339.9                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC4
#                    Df Sum Sq Mean Sq F value  Pr(>F)   
# Site                2   1692   846.2   5.486 0.00795 **
# Storage             1    156   155.8   1.010 0.32103   
# Scion               6    166    27.7   0.179 0.98087   
# Site:Storage        2   1216   608.2   3.943 0.02757 * 
# Site:Scion         12   1399   116.6   0.756 0.68963   
# Storage:Scion       6    637   106.1   0.688 0.66035   
# Site:Storage:Scion 12   1578   131.5   0.853 0.59855   
# Residuals          39   6016   154.2                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC5
#                    Df Sum Sq Mean Sq F value Pr(>F)
# Site                2    238   119.1   0.739  0.484
# Storage             1    457   456.7   2.835  0.100
# Scion               6    326    54.4   0.338  0.913
# Site:Storage        2    313   156.3   0.970  0.388
# Site:Scion         12   1292   107.7   0.669  0.770
# Storage:Scion       6    534    88.9   0.552  0.765
# Site:Storage:Scion 12    680    56.7   0.352  0.973
# Residuals          39   6281   161.1               
# 
# $PC6
#                    Df Sum Sq Mean Sq F value Pr(>F)
# Site                2    141   70.29   0.602  0.553
# Storage             1    279  278.51   2.384  0.131
# Scion               6    225   37.51   0.321  0.922
# Site:Storage        2    262  130.85   1.120  0.337
# Site:Scion         12    428   35.65   0.305  0.985
# Storage:Scion       6    582   96.92   0.830  0.554
# Site:Storage:Scion 12   1540  128.37   1.099  0.388
# Residuals          39   4557  116.85               
# 
# $PC7
#                    Df Sum Sq Mean Sq F value Pr(>F)
# Site                2   44.9   22.47   0.283  0.755
# Storage             1   16.0   16.02   0.202  0.655
# Scion               6  689.5  114.92   1.450  0.221
# Site:Storage        2   13.3    6.67   0.084  0.919
# Site:Scion         12 1060.5   88.37   1.115  0.376
# Storage:Scion       6  265.0   44.16   0.557  0.761
# Site:Storage:Scion 12  838.4   69.87   0.881  0.572
# Residuals          39 3091.5   79.27               
# 
# $PC8
#                    Df Sum Sq Mean Sq F value Pr(>F)   
# Site                2    1.1     0.6   0.008 0.9919   
# Storage             1  425.5   425.5   6.057 0.0184 * 
# Scion               6  147.8    24.6   0.351 0.9051   
# Site:Storage        2  769.7   384.9   5.478 0.0080 **
# Site:Scion         12  613.2    51.1   0.727 0.7163   
# Storage:Scion       6  594.6    99.1   1.411 0.2352   
# Site:Storage:Scion 12  401.4    33.5   0.476 0.9166   
# Residuals          39 2739.9    70.3                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC9
#                    Df Sum Sq Mean Sq F value   Pr(>F)    
# Site                2   74.9    37.4   0.648   0.5287    
# Storage             1  210.8   210.8   3.647   0.0635 .  
# Scion               6  246.7    41.1   0.711   0.6425    
# Site:Storage        2 1500.6   750.3  12.982 4.77e-05 ***
# Site:Scion         12  418.6    34.9   0.604   0.8258    
# Storage:Scion       6  238.5    39.8   0.688   0.6606    
# Site:Storage:Scion 12  206.8    17.2   0.298   0.9860    
# Residuals          39 2254.0    57.8                     
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC10
#                    Df Sum Sq Mean Sq F value Pr(>F)
# Site                2   17.2    8.62   0.107  0.898
# Storage             1   13.0   12.98   0.162  0.690
# Scion               6  279.3   46.55   0.580  0.744
# Site:Storage        2  112.7   56.35   0.702  0.502
# Site:Scion         12  544.6   45.38   0.566  0.856
# Storage:Scion       6   80.8   13.46   0.168  0.984
# Site:Storage:Scion 12  646.4   53.87   0.671  0.767
# Residuals          39 3129.1   80.23
```

#### Percent variation in first 10 PCs for each factor


```r
# Extract PC scores as a list of dataframes
pcas <- lapply(pca_summary, function(i) data.frame(unclass(i)))

# Merge into single dataframe
pcs_factors_tidy <- lapply(
  names(pcas),
  function(name) {
    pcas[[name]] %>%
    mutate(
      PC = name, #substring(name, 3),
      Factor = gsub(" ", "", rownames(pcas[[name]])),
      var = Sum.Sq / sum(pcas[[name]]$Sum.Sq) * 100,
      pc_var = subset(pca_var, PC == name)$"perc_var",
      total_var = var * pc_var / 100,
      sig = case_when(
        is.na(Pr..F.) ~ "",
        Pr..F. < 0.001 ~ "***",
        Pr..F. < 0.01 ~ "**",
        Pr..F. < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variance = ifelse(
        total_var < 0.01, paste0("<0.01", sig),
        paste0(round(total_var, 2), sig)
      )
    )
  }
) %>% bind_rows() %>% data.table()

# Order PCs and factors
pcs_factors_tidy$PC <- factor(pcs_factors_tidy$PC, levels = paste0("PC", 1:n_pcs))
pcs_factors_tidy$Factor <- factor(pcs_factors_tidy$Factor, levels = unique(pcs_factors_tidy$Factor))

# Significant factors
pcs_factors_tidy[
  Pr..F. < 0.05, 
  c("PC", "Factor", "Df", "F.value", "Pr..F.", "var", "pc_var", "total_var")
] %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:left;"> Factor </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> F.value </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> var </th>
   <th style="text-align:right;"> pc_var </th>
   <th style="text-align:right;"> total_var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 238.165810 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 84.810927 </td>
   <td style="text-align:right;"> 27.1 </td>
   <td style="text-align:right;"> 22.9837611 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 7.479562 </td>
   <td style="text-align:right;"> 0.0017803 </td>
   <td style="text-align:right;"> 2.663475 </td>
   <td style="text-align:right;"> 27.1 </td>
   <td style="text-align:right;"> 0.7218016 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 62.212524 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 63.678672 </td>
   <td style="text-align:right;"> 21.2 </td>
   <td style="text-align:right;"> 13.4998786 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8.068552 </td>
   <td style="text-align:right;"> 0.0011684 </td>
   <td style="text-align:right;"> 8.258702 </td>
   <td style="text-align:right;"> 21.2 </td>
   <td style="text-align:right;"> 1.7508448 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 5.485916 </td>
   <td style="text-align:right;"> 0.0079546 </td>
   <td style="text-align:right;"> 13.159308 </td>
   <td style="text-align:right;"> 4.2 </td>
   <td style="text-align:right;"> 0.5526910 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3.942945 </td>
   <td style="text-align:right;"> 0.0275701 </td>
   <td style="text-align:right;"> 9.458115 </td>
   <td style="text-align:right;"> 4.2 </td>
   <td style="text-align:right;"> 0.3972408 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 6.057208 </td>
   <td style="text-align:right;"> 0.0183793 </td>
   <td style="text-align:right;"> 7.474228 </td>
   <td style="text-align:right;"> 1.8 </td>
   <td style="text-align:right;"> 0.1345361 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 5.478399 </td>
   <td style="text-align:right;"> 0.0080014 </td>
   <td style="text-align:right;"> 13.520024 </td>
   <td style="text-align:right;"> 1.8 </td>
   <td style="text-align:right;"> 0.2433604 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 12.982328 </td>
   <td style="text-align:right;"> 0.0000477 </td>
   <td style="text-align:right;"> 29.133773 </td>
   <td style="text-align:right;"> 1.7 </td>
   <td style="text-align:right;"> 0.4952741 </td>
  </tr>
</tbody>
</table>

```r
# Table with factors as columns and PCs as rows
# pcs_factors <- dcast(pcs_factors_tidy, PC ~ Factor, value.var = "variance")
pcs_factors <- pcs_factors_tidy %>%
  select(PC, pc_var, Factor, variance) %>%
  spread(key = Factor, value = variance)

# Add a row with sum of each column by removing '<' and '*' and converting to numeric
pc_sums <- pcs_factors %>%
  mutate_all(~ gsub("[<|*]", "", .)) %>%
  mutate_all(as.numeric) %>%
  colSums()

pcs_factors %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:right;"> pc_var </th>
   <th style="text-align:left;"> Site </th>
   <th style="text-align:left;"> Storage </th>
   <th style="text-align:left;"> Scion </th>
   <th style="text-align:left;"> Site:Storage </th>
   <th style="text-align:left;"> Site:Scion </th>
   <th style="text-align:left;"> Storage:Scion </th>
   <th style="text-align:left;"> Site:Storage:Scion </th>
   <th style="text-align:left;"> Residuals </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:right;"> 27.1 </td>
   <td style="text-align:left;"> 22.98*** </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> 0.28 </td>
   <td style="text-align:left;"> 0.72** </td>
   <td style="text-align:left;"> 0.42 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> 0.58 </td>
   <td style="text-align:left;"> 1.88 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:right;"> 21.2 </td>
   <td style="text-align:left;"> 13.5*** </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 1.75** </td>
   <td style="text-align:left;"> 0.74 </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> 0.67 </td>
   <td style="text-align:left;"> 4.23 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:right;"> 8.4 </td>
   <td style="text-align:left;"> 0.68 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.32 </td>
   <td style="text-align:left;"> 0.52 </td>
   <td style="text-align:left;"> 1.03 </td>
   <td style="text-align:left;"> 0.45 </td>
   <td style="text-align:left;"> 1.03 </td>
   <td style="text-align:left;"> 4.29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:right;"> 4.2 </td>
   <td style="text-align:left;"> 0.55** </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.4* </td>
   <td style="text-align:left;"> 0.46 </td>
   <td style="text-align:left;"> 0.21 </td>
   <td style="text-align:left;"> 0.52 </td>
   <td style="text-align:left;"> 1.96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:right;"> 3.3 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.11 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.42 </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> 2.05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:right;"> 2.6 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> 0.5 </td>
   <td style="text-align:left;"> 1.48 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:right;"> 1.9 </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.33 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.26 </td>
   <td style="text-align:left;"> 0.98 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:right;"> 1.8 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.13* </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.24** </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> 0.87 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:right;"> 1.7 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.5*** </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.74 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:right;"> 1.6 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> 0.21 </td>
   <td style="text-align:left;"> 1.04 </td>
  </tr>
</tbody>
</table>

```r
# save_kable(kable(pcs_factors, format = "pandoc"), "tables/FUN_pca_factors.docx")
```

#### PCA plot


```r
fun_pca_plot <- plotOrd(
  fun_pca,
  colData(dds),
  design = "Site",
  shapes = "Storage",
  axes = c(1, 2),
  cbPalette = T,
  alpha = 0.75,
) # + facet_wrap(~facet) 
#   geom_line(aes(group=facet),alpha=0.25,linetype=3,colour="#000000") + 
#   theme(text = element_text(size=14))

ggsave(filename = "fun_pca_plot.png", plot = fun_pca_plot, path = "figures/")

fun_pca_plot
```

![](root_endophytes_files/figure-html/FUN PCA plot-1.png)<!-- -->

#### PCA sum of squares (% var)


```r
sum_squares <- apply(mypca$x, 2 ,function(x) 
  summary(aov(update(formula, x ~ .), data = cbind(x, colData(dds))))[[1]][2]
)
sum_squares <- do.call(cbind, sum_squares)
x <- t(apply(sum_squares, 2, prop.table))
perVar <- x * mypca$percentVar
#colSums(perVar)
round(colSums(perVar) / sum(colSums(perVar)) * 100, 3)
```

```
# Site               Storage            Scion              Site:Storage       
#             38.022              1.000              3.083              4.878 
# Site:Scion         Storage:Scion      Site:Storage:Scion Residuals          
#              7.803              4.057              8.281             32.876
```

### PCA with canker lesion counts


```r
pca_summary <- apply(
  mypca$x[, 1:n_pcs], 2, 
  function(x){
    summary(aov(update(design_with_canker, x ~ .), data = as.data.frame(cbind(x, colData(dds)))))
  }
)

pca_summary
```

```
# $PC1
#                            Df Sum Sq Mean Sq F value  Pr(>F)   
# Site                        2  64511   32256 120.622 0.00822 **
# Storage                     1    199     199   0.743 0.47948   
# Scion                       6   1226     204   0.764 0.66243   
# Cankers                     1    470     470   1.758 0.31602   
# Site:Storage                2   1290     645   2.413 0.29303   
# Site:Scion                 12   1353     113   0.422 0.86439   
# Storage:Scion               6    990     165   0.617 0.72629   
# Site:Cankers                2     45      22   0.084 0.92281   
# Storage:Cankers             1    113     113   0.423 0.58200   
# Scion:Cankers               6    272      45   0.169 0.96170   
# Site:Storage:Scion         11   1884     171   0.641 0.74696   
# Site:Storage:Cankers        2    162      81   0.304 0.76715   
# Site:Scion:Cankers         11    891      81   0.303 0.92475   
# Storage:Scion:Cankers       5    532     106   0.398 0.82425   
# Site:Storage:Scion:Cankers  2   1341     670   2.507 0.28515   
# Residuals                   2    535     267                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 8 observations deleted due to missingness
# 
# $PC2
#                            Df Sum Sq Mean Sq F value Pr(>F)  
# Site                        2  43789   21894  55.915 0.0176 *
# Storage                     1   1322    1322   3.376 0.2075  
# Scion                       6    286      48   0.122 0.9809  
# Cankers                     1    451     451   1.153 0.3953  
# Site:Storage                2   3773    1886   4.818 0.1719  
# Site:Scion                 12   1240     103   0.264 0.9470  
# Storage:Scion               6    725     121   0.309 0.8888  
# Site:Cankers                2    492     246   0.628 0.6142  
# Storage:Cankers             1     73      73   0.188 0.7072  
# Scion:Cankers               6   1781     297   0.758 0.6649  
# Site:Storage:Scion         11   1553     141   0.361 0.8941  
# Site:Storage:Cankers        2    235     117   0.300 0.7694  
# Site:Scion:Cankers         11   2549     232   0.592 0.7709  
# Storage:Scion:Cankers       5   2955     591   1.509 0.4445  
# Site:Storage:Scion:Cankers  2    749     374   0.956 0.5112  
# Residuals                   2    783     392                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 8 observations deleted due to missingness
# 
# $PC3
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2   2531  1265.4   1.547  0.393
# Storage                     1     65    64.6   0.079  0.805
# Scion                       6    868   144.7   0.177  0.958
# Cankers                     1     28    27.7   0.034  0.871
# Site:Storage                2   2519  1259.3   1.540  0.394
# Site:Scion                 12   2398   199.8   0.244  0.956
# Storage:Scion               6    885   147.5   0.180  0.957
# Site:Cankers                2    361   180.3   0.220  0.819
# Storage:Cankers             1    652   652.1   0.797  0.466
# Scion:Cankers               6   1681   280.1   0.342  0.870
# Site:Storage:Scion         11   5346   486.0   0.594  0.770
# Site:Storage:Cankers        2    534   267.1   0.326  0.754
# Site:Scion:Cankers         11   2589   235.3   0.288  0.932
# Storage:Scion:Cankers       5   1012   202.4   0.247  0.910
# Site:Storage:Scion:Cankers  2    328   163.8   0.200  0.833
# Residuals                   2   1636   818.0               
# 8 observations deleted due to missingness
# 
# $PC4
#                            Df Sum Sq Mean Sq F value Pr(>F)  
# Site                        2  671.3   335.7  16.430 0.0574 .
# Storage                     1    0.6     0.6   0.028 0.8828  
# Scion                       6  335.8    56.0   2.740 0.2914  
# Cankers                     1    2.5     2.5   0.124 0.7583  
# Site:Storage                2  355.0   177.5   8.688 0.1032  
# Site:Scion                 12  826.8    68.9   3.372 0.2513  
# Storage:Scion               6  423.5    70.6   3.455 0.2414  
# Site:Cankers                2  254.5   127.2   6.228 0.1384  
# Storage:Cankers             1  204.8   204.8  10.024 0.0869 .
# Scion:Cankers               6  793.9   132.3   6.476 0.1398  
# Site:Storage:Scion         11  864.2    78.6   3.845 0.2244  
# Site:Storage:Cankers        2  336.1   168.0   8.225 0.1084  
# Site:Scion:Cankers         11 1176.3   106.9   5.234 0.1712  
# Storage:Scion:Cankers       5  142.5    28.5   1.395 0.4676  
# Site:Storage:Scion:Cankers  2   61.0    30.5   1.493 0.4012  
# Residuals                   2   40.9    20.4                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 8 observations deleted due to missingness
# 
# $PC5
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2  628.5  314.27   5.472  0.155
# Storage                     1   89.9   89.86   1.565  0.337
# Scion                       6   64.4   10.74   0.187  0.954
# Cankers                     1   28.4   28.37   0.494  0.555
# Site:Storage                2   69.3   34.66   0.603  0.624
# Site:Scion                 12  400.2   33.35   0.581  0.780
# Storage:Scion               6  292.8   48.79   0.850  0.630
# Site:Cankers                2  215.4  107.72   1.876  0.348
# Storage:Cankers             1   43.9   43.86   0.764  0.474
# Scion:Cankers               6  177.0   29.51   0.514  0.777
# Site:Storage:Scion         11  378.5   34.41   0.599  0.767
# Site:Storage:Cankers        2    8.1    4.04   0.070  0.934
# Site:Scion:Cankers         11  725.0   65.90   1.148  0.555
# Storage:Scion:Cankers       5  528.8  105.77   1.842  0.388
# Site:Storage:Scion:Cankers  2   87.3   43.65   0.760  0.568
# Residuals                   2  114.9   57.43               
# 8 observations deleted due to missingness
# 
# $PC6
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2  175.4    87.7   0.236  0.809
# Storage                     1   56.7    56.7   0.153  0.734
# Scion                       6  217.8    36.3   0.098  0.988
# Cankers                     1  186.7   186.7   0.503  0.552
# Site:Storage                2  463.1   231.6   0.623  0.616
# Site:Scion                 12  353.4    29.5   0.079  0.999
# Storage:Scion               6  643.3   107.2   0.289  0.900
# Site:Cankers                2  579.0   289.5   0.779  0.562
# Storage:Cankers             1    7.4     7.4   0.020  0.900
# Scion:Cankers               6  268.3    44.7   0.120  0.981
# Site:Storage:Scion         11  841.3    76.5   0.206  0.969
# Site:Storage:Cankers        2  162.7    81.3   0.219  0.820
# Site:Scion:Cankers         11 1487.4   135.2   0.364  0.892
# Storage:Scion:Cankers       5  251.0    50.2   0.135  0.968
# Site:Storage:Scion:Cankers  2  268.2   134.1   0.361  0.735
# Residuals                   2  742.9   371.5               
# 8 observations deleted due to missingness
# 
# $PC7
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2   26.4   13.19   0.099  0.910
# Storage                     1    5.3    5.26   0.039  0.861
# Scion                       6  727.9  121.32   0.909  0.608
# Cankers                     1    3.2    3.21   0.024  0.891
# Site:Storage                2    8.3    4.15   0.031  0.970
# Site:Scion                 12 1187.2   98.93   0.742  0.704
# Storage:Scion               6  138.2   23.04   0.173  0.960
# Site:Cankers                2    8.5    4.27   0.032  0.969
# Storage:Cankers             1    9.9    9.86   0.074  0.811
# Scion:Cankers               6  185.4   30.89   0.232  0.931
# Site:Storage:Scion         11  958.7   87.15   0.653  0.741
# Site:Storage:Cankers        2  133.5   66.74   0.500  0.667
# Site:Scion:Cankers         11 1856.1  168.73   1.265  0.522
# Storage:Scion:Cankers       5  272.1   54.43   0.408  0.819
# Site:Storage:Scion:Cankers  2   89.3   44.64   0.335  0.749
# Residuals                   2  266.8  133.41               
# 8 observations deleted due to missingness
# 
# $PC8
#                            Df Sum Sq Mean Sq F value Pr(>F)  
# Site                        2   11.0     5.5   0.307 0.7651  
# Storage                     1  574.2   574.2  31.914 0.0299 *
# Scion                       6  155.4    25.9   1.439 0.4647  
# Cankers                     1   67.8    67.8   3.770 0.1917  
# Site:Storage                2  595.1   297.6  16.538 0.0570 .
# Site:Scion                 12  555.2    46.3   2.571 0.3140  
# Storage:Scion               6  579.9    96.7   5.372 0.1652  
# Site:Cankers                2    6.6     3.3   0.183 0.8455  
# Storage:Cankers             1   14.1    14.1   0.781 0.4700  
# Scion:Cankers               6  627.4   104.6   5.812 0.1541  
# Site:Storage:Scion         11  326.5    29.7   1.650 0.4373  
# Site:Storage:Cankers        2  488.9   244.4  13.585 0.0686 .
# Site:Scion:Cankers         11 1354.2   123.1   6.842 0.1343  
# Storage:Scion:Cankers       5  160.5    32.1   1.784 0.3970  
# Site:Storage:Scion:Cankers  2   71.7    35.8   1.992 0.3342  
# Residuals                   2   36.0    18.0                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 8 observations deleted due to missingness
# 
# $PC9
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2   38.0    19.0   0.206  0.829
# Storage                     1  213.1   213.1   2.305  0.268
# Scion                       6  352.9    58.8   0.636  0.717
# Cankers                     1    7.0     7.0   0.076  0.809
# Site:Storage                2 1474.6   737.3   7.974  0.111
# Site:Scion                 12  421.8    35.1   0.380  0.887
# Storage:Scion               6  272.5    45.4   0.491  0.789
# Site:Cankers                2   18.1     9.1   0.098  0.911
# Storage:Cankers             1   75.0    75.0   0.811  0.463
# Scion:Cankers               6  318.9    53.2   0.575  0.746
# Site:Storage:Scion         11  419.8    38.2   0.413  0.866
# Site:Storage:Cankers        2  113.9    56.9   0.616  0.619
# Site:Scion:Cankers         11  368.6    33.5   0.362  0.893
# Storage:Scion:Cankers       5  414.5    82.9   0.897  0.602
# Site:Storage:Scion:Cankers  2  197.2    98.6   1.066  0.484
# Residuals                   2  184.9    92.5               
# 8 observations deleted due to missingness
# 
# $PC10
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2    2.9     1.4   0.004  0.996
# Storage                     1   47.5    47.5   0.121  0.762
# Scion                       6  339.1    56.5   0.144  0.973
# Cankers                     1   13.4    13.4   0.034  0.871
# Site:Storage                2  206.0   103.0   0.262  0.793
# Site:Scion                 12  568.8    47.4   0.120  0.995
# Storage:Scion               6  170.3    28.4   0.072  0.994
# Site:Cankers                2   27.2    13.6   0.035  0.967
# Storage:Cankers             1    1.8     1.8   0.005  0.952
# Scion:Cankers               6  642.0   107.0   0.272  0.909
# Site:Storage:Scion         11  713.8    64.9   0.165  0.983
# Site:Storage:Cankers        2  204.6   102.3   0.260  0.794
# Site:Scion:Cankers         11  641.4    58.3   0.148  0.988
# Storage:Scion:Cankers       5  196.4    39.3   0.100  0.982
# Site:Storage:Scion:Cankers  2  104.2    52.1   0.132  0.883
# Residuals                   2  787.6   393.8               
# 8 observations deleted due to missingness
```

```r
# Extract PC scores as a list of dataframes
pcas <- lapply(pca_summary, function(i) data.frame(unclass(i)))

# Merge into single dataframe
pcs_factors_tidy <- lapply(
  names(pcas),
  function(name) {
    pcas[[name]] %>%
    mutate(
      PC = name, #substring(name, 3),
      Factor = gsub(" ", "", rownames(pcas[[name]])),
      var = Sum.Sq / sum(pcas[[name]]$Sum.Sq) * 100,
      pc_var = subset(pca_var, PC == name)$"perc_var",
      total_var = var * pc_var / 100,
      sig = case_when(
        is.na(Pr..F.) ~ "",
        Pr..F. < 0.001 ~ "***",
        Pr..F. < 0.01 ~ "**",
        Pr..F. < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variance = ifelse(
        total_var < 0.01, paste0("<0.01", sig),
        paste0(round(total_var, 2), sig)
      )
    )
  }
) %>% bind_rows() %>% data.table()

# Order PCs and factors
pcs_factors_tidy$PC <- factor(pcs_factors_tidy$PC, levels = paste0("PC", 1:n_pcs))
pcs_factors_tidy$Factor <- factor(pcs_factors_tidy$Factor, levels = unique(pcs_factors_tidy$Factor))

# Significant factors
pcs_factors_tidy[
  Pr..F. < 0.05, 
  c("PC", "Factor", "Df", "F.value", "Pr..F.", "var", "pc_var", "total_var")
] %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:left;"> Factor </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> F.value </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> var </th>
   <th style="text-align:right;"> pc_var </th>
   <th style="text-align:right;"> total_var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 120.62182 </td>
   <td style="text-align:right;"> 0.0082222 </td>
   <td style="text-align:right;"> 85.09021 </td>
   <td style="text-align:right;"> 27.1 </td>
   <td style="text-align:right;"> 23.0594461 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 55.91478 </td>
   <td style="text-align:right;"> 0.0175701 </td>
   <td style="text-align:right;"> 69.77654 </td>
   <td style="text-align:right;"> 21.2 </td>
   <td style="text-align:right;"> 14.7926268 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 31.91434 </td>
   <td style="text-align:right;"> 0.0299341 </td>
   <td style="text-align:right;"> 10.20915 </td>
   <td style="text-align:right;"> 1.8 </td>
   <td style="text-align:right;"> 0.1837646 </td>
  </tr>
</tbody>
</table>

```r
# Table with factors as columns and PCs as rows
# pcs_factors <- dcast(pcs_factors_tidy, PC ~ Factor, value.var = "variance")
pcs_factors <- pcs_factors_tidy %>%
  select(PC, pc_var, Factor, variance) %>%
  spread(key = Factor, value = variance)

# Add a row with sum of each column by removing '<' and '*' and converting to numeric
pc_sums <- pcs_factors %>%
  mutate_all(~ gsub("[<|*]", "", .)) %>%
  mutate_all(as.numeric) %>%
  colSums()

pcs_factors %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:right;"> pc_var </th>
   <th style="text-align:left;"> Site </th>
   <th style="text-align:left;"> Storage </th>
   <th style="text-align:left;"> Scion </th>
   <th style="text-align:left;"> Cankers </th>
   <th style="text-align:left;"> Site:Storage </th>
   <th style="text-align:left;"> Site:Scion </th>
   <th style="text-align:left;"> Storage:Scion </th>
   <th style="text-align:left;"> Site:Cankers </th>
   <th style="text-align:left;"> Storage:Cankers </th>
   <th style="text-align:left;"> Scion:Cankers </th>
   <th style="text-align:left;"> Site:Storage:Scion </th>
   <th style="text-align:left;"> Site:Storage:Cankers </th>
   <th style="text-align:left;"> Site:Scion:Cankers </th>
   <th style="text-align:left;"> Storage:Scion:Cankers </th>
   <th style="text-align:left;"> Site:Storage:Scion:Cankers </th>
   <th style="text-align:left;"> Residuals </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:right;"> 27.1 </td>
   <td style="text-align:left;"> 23.06** </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.44 </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> 0.46 </td>
   <td style="text-align:left;"> 0.48 </td>
   <td style="text-align:left;"> 0.35 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.67 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.32 </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> 0.48 </td>
   <td style="text-align:left;"> 0.19 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:right;"> 21.2 </td>
   <td style="text-align:left;"> 14.79* </td>
   <td style="text-align:left;"> 0.45 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 1.27 </td>
   <td style="text-align:left;"> 0.42 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.6 </td>
   <td style="text-align:left;"> 0.52 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.86 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> 0.26 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:right;"> 8.4 </td>
   <td style="text-align:left;"> 0.91 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.31 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.9 </td>
   <td style="text-align:left;"> 0.86 </td>
   <td style="text-align:left;"> 0.32 </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> 0.6 </td>
   <td style="text-align:left;"> 1.92 </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> 0.93 </td>
   <td style="text-align:left;"> 0.36 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.59 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:right;"> 4.2 </td>
   <td style="text-align:left;"> 0.43 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> 0.54 </td>
   <td style="text-align:left;"> 0.27 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> 0.51 </td>
   <td style="text-align:left;"> 0.56 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> 0.76 </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.03 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:right;"> 3.3 </td>
   <td style="text-align:left;"> 0.54 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.34 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.32 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.62 </td>
   <td style="text-align:left;"> 0.45 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:right;"> 2.6 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.33 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.58 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.29 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:right;"> 1.9 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.24 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.38 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.31 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.6 </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> 0.09 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:right;"> 1.8 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.18* </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.43 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.01 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:right;"> 1.7 </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.51 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> 0.11 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:right;"> 1.6 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> 0.24 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.27 </td>
  </tr>
</tbody>
</table>

### ADONIS


```r
# Calculate Bray-Curtis distance matrix
vg <- vegdist(t(counts(dds, normalize = T)), method = "bray")

formula <- update(FULL_DESIGN, vg ~ .)

set.seed(sum(utf8ToInt("Hamish McLean")))
result <- adonis2(formula, colData(dds), permutations = 1000)
result
```

```
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = formula, data = colData(dds), permutations = 1000)
#          Df SumOfSqs     R2      F   Pr(>F)    
# Model    41  13.9181 0.6707 1.9374 0.000999 ***
# Residual 39   6.8334 0.3293                    
# Total    80  20.7514 1.0000                    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% data.frame()
df$Perc.Var <- df$SumOfSqs / df["Total", "SumOfSqs"] * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> SumOfSqs </th>
   <th style="text-align:right;"> R2 </th>
   <th style="text-align:right;"> F </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Model </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 13.918056 </td>
   <td style="text-align:right;"> 0.6707038 </td>
   <td style="text-align:right;"> 1.937424 </td>
   <td style="text-align:right;"> 0.000999 </td>
   <td style="text-align:right;"> 67.07038 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residual </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 6.833365 </td>
   <td style="text-align:right;"> 0.3292962 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 32.92962 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Total </td>
   <td style="text-align:right;"> 80 </td>
   <td style="text-align:right;"> 20.751421 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 100.00000 </td>
  </tr>
</tbody>
</table>

### NMDS ordination


```r
set.seed(sum(utf8ToInt("Hamish McLean")))
ord <- metaMDS(vg,trace=0) 
#sratmax=20000,maxit=20000,try = 177, trymax = 177

fun_nmds <- scores(ord)

fun_nmds_plot <- plotOrd(
  fun_nmds, colData(dds), 
  design = "Site", 
  shape = "Storage", 
  alpha = 0.75, cbPalette = T
) #+ theme(text = element_text(size = 14))

ggsave(filename = "fun_nmds_plot.png", plot = fun_nmds_plot, path = "figures/")

fun_nmds_plot
```

![](root_endophytes_files/figure-html/FUN NMDS-1.png)<!-- -->

## ASV abundance

### Explore distribution of ASV read counts


```r
# Extract normalised counts from DESeq object
asv_counts <- counts(dds, normalize = T) %>% as.data.frame()

# Sum ASV counts across samples
total_asv_counts <- rowSums(asv_counts)

# Sort ASVs by abundance
total_asv_counts <- total_asv_counts[order(total_asv_counts, decreasing = T)]

# Caculate cumulative percentage
cumulative <- data.frame(
  cumulative = cumsum(total_asv_counts) / sum(total_asv_counts) * 100,
  no = seq_along(total_asv_counts)
)

# Plot cumulative percentage of ASVs
fun_cum_asv <- ggline(
  data = cumulative, x = "no", y = "cumulative", 
  plot_type = "l", palette = cbPalette,
  title = "Cumulative percentage of fungal ASVs", xlab = "Number of ASVs", 
  ylab = "Cumulative percentage of reads"
)
ggsave(filename = "fun_cum_asv.png", plot = fun_cum_asv, path = "figures/")
fun_cum_asv
```

![](root_endophytes_files/figure-html/FUN top ASVs-1.png)<!-- -->

```r
# Find the number of ASVs that account for 50%, 80%, and 99% of total reads
cat(
  "Number of ASVs that account for 50%, 80%, 90%, and 99% of total reads", "\n\n",
  "50%:", sum(cumulative <= 50), "\n",
  "80%:", sum(cumulative <= 80), "\n",
  "90%:", sum(cumulative <= 90), "\n",
  "99%:", sum(cumulative <= 99), "\n"
)
```

```
# Number of ASVs that account for 50%, 80%, 90%, and 99% of total reads 
# 
#  50%: 57 
#  80%: 140 
#  90%: 238 
#  99%: 741
```

```r
# Find the cumulative percentage accounted for by top x ASVs
cat(
  "Percentage of total reads accounted for by the top 100, 200,and 500 ASVs:", "\n\n",
  "100:", round(cumulative[cumulative$no == 100, "cumulative"], 1) , "\n",
  "200:", round(cumulative[cumulative$no == 200, "cumulative"], 1) , "\n",
  "500:", round(cumulative[cumulative$no == 500, "cumulative"], 1) , "\n"
)
```

```
# Percentage of total reads accounted for by the top 100, 200,and 500 ASVs: 
# 
#  100: 85.9 
#  200: 92.6 
#  500: 98.1
```

```r
# Average ASV counts in order
mean_asv_counts <- rowMeans(asv_counts)
mean_asv_counts <- mean_asv_counts[order(mean_asv_counts, decreasing = T)]

# Plot read count distribution
fun_asv_counts <- ggline(
  data = data.frame(ASV = seq_along(mean_asv_counts), counts = mean_asv_counts),
  x = "ASV", y = "counts", plot_type = "l",
  title = "Fungal ASV read count distribution", xlab = "ASV", ylab = "Mean read count"
)
ggsave(filename = "fun_asv_counts.png", plot = fun_asv_counts, path = "figures/")
fun_asv_counts
```

![](root_endophytes_files/figure-html/FUN top ASVs-2.png)<!-- -->

```r
# Number of ASVs with mean read count > 100, 200, and 500
cat(
  "Number of ASVs with mean read count > 100, 200, and 500", "\n\n",
  "100:", sum(rowMeans(asv_counts) > 100), "\n",
  "200:", sum(rowMeans(asv_counts) > 200), "\n",
  "500:", sum(rowMeans(asv_counts) > 500), "\n"
)
```

```
# Number of ASVs with mean read count > 100, 200, and 500 
# 
#  100: 147 
#  200: 86 
#  500: 49
```

### Filter top ASVs with 100 % of reads


```r
# Filter the top x abundant ASVs by the sum of their normalised counts
# top_asvs <- asv_counts[order(rowSums(asv_counts), decreasing = T)[1:DIFFOTU], ]

# Filter ASVs with mean read count > 100
# top_asvs <- asv_counts[rowMeans(asv_counts) > 100, ]

# Filter top ASVs with 90% of reads
# top_asvs <- asv_counts[filter_otus(asv_counts, FUNASVFILTER), ]

# No filtering
top_asvs <- asv_counts

# Check that sample names match
identical(names(top_asvs), rownames(colData))
```

```
# [1] TRUE
```

```r
# Extract taxonomic data for top ASVs
top_taxa <- taxData[rownames(top_asvs), ]

# Log transform normalised counts
# top_asvs <- log10(top_asvs + 1)

top_asv_data <- data.frame(t(top_asvs))
top_asv_ids <- rownames(top_asvs)

# Check that sample names match
identical(rownames(top_asv_data), rownames(colData))
```

```
# [1] TRUE
```

```r
# Add sample metadata to top ASV data
top_asv_data <- merge(top_asv_data, colData, by = 0) %>% column_to_rownames("Row.names")
```

### Effect of design factors on abundance of top ASVs

Effect of Site, Scion, and Storage on abundance of top ASVs


```r
# ANOVA of top ASVs
asv_lm_anova <- function(asv, formula, data) {
  f = update(formula, paste0("log(", asv, " + 1) ~ ."))
  a = aov(f, data = data) 
  a = a %>% summary() %>% unclass() %>% data.frame()
  return(a)
}

# Negative binomial regression model
asv_negbin_anova <- function(asv, formula, data) {
  f = update(formula, paste0(asv, " ~ ."))
  m = glm.nb(f, data = data)
  a = anova(m, test = "Chisq") %>% data.frame()
  return(a)
}

formula <- FULL_DESIGN

# Full design model does not converge
# formula <- y ~ site + Scion + Storage + site:Scion + site:Storage

# Extend ANOVA results with ASV metadata
extend_asv_anova <- function(anova_result, asv) {
  anova_result %>% mutate(
    ASV = asv,
    Taxonomy = taxData[asv, "rank"],
    Abundance = round(mean(top_asv_data[[asv]]), 3),
    Factor = gsub(" ", "", rownames(anova_result)),
    var = Sum.Sq / sum(anova_result$Sum.Sq) * 100,
    # sig = case_when(
    #   is.na(Pr..F.) ~ "",
    #   Pr..F. < 0.001 ~ "***",
    #   Pr..F. < 0.01 ~ "**",
    #   Pr..F. < 0.05 ~ "*",
    #   TRUE ~ ""
    # ),
    # Variance = ifelse(
    #   perc_var < 0.01, paste0("<0.01", sig),
    #   paste0(round(perc_var, 2), sig)
    # )
  )
}

# Perform ANOVA on list of top ASVs
top_asvs_anova_results <- lapply(
  top_asv_ids, 
  function(asv) {
    asv_lm_anova(asv, formula, top_asv_data) %>%
    extend_asv_anova(asv)
  }
) %>% bind_rows() %>% data.table()

# Group by factor and adjust p-values
top_asvs_anova_results <- top_asvs_anova_results %>% 
  group_by(Factor) %>% 
  mutate(p.adj = p.adjust(`Pr..F.`, method = "BH")) %>% 
  data.table()

# Order factors by original order
top_asvs_anova_results$Factor <- factor(top_asvs_anova_results$Factor, levels = unique(top_asvs_anova_results$Factor))

# Summary of top ASV ANOVA results
top_asvs_anova_summary <- top_asvs_anova_results %>% 
  select(ASV, Taxonomy, Abundance, Factor, var, p.adj) %>% 
  pivot_wider(names_from = Factor, values_from = c(var, p.adj), names_glue = "{Factor}_{.value}") %>%
  select(
    ASV, Taxonomy, Abundance, Site_var, Site_p.adj, Storage_var, Storage_p.adj, 
    Scion_var, Scion_p.adj, "Site:Storage_var", "Site:Storage_p.adj", "Site:Scion_var", 
    "Site:Scion_p.adj", "Storage:Scion_var", "Storage:Scion_p.adj", "Site:Storage:Scion_var", 
    "Site:Storage:Scion_p.adj"
  ) %>%
  data.table()

# top_asvs_anova_summary %>%
#   kbl() %>%
#   kable_styling("striped", full_width = F)

cat(
  "Number of ASVs with statistically significant (*P* < 0.05) adjusted p-values", "\n\n",
  "Site:", nrow(top_asvs_anova_results[Factor == "Site" & p.adj < 0.05, ]), "\n",
  "Storage:", nrow(top_asvs_anova_results[Factor == "Storage" & p.adj < 0.05, ]), "\n",
  "Scion:", nrow(top_asvs_anova_results[Factor == "Scion" & p.adj < 0.05, ]), "\n",
  "Site:Storage:", nrow(top_asvs_anova_results[Factor == "Site:Storage" & p.adj < 0.05, ]), "\n",
  "Site:Scion:", nrow(top_asvs_anova_results[Factor == "Site:Scion" & p.adj < 0.05, ]), "\n",
  "Storage:Scion:", nrow(top_asvs_anova_results[Factor == "Storage:Scion" & p.adj < 0.05, ]), "\n",
  "Site:Storage:Scion:", nrow(top_asvs_anova_results[Factor == "Site:Storage:Scion" & p.adj < 0.05, ]), "\n\n",
  "Total ASVs:", length(unique(top_asvs_anova_results$ASV)), "\n\n"
)

# Filter by significant effect of scion and its interactions
scion_asvs <- top_asvs_anova_results[grepl("Scion", Factor) & p.adj < 0.05, ]

# scion_asvs %>%
#   kbl() %>%
#   kable_styling("striped", full_width = F)

cat(
  length(scion_asvs$ASV), 
  "ASVs with significant (*P* < 0.05) adjusted p-values for the effect of Scion and its interactions.", "\n\n"
)

# Summary of ASVs with significant Scion effect
fun_significant_scion <- top_asvs_anova_summary[ASV %in% scion_asvs$ASV, ] %>%
  arrange(desc(Abundance))

# kbl(fun_significant_scion) %>%
#   kable_styling("striped")

kbl(fun_significant_scion) %>%
  kable_styling("striped") %>%
  save_kable("tables/FUN_scion_asvs.html")

top_asvs_anova_summary %>%
  arrange(desc(Abundance)) %>%
  kbl() %>%
  kable_styling("striped") %>%
  save_kable("tables/FUN_all_asvs.html")
```


```r
# Export significant ASVs as fasta

# Read FUN ASVs
FUN_asvs <- read.fasta("data/FUN.zotus.fa")
# Replace 'OTU' with 'ASV' in sequence names
names(FUN_asvs) <- gsub("OTU", "ASV", names(FUN_asvs))

# Write significant ASVs to fasta
write.fasta(
  sequences = FUN_asvs[scion_asvs$ASV],
  names = paste(scion_asvs$ASV, taxData[scion_asvs$ASV, "rank"]),
  file = "fasta/FUN_scion_asvs.fasta"
)
```

## Canker counts

Testing the effects of of total abundance, ASV abundance, α-diversity, and β-diversity on canker counts.

This uses a nested negative binomial regression model.

The base model for canker counts uses the formula: Cankers ~ Site * Storage * Scion.


```r
# Filter out samples with missing canker count
canker_abundance_data <- colData[complete.cases(colData$Cankers), ]

# Base model
canker_design <- "Cankers ~ Site * Storage * Scion"
base_model <- glm.nb(canker_design, data = canker_abundance_data)

# Abundance model
abundance_design <- paste(canker_design, "+ log(copy_number)")
abundance_model <- glm.nb(abundance_design, data = canker_abundance_data)

# ANOVA of abundance with canker count
anova(base_model, abundance_model) %>%
  kbl() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Model </th>
   <th style="text-align:right;"> theta </th>
   <th style="text-align:right;"> Resid. df </th>
   <th style="text-align:right;">    2 x log-lik. </th>
   <th style="text-align:left;"> Test </th>
   <th style="text-align:right;">    df </th>
   <th style="text-align:right;"> LR stat. </th>
   <th style="text-align:right;"> Pr(Chi) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site * Storage * Scion </td>
   <td style="text-align:right;"> 2.828262 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> -546.8670 </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site * Storage * Scion + log(copy_number) </td>
   <td style="text-align:right;"> 2.831854 </td>
   <td style="text-align:right;"> 31 </td>
   <td style="text-align:right;"> -546.8376 </td>
   <td style="text-align:left;"> 1 vs 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0293702 </td>
   <td style="text-align:right;"> 0.863927 </td>
  </tr>
</tbody>
</table>

### Effect of ASV abundance on canker count


```r
# Filter out samples with missing canker count
canker_top_asv_data <- top_asv_data[complete.cases(top_asv_data$Cankers), ]
# all.equal(top_asv_data[c("Site", "Storage", "Scion", "Cankers")],cankers)

# Base model design
canker_design <- "Cankers ~ Site * Storage * Scion"

# Base model with ASV abundance data
base_model <- glm.nb(canker_design, data = canker_top_asv_data)

# Fits glm.nb model and returns a list of the model and warnings
glm.nb_with_warnings <- function(f, data) {
  m = tryCatch(
    {
      list(
        fit = glm.nb(f, data = data),
        warning = NA
      )
    },
    warning = function(w) {
      list(
        fit = glm.nb(f, data = data),
        warning = paste(conditionMessage(w), collapse = ", ")
      )
    }
  )
  return(m)
}

# ANOVA of top ASVs with canker count
asv_canker_anova <- function(asv, design, base_model, data) {
  tryCatch({
    log_asv = paste0("log(", asv, " + 1)")
    f = paste(design, "+", log_asv)#, "+", log_asv, ":Site")
    m = glm.nb_with_warnings(f, data)
    a = anova(base_model, m$fit) %>% data.frame()
    b = suppressWarnings(anova(m$fit)) %>% data.frame()
    total_deviance = sum(b$Deviance, na.rm = T) + tail(b$Resid..Dev, 1)
    d = data.frame(
      ASV = asv,
      Taxonomy = taxData[asv, "rank"],
      Abundance = round(mean(data[[asv]]), 3),
      coef = round(m$fit$coefficients[log_asv], 3),
      var = round(b[log_asv, 'Deviance'] / total_deviance * 100, 3),
      p = a[2, 'Pr.Chi.'],
      warning = m$warning
    )
    return(d)
  }, warning = function(w) {
    return(data.frame(
      ASV = asv,
      Taxonomy = taxData[asv, "rank"],
      Abundance = round(mean(data[[asv]])),
      coef = NA,
      var = NA,
      p = NA,
      warning = paste("Warning: ", conditionMessage(w))
    ))
  }, error = function(e) {
    return(data.frame(
      ASV = asv,
      Taxonomy = taxData[asv, "rank"],
      Abundance = round(mean(data[[asv]])),
      coef = NA,
      var = NA,
      p = NA,
      warning = paste("Error: ", conditionMessage(e))
    ))
  })
}

# Effect of ASV abundance on canker count for top ASVs
asv_canker_results <- lapply(
  top_asv_ids, 
  function(x) asv_canker_anova(x, canker_design, base_model, canker_top_asv_data)
) %>% bind_rows() %>% data.table()

# Adjust p-values for multiple testing
asv_canker_results$p_adjusted <- p.adjust(asv_canker_results$p, method = "BH")

# Summary of ASVs with statistically significant (*P* < 0.05) adjusted p-values
cat(
  nrow(asv_canker_results[p_adjusted < 0.05, ]), "of", nrow(asv_canker_results),
  "ASVs have statistically significant (*P* < 0.05) adjusted p-values\n\n"
)

fun_sig_canker_asvs <- asv_canker_results[p_adjusted < 0.05, ] %>%
  arrange(desc(Abundance))

# kbl(fun_sig_canker_asvs) %>%
#   kable_styling("striped")

kbl(fun_sig_canker_asvs) %>%
  kable_styling("striped") %>%
  save_kable("tables/FUN_canker_asvs.html")

cat(
  "\nSignificant canker bacterial ASVs\n",
  "\nTotal ASVs: ", nrow(asv_canker_results),
  "\nWarnings: ", nrow(asv_canker_results[!is.na(warning)]),
  "\nTotal ASVs without warnings: ", nrow(asv_canker_results[is.na(warning)]),
  "\nTotal significant ASVs: ", nrow(asv_canker_results[p_adjusted < 0.05]),
  "\nASVs with negative effect: ", nrow(asv_canker_results[p_adjusted < 0.05 & coef < 0]),
  "\nASVs with positive effect: ", nrow(asv_canker_results[p_adjusted < 0.05 & coef > 0]), "\n"
)
```

#### Effect of ASV abundance on canker count per site

Filter top ASVs with 100 % of reads per site
and test the effect of ASV abundance on canker count per site.


```r
# For each site, select top ASVs with 90% of reads
top_asvs_per_site <- lapply(
  unique(colData$Site),
  function(site) {
    samples <- filter(colData, Site == site)
    top_asv_data <- select(asv_counts, rownames(samples))
    # top_asvs <- filter(top_asv_data, rowMeans(top_asv_data) > 100)
    # top_asvs <- top_asv_data[filter_otus(top_asv_data, FUNASVFILTER), ]
    top_asvs <- top_asv_data[rowSums(top_asv_data) > 0, ]
    top_asv_ids <- rownames(top_asvs)
    top_asvs <- data.frame(t(top_asvs)) %>% merge(samples, by = 0) %>% column_to_rownames("Row.names")
    top_asvs <- top_asvs[complete.cases(top_asvs$Cankers), ]
    return(list(asvs = top_asv_ids, data = top_asvs))
  }
)

# Summary of top ASVs per site
data.table(
  Site = 1:3,
  topASVs = sapply(top_asvs_per_site, function(x) length(x$asvs)),
  totalASVs = sapply(top_asvs_per_site, function(x) sum(rowSums(asv_counts[, rownames(x$data)]) > 0, na.rm = TRUE)),
  topreads = sapply(top_asvs_per_site, function(x) sum(select(x$data, starts_with("ASV")))),
  totalreads = sapply(top_asvs_per_site, function(x) sum(asv_counts[, rownames(x$data)]))
) %>% mutate(
  "ASV%" = round(topASVs / totalASVs * 100),
  "Read%" = round(topreads / totalreads * 100)
) %>% kbl() %>% kable_styling("striped")  

canker_site_design <- "Cankers ~ Storage * Scion"

# ANOVA of ASV abundance with canker count per ASV
asv_canker_site_anova <- function(asvs, data) {
  base_model <- glm.nb(canker_site_design, data = data)
  results <- lapply(
    asvs, 
    function(asv) asv_canker_anova(asv, canker_site_design, base_model, data)
  ) %>% bind_rows() %>% data.table()
  results$p_adjusted <- p.adjust(results$p, method = "BH")
  return(results)
}

# Run ANOVA per site
asv_canker_site_results <- lapply(
  top_asvs_per_site,
  function(x) asv_canker_site_anova(x$asvs, x$data)
)

# Add site to each result as new column and merge into single data.table
asv_canker_site_results <- lapply(
  1:3, 
  function(site) {
    result <- asv_canker_site_results[[site]]
    result$Site <- site
    result
  }
) %>% bind_rows()

# Significant ASVs
significant_asvs <- asv_canker_site_results[p_adjusted < 0.05 & is.na(warning), ]

fun_canker_site_asvs <- significant_asvs[, c("Site", "ASV", "Taxonomy", "Abundance", "coef", "var", "p_adjusted", "warning")] %>%
  arrange(Site, desc(Abundance))

# kbl(fun_canker_site_asvs) %>%
#   kable_styling("striped")

kbl(fun_canker_site_asvs) %>%
  kable_styling("striped") %>%
  save_kable("tables/FUN_canker_site_asvs.html")

for(site in 1:3){
  cat(
    "\nSignificant canker fungal ASVs for site ", site, "\n",
    "\nTotal ASVs: ", nrow(asv_canker_site_results[Site == site]),
    "\nWarnings: ", nrow(asv_canker_site_results[Site == site & !is.na(warning)]),
    "\nTotal ASVs without warnings: ", nrow(asv_canker_site_results[Site == site & is.na(warning)]),
    "\nTotal significant ASVs: ", nrow(asv_canker_site_results[Site == site & p_adjusted < 0.05]),
    "\nASVs with negative effect: ", nrow(asv_canker_site_results[Site == site & p_adjusted < 0.05 & coef < 0]),
    "\nASVs with positive effect: ", nrow(asv_canker_site_results[Site == site & p_adjusted < 0.05 & coef > 0]), "\n"
  )
}
```


```r
# Export significant ASVs as FASTA
write.fasta(
  sequences = FUN_asvs[as.character(significant_asvs$ASV)],
  names = paste(significant_asvs$ASV, taxData[as.character(significant_asvs$ASV), "rank"]),
  file = "fasta/FUN_canker_asvs.fasta"
)
```

```
# Error in eval(expr, envir, enclos): object 'FUN_asvs' not found
```

##### Plot of ASV abundance against canker count


```r
# List of significant ASVs
significant_asv_list <- significant_asvs$ASV %>% unlist()

significant_asv_data <- asv_counts[significant_asv_list, ] %>% 
  t() %>% 
  data.frame() %>% 
  merge(colData, by = 0) %>% 
  column_to_rownames("Row.names") %>%
  select(c(significant_asv_list, "Site", "Storage", "Scion", "Cankers"))

# Melt data for ggplot
significant_asv_long_data <- significant_asv_data %>% reshape2::melt(
  id.vars = c("Site", "Storage", "Scion", "Cankers"), variable.name = "ASV", value.name = "Abundance"
)

# Log trasnform abundance
significant_asv_long_data$log10_abundance <- log10(significant_asv_long_data$Abundance + 1)

fun_asv_canker_plot <- ggscatter(
  data = significant_asv_long_data, x = "log10_abundance", y = "Cankers", 
  color = "Storage", facet.by = c("ASV", "Site"),
  xlab = "ASV abundance (log10)", ylab = "Canker count",
  palette = cbPalette, legend = "bottom"
)

ggsave(
  filename = "fun_asv_canker_plot.png", plot = fun_asv_canker_plot, path = "figures/",
  height = 40, width = 20, units = "cm"
)

fun_asv_canker_plot
```

#### Effect of aggregated genera on canker counts


```r
# Add genus from taxData to countData
fun_genus_data <- counts(dds, normalize = T) %>% as.data.frame() %>% mutate(
  genus = taxData[rownames(countData), "genus"]
)

# Group by genus
fun_genus_data <- fun_genus_data %>% group_by(genus) %>% summarise_all(sum) %>% as.data.frame()

# Set rownames as genus
rownames(fun_genus_data) <- fun_genus_data$genus
fun_genus_data <- dplyr::select(fun_genus_data, -genus)

# Filter genera with mean abundance < 100
fun_genus_data <- fun_genus_data[rowMeans(fun_genus_data) > 10, ]

# Rank not genus
not_genus <- rownames(fun_genus_data)[grep("\\([a-z]\\)", rownames(fun_genus_data))]
# Remove rows with genus in not_genus
fun_genus_data <- fun_genus_data[!rownames(fun_genus_data) %in% not_genus, ]
cat(
  length(not_genus), "non-genus ranks removed:\n\n",
  not_genus, "\n"
)

# Final genus list
fun_genera <- rownames(fun_genus_data)

# Transpose and add metadata from colData
fun_genus_data <- t(fun_genus_data) %>% as.data.frame() %>% mutate(
  Site = colData[rownames(.), "Site"],
  Storage = colData[rownames(.), "Storage"],
  Scion = colData[rownames(.), "Scion"],
  Cankers = colData[rownames(.), "Cankers"]
)

# Filter out samples with missing canker count
fun_genus_data <- fun_genus_data[complete.cases(fun_genus_data$Cankers), ]
```


```r
# Base model design
canker_design = "Cankers ~ Site * Storage * Scion"

# Base model with genus abundance data
base_model <- glm.nb(canker_design, data = fun_genus_data)

# ANOVA of genus abundance with canker count
genus_canker_anova <- function(genus, design, base_model, data) {
  log_genus = paste0("log(", genus, " + 1)")
  f = paste(design, "+", log_genus)#, "+", log_genus, ":Site")
  m = glm.nb(f, data = data)
  a = anova(base_model, m) %>% data.frame()
  b = suppressWarnings(anova(m)) %>% data.frame()
  total_deviance = sum(b$Deviance, na.rm = T) + tail(b$Resid..Dev, 1)
  d = data.table(
    Genus = genus,
    Abundance = mean(data[[genus]]),
    coef = m$coefficients[log_genus],
    var = b[log_genus, 'Deviance'] / total_deviance * 100,
    p = a[2, 'Pr.Chi.']
  )
  return(d)
}

# genus_canker_anova(fun_genera[1], canker_design, base_model, fun_genus_data)

# Effect of genera abundance on canker counts
genus_canker_results <- lapply(
  fun_genera, 
  function(x) genus_canker_anova(x, canker_design, base_model, fun_genus_data)
) %>% bind_rows() %>% data.table()

# Adjust p-values for multiple testing
genus_canker_results$p_adjusted <- p.adjust(genus_canker_results$p, method = "BH")

# Summary of genera with statistically significant (*P* < 0.05) adjusted p-values
cat(
  nrow(genus_canker_results[p_adjusted < 0.05, ]), "of", nrow(genus_canker_results),
  "genera have statistically significant (*P* < 0.05) adjusted p-values\n\n"
)
if(nrow(genus_canker_results[p_adjusted < 0.05, ]) > 0) {
  genus_canker_results[p_adjusted < 0.05, ] %>%
    arrange(desc(Abundance)) %>%
    kbl() %>%
    kable_styling("striped")
}

# for (genus in sig_genus) {
#   log_genus = paste0("log(", genus, " + 1)")
#   f = paste(canker_design, "+", log_genus)
#   m = glm.nb(f, data = fun_genus_data)
#   print(anova(base_model, m))
# }
```

##### Plot of genus abundance against canker count


```r
significant_genera <- genus_canker_results[p_adjusted < 0.05]$Genus %>% unlist()

significant_genera_data <- fun_genus_data[, c(significant_genera, FACTORS, "Cankers")]

# Melt data for ggplot
significant_genera_data <- significant_genera_data %>% reshape2::melt(
  id.vars = c(FACTORS, "Cankers"), variable.name = "Genus", value.name = "Abundance"
)

# Log transform abundance
significant_genera_data$log10_abundance <- log(significant_genera_data$Abundance + 1)
# significant_genera_data$log10_cankers <- log(significant_genera_data$Cankers + 1)

# Plot of genus abundance against canker count
fun_genus_canker_plot <- ggscatter(
  data = significant_genera_data, x = "log10_abundance", y = "Cankers", 
  color = "Site", shape = "Storage",
  xlab = "Genus abundance (log10)", ylab = "Canker count",
  free_x = T, free_y = T, palette = cbPalette
) + facet_wrap(~ Genus, scales = "free")

ggsave(
  filename = "fun_genus_canker_plot.png", plot = fun_genus_canker_plot, path = "figures/",
  height = 10, width = 20, units = "cm"
)
```

### Effect of α-diversity on canker count


```r
# ANOVA of α-diversity with canker count

# Base model with α-diversity data
base_model <- glm.nb(canker_design, data = all_alpha_ord)

measures <- c("S.chao1", "shannon", "simpson")

# ANOVA of α-diversity with canker count
alpha_canker_anova <- function(measure, data) {
  f = paste(canker_design, "+", measure)
  m = glm.nb(f, data = data)
  a = anova(base_model, m) %>% data.frame()
  b = anova(m) %>% data.frame()
  total_deviance = sum(b$Deviance, na.rm = T) + tail(b$Resid..Dev, 1)
  d = data.frame(
    measure = measure,
    df = a[2, 'X...df'],
    coef = m$coefficients[measure],
    lr = a[2, 'LR.stat.'],
    p = a[2, 'Pr.Chi.'],
    var = b[measure, 'Deviance'] / total_deviance * 100
  )
  return(d)
}

# alpha_canker_anova("shannon", all_alpha_ord)

# Effect of α-diversity on canker count for each measure
alpha_canker_results <- data.table(bind_rows(lapply(measures, function(x) alpha_canker_anova(x, all_alpha_ord))))

alpha_canker_results %>%
  kbl() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> measure </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> coef </th>
   <th style="text-align:right;"> lr </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:right;"> var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> S.chao1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0033377 </td>
   <td style="text-align:right;"> 2.9113117 </td>
   <td style="text-align:right;"> 0.0879603 </td>
   <td style="text-align:right;"> 0.0384239 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shannon </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0830027 </td>
   <td style="text-align:right;"> 0.0720829 </td>
   <td style="text-align:right;"> 0.7883279 </td>
   <td style="text-align:right;"> 0.2109241 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> simpson </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -4.0314477 </td>
   <td style="text-align:right;"> 4.1404456 </td>
   <td style="text-align:right;"> 0.0418701 </td>
   <td style="text-align:right;"> 1.8197613 </td>
  </tr>
</tbody>
</table>

```r
# ANOVA results
for (measure in measures) {
  f = paste(canker_design, "+", measure)
  m = glm.nb(f, data = all_alpha_ord)
  print(anova(base_model, m))
}
```

```
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: Cankers
#                              Model    theta Resid. df    2 x log-lik.   Test
# 1           Site * Storage * Scion 2.828262        32       -546.8670       
# 2 Site * Storage * Scion + S.chao1 2.934597        31       -543.9557 1 vs 2
#      df LR stat.    Pr(Chi)
# 1                          
# 2     1 2.911312 0.08796031
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: Cankers
#                              Model    theta Resid. df    2 x log-lik.   Test
# 1           Site * Storage * Scion 2.828262        32       -546.8670       
# 2 Site * Storage * Scion + shannon 2.836456        31       -546.7949 1 vs 2
#      df   LR stat.   Pr(Chi)
# 1                           
# 2     1 0.07208289 0.7883279
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: Cankers
#                              Model    theta Resid. df    2 x log-lik.   Test
# 1           Site * Storage * Scion 2.828262        32       -546.8670       
# 2 Site * Storage * Scion + simpson 3.056752        31       -542.7265 1 vs 2
#      df LR stat.    Pr(Chi)
# 1                          
# 2     1 4.140446 0.04187015
```

### Effect of β-diversity on canker count


```r
no_pcs <- 10

# Merge PC scores with canker data
pc_scores <- merge(colData, data.frame(mypca$x[, 1:no_pcs]), by = "row.names") %>% 
  column_to_rownames("Row.names")

pcs <- tail(colnames(pc_scores), no_pcs)

# Base model with β-diversity data
base_model <- glm.nb(canker_design, data = pc_scores)

# ANOVA of β-diversity with canker count
beta_canker_anova <- function(pc, data) {
  f = paste0(canker_design, "+", pc)
  m = glm.nb(f, data = data)
  a = anova(base_model, m) %>% data.frame()
  b = anova(m) %>% data.frame()
  total_deviance = sum(b$Deviance, na.rm = T) + tail(b$Resid..Dev, 1)
  d = data.frame(
    PC = pc,
    df = a[2, 'X...df'],
    coef = m$coefficients[pc],
    lr = a[2, 'LR.stat.'],
    p = a[2, 'Pr.Chi.'],
    var = b[pc, 'Deviance'] / total_deviance * 100
  )
  return(d)
}

# Effect of β-diversity on canker count for each PC
beta_canker_results <- data.table(bind_rows(lapply(pcs, function(x) beta_canker_anova(x, pc_scores))))

beta_canker_results %>%
  kbl() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> coef </th>
   <th style="text-align:right;"> lr </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:right;"> var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0128138 </td>
   <td style="text-align:right;"> 0.9851543 </td>
   <td style="text-align:right;"> 0.3209296 </td>
   <td style="text-align:right;"> 4.6547052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0025164 </td>
   <td style="text-align:right;"> 0.1365507 </td>
   <td style="text-align:right;"> 0.7117345 </td>
   <td style="text-align:right;"> 0.0513068 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0119703 </td>
   <td style="text-align:right;"> 2.5742584 </td>
   <td style="text-align:right;"> 0.1086150 </td>
   <td style="text-align:right;"> 0.1146305 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0263314 </td>
   <td style="text-align:right;"> 3.3395217 </td>
   <td style="text-align:right;"> 0.0676343 </td>
   <td style="text-align:right;"> 1.1475272 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0047435 </td>
   <td style="text-align:right;"> 0.0719010 </td>
   <td style="text-align:right;"> 0.7885887 </td>
   <td style="text-align:right;"> 0.5930129 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0264086 </td>
   <td style="text-align:right;"> 4.3863504 </td>
   <td style="text-align:right;"> 0.0362278 </td>
   <td style="text-align:right;"> 1.4714388 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0058366 </td>
   <td style="text-align:right;"> 0.1431624 </td>
   <td style="text-align:right;"> 0.7051572 </td>
   <td style="text-align:right;"> 0.0274986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0093899 </td>
   <td style="text-align:right;"> 0.3581093 </td>
   <td style="text-align:right;"> 0.5495581 </td>
   <td style="text-align:right;"> 0.8596759 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0185098 </td>
   <td style="text-align:right;"> 1.0450939 </td>
   <td style="text-align:right;"> 0.3066397 </td>
   <td style="text-align:right;"> 0.1436870 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0065698 </td>
   <td style="text-align:right;"> 0.1379981 </td>
   <td style="text-align:right;"> 0.7102794 </td>
   <td style="text-align:right;"> 0.7028716 </td>
  </tr>
</tbody>
</table>


```r
# Save environment
save.image("FUN.RData")
```

<!-- #=============================================================================== -->
# **Bacteria**
<!-- #=============================================================================== -->


```r
# Unpack bacteria data
invisible(mapply(assign, names(ubiome_BAC), ubiome_BAC, MoreArgs = list(envir = globalenv())))
```

## ASV and sample summary

### Read and sample summary


```r
cat(
  "Raw reads", "\n\n",
  "Total raw reads:\t\t", sum(countData), "\n",
  "Mean raw reads per sample:\t", mean(colSums(countData)), "\n",
  "Median raw reads per sample:\t", median(colSums(countData)), "\n",
  "Max raw reads per sample:\t", max(colSums(countData)), "\n",
  "Min raw reads per sample:\t", min(colSums(countData)), "\n\n"
)
```

```
# Raw reads 
# 
#  Total raw reads:		 3365816 
#  Mean raw reads per sample:	 41046.54 
#  Median raw reads per sample:	 40406.5 
#  Max raw reads per sample:	 89023 
#  Min raw reads per sample:	 15049
```

```r
#colSums(countData)

nct <- counts(dds, normalize = T)
cat("Normalised reads", "\n\n",
  "Total normalised reads:\t\t", sum(nct), "\n",
  "Mean normalised reads per sample:\t", mean(colSums(nct)), "\n",
  "Median normalised reads per sample:\t", median(colSums(nct)), "\n",
  "Min normalised reads per sample:\t", min(colSums(nct)), "\n",
  "Max normalised reads per sample:\t", max(colSums(nct)), "\n\n"
)
```

```
# Normalised reads 
# 
#  Total normalised reads:		 4585124 
#  Mean normalised reads per sample:	 55916.14 
#  Median normalised reads per sample:	 53407.32 
#  Min normalised reads per sample:	 9940.045 
#  Max normalised reads per sample:	 139399.1
```

```r
#round(colSums(counts(dds,normalize = T)),0)
```

### ASV summary 


```r
cat(
  "Total ASVs:\t\t", nrow(taxData),"\n\n",
  "Raw reads per ASV summary", "\n\n",
  "Mean raw reads per ASV:\t", mean(rowSums(countData)),"\n",
  "Median raw per ASV:\t\t", median(rowSums(countData)),"\n",
  "ASV raw Min reads:\t\t", min(rowSums(countData)),"\n",
  "ASV raw Max reads:\t\t", max(rowSums(countData)),"\n\n"
)
```

```
# Total ASVs:		 5883 
# 
#  Raw reads per ASV summary 
# 
#  Mean raw reads per ASV:	 572.1258 
#  Median raw per ASV:		 120 
#  ASV raw Min reads:		 34 
#  ASV raw Max reads:		 106398
```

```r
cat(
  "Normalised reads per ASV summary","\n\n",
  "Mean normalised reads per ASV:\t\t", mean(rowSums(nct)),"\n",
  "Median normalised reads per ASV:\t", median(rowSums(nct)),"\n",
  "ASV normalised Min reads:\t\t", min(rowSums(nct)),"\n",
  "ASV normalised Max reads:\t\t", max(rowSums(nct)),"\n\n"
)
```

```
# Normalised reads per ASV summary 
# 
#  Mean normalised reads per ASV:		 779.3853 
#  Median normalised reads per ASV:	 156.6744 
#  ASV normalised Min reads:		 20.59746 
#  ASV normalised Max reads:		 139400.3
```

```r
y <- rowSums(nct)
y <- y[order(y, decreasing = T)]
# proportion
xy <- y/sum(y)

cat("Top ", TOPOTU, "ASVs:\n")
```

```
# Top  10 ASVs:
```

```r
data.frame(
  counts = y[1:TOPOTU], 
  proportion = xy[1:TOPOTU], 
  rank = taxData[names(y)[1:TOPOTU],]$rank
) %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> counts </th>
   <th style="text-align:right;"> proportion </th>
   <th style="text-align:left;"> rank </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> ASV1 </td>
   <td style="text-align:right;"> 139400.30 </td>
   <td style="text-align:right;"> 0.0304027 </td>
   <td style="text-align:left;"> Streptomyces(g) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV2 </td>
   <td style="text-align:right;"> 124366.87 </td>
   <td style="text-align:right;"> 0.0271240 </td>
   <td style="text-align:left;"> Kineosporiales(o) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV3 </td>
   <td style="text-align:right;"> 113814.51 </td>
   <td style="text-align:right;"> 0.0248226 </td>
   <td style="text-align:left;"> Kineosporiaceae(f) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV4 </td>
   <td style="text-align:right;"> 99723.88 </td>
   <td style="text-align:right;"> 0.0217494 </td>
   <td style="text-align:left;"> Streptomyces(g) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV5 </td>
   <td style="text-align:right;"> 96979.25 </td>
   <td style="text-align:right;"> 0.0211508 </td>
   <td style="text-align:left;"> Streptomyces(g) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV6 </td>
   <td style="text-align:right;"> 79618.94 </td>
   <td style="text-align:right;"> 0.0173646 </td>
   <td style="text-align:left;"> Streptomyces(g) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV7 </td>
   <td style="text-align:right;"> 51308.93 </td>
   <td style="text-align:right;"> 0.0111903 </td>
   <td style="text-align:left;"> Bradyrhizobium(g) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV8 </td>
   <td style="text-align:right;"> 49535.82 </td>
   <td style="text-align:right;"> 0.0108036 </td>
   <td style="text-align:left;"> Actinoplanes(g) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV20 </td>
   <td style="text-align:right;"> 33817.55 </td>
   <td style="text-align:right;"> 0.0073755 </td>
   <td style="text-align:left;"> Actinobacteria(c) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ASV9 </td>
   <td style="text-align:right;"> 33580.42 </td>
   <td style="text-align:right;"> 0.0073238 </td>
   <td style="text-align:left;"> Nonomuraea(g) </td>
  </tr>
</tbody>
</table>

## Taxonomy Summary

### Taxonomy identifiable

Proportion of ASVs which can be assigned (with the given confidence) at each taxonomic rank


```r
# Proportion of ASVs which can be assigned (with the given confidence) at each taxonomic rank

tx <- copy(taxData)
setDT(tx)
cols <- names(tx)[9:15]

tx[, (cols) := lapply(.SD, as.factor), .SDcols = cols]

data.table(
  rank = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
  "0.8" = round(unlist(lapply(cols, function(col) sum(as.number(tx[[col]]) >= 0.8) / nrow(tx))), 2),
  "0.65" = round(unlist(lapply(cols, function(col) sum(as.number(tx[[col]]) >= 0.65) / nrow(tx))), 2),
  "0.5" = round(unlist(lapply(cols, function(col) sum(as.number(tx[[col]]) >= 0.5) / nrow(tx))), 2)
) %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> rank </th>
   <th style="text-align:right;"> 0.8 </th>
   <th style="text-align:right;"> 0.65 </th>
   <th style="text-align:right;"> 0.5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> kingdom </td>
   <td style="text-align:right;"> 1.00 </td>
   <td style="text-align:right;"> 1.00 </td>
   <td style="text-align:right;"> 1.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> phylum </td>
   <td style="text-align:right;"> 0.94 </td>
   <td style="text-align:right;"> 0.97 </td>
   <td style="text-align:right;"> 0.99 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> class </td>
   <td style="text-align:right;"> 0.84 </td>
   <td style="text-align:right;"> 0.90 </td>
   <td style="text-align:right;"> 0.93 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> order </td>
   <td style="text-align:right;"> 0.65 </td>
   <td style="text-align:right;"> 0.72 </td>
   <td style="text-align:right;"> 0.79 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> family </td>
   <td style="text-align:right;"> 0.51 </td>
   <td style="text-align:right;"> 0.57 </td>
   <td style="text-align:right;"> 0.64 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> genus </td>
   <td style="text-align:right;"> 0.35 </td>
   <td style="text-align:right;"> 0.43 </td>
   <td style="text-align:right;"> 0.51 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> species </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
  </tr>
</tbody>
</table>

% of reads which can be assigned to each taxonomic ranks


```r
tx <-taxData[rownames(dds),]
nc <- counts(dds, normalize = T)
ac <- sum(nc)

data.table(
  rank = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
  "0.8" = round(unlist(lapply(cols, function(col)(sum(nc[which(as.numeric(tx[[col]]) >= 0.8),]) / ac * 100))), 2),
  "0.65" = round(unlist(lapply(cols, function(col)(sum(nc[which(as.numeric(tx[[col]]) >= 0.65),]) / ac * 100))), 2),
  "0.5" = round(unlist(lapply(cols, function(col)(sum(nc[which(as.numeric(tx[[col]]) >= 0.5),]) / ac * 100))), 2)
) %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> rank </th>
   <th style="text-align:right;"> 0.8 </th>
   <th style="text-align:right;"> 0.65 </th>
   <th style="text-align:right;"> 0.5 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> kingdom </td>
   <td style="text-align:right;"> 99.95 </td>
   <td style="text-align:right;"> 100.00 </td>
   <td style="text-align:right;"> 100.00 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> phylum </td>
   <td style="text-align:right;"> 97.78 </td>
   <td style="text-align:right;"> 99.43 </td>
   <td style="text-align:right;"> 99.75 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> class </td>
   <td style="text-align:right;"> 92.63 </td>
   <td style="text-align:right;"> 95.82 </td>
   <td style="text-align:right;"> 98.05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> order </td>
   <td style="text-align:right;"> 71.84 </td>
   <td style="text-align:right;"> 81.25 </td>
   <td style="text-align:right;"> 87.03 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> family </td>
   <td style="text-align:right;"> 58.63 </td>
   <td style="text-align:right;"> 67.80 </td>
   <td style="text-align:right;"> 73.01 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> genus </td>
   <td style="text-align:right;"> 44.35 </td>
   <td style="text-align:right;"> 51.20 </td>
   <td style="text-align:right;"> 58.33 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> species </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
   <td style="text-align:right;"> 0.00 </td>
  </tr>
</tbody>
</table>

### Taxonomy plots

Plots of proportion of normalised reads assigned to members of phylum and class.


```r
dat <- list(as.data.frame(counts(dds, normalize = T)), taxData, as.data.frame(colData(dds)))

design <- c("Site", "Storage")

# md1 <- getSummedTaxa(dat, conf = TAXCONF, design = design, cutoff = 0.1)
md1 <- getSummedTaxa(dat, conf = TAXCONF, design = design, taxon = "phylum", cutoff = 0.1)

md1[, Site := factor(Site, levels = c(1, 2, 3))]
md1[, Storage := factor(Storage, levels = c("no", "yes"))]
md1[, taxon := factor(taxon, levels = unique(taxon[order(value, decreasing = T)]))]

removals <- md1[, .(value = mean(value)), by = "taxon"][value < 0.5, taxon]
md1 <- md1[!taxon %in% removals, ]

bac_phylum_plot <- plotfun1(md1, x = "taxon", fill = "Site") +
  facet_wrap(~ Storage)

ggsave("figures/bac_phylum.png", bac_phylum_plot, width = 25, height = 15, units = "cm")

bac_phylum_plot
```

![](root_endophytes_files/figure-html/BAC taxonomy plots-1.png)<!-- -->

```r
md2 <- getSummedTaxa(dat, conf = TAXCONF, design = design, taxon = "class", cutoff = 0.1, topn = 9)

md2[, Site := factor(Site, levels = c(1, 2, 3))]
md2[, Storage := factor(Storage, levels = c("no", "yes"))]
md2[, taxon := factor(taxon, levels = unique(taxon[order(value, decreasing = T)]))]

removals <- md2[, .(value = mean(value)), by = "taxon"][value < 0.5, taxon]
md2 <- md2[!taxon %in% removals, ]

bac_class_plot <- plotfun1(md2, x = "taxon", fill = "Site") +
  facet_wrap(~ Storage)

ggsave("figures/bac_class.png", bac_class_plot, width = 25, height = 15, units = "cm")

bac_class_plot
```

![](root_endophytes_files/figure-html/BAC taxonomy plots-2.png)<!-- -->

## Community size


```r
abundance_plot <- ggplot(
  data = as.data.frame(colData(dds)), 
  aes(x = Site, y = log_copy_number, colour = Scion, shape = Storage)
) + geom_jitter() + 
  scale_colour_manual(values = cbPalette)

abundance_plot <- ggboxplot(
  data = as.data.frame(colData(dds)), x = "Site", y = "log_copy_number", 
  color = "Scion", add = "jitter", legend = "top", 
  title = "Bacterial abundance", xlab = "Site", ylab = "log10 copy number"
)

ggsave(
  filename = "bac_abundance.png", plot = abundance_plot, path = "figures/", 
  height = 20, width = 20, units = "cm"
)

abundance_plot
```

![](root_endophytes_files/figure-html/BAC size-1.png)<!-- -->

```r
# Formula for ANOVA
formula <- update(FULL_DESIGN, log_copy_number ~ .)

abundance_anova <- aov(formula, data = as.data.frame(colData(dds)))

# Normality check
par(mfrow = c(2, 2))
plot(abundance_anova)
```

![](root_endophytes_files/figure-html/BAC size-2.png)<!-- -->

```r
png("figures/bac_abundance_norm.png", width = 800, height = 600)
par(mfrow = c(2, 2))
plot(abundance_anova)
dev.off()
```

```
# png 
#   2
```

```r
# Results
summary(abundance_anova)
```

```
#                    Df Sum Sq Mean Sq F value   Pr(>F)    
# Site                2 1.9657  0.9828  25.302 7.91e-08 ***
# Storage             1 0.0798  0.0798   2.056   0.1594    
# Scion               6 0.5131  0.0855   2.202   0.0628 .  
# Site:Storage        2 0.0768  0.0384   0.989   0.3809    
# Site:Scion         12 0.2677  0.0223   0.574   0.8494    
# Storage:Scion       6 0.1640  0.0273   0.704   0.6484    
# Site:Storage:Scion 12 0.1889  0.0157   0.405   0.9530    
# Residuals          40 1.5538  0.0388                     
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
abundance_results <- abundance_anova %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(abundance_results$Sum.Sq)
abundance_results$Perc.Var <- abundance_results$Sum.Sq / total_variance * 100

abundance_results %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> Sum.Sq </th>
   <th style="text-align:right;"> Mean.Sq </th>
   <th style="text-align:right;"> F.value </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1.9656600 </td>
   <td style="text-align:right;"> 0.9828300 </td>
   <td style="text-align:right;"> 25.3015018 </td>
   <td style="text-align:right;"> 0.0000001 </td>
   <td style="text-align:right;"> 40.867065 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0798455 </td>
   <td style="text-align:right;"> 0.0798455 </td>
   <td style="text-align:right;"> 2.0555037 </td>
   <td style="text-align:right;"> 0.1594283 </td>
   <td style="text-align:right;"> 1.660028 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.5131282 </td>
   <td style="text-align:right;"> 0.0855214 </td>
   <td style="text-align:right;"> 2.2016206 </td>
   <td style="text-align:right;"> 0.0628090 </td>
   <td style="text-align:right;"> 10.668194 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0768326 </td>
   <td style="text-align:right;"> 0.0384163 </td>
   <td style="text-align:right;"> 0.9889711 </td>
   <td style="text-align:right;"> 0.3808702 </td>
   <td style="text-align:right;"> 1.597389 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.2677207 </td>
   <td style="text-align:right;"> 0.0223101 </td>
   <td style="text-align:right;"> 0.5743394 </td>
   <td style="text-align:right;"> 0.8494197 </td>
   <td style="text-align:right;"> 5.566049 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.1639852 </td>
   <td style="text-align:right;"> 0.0273309 </td>
   <td style="text-align:right;"> 0.7035927 </td>
   <td style="text-align:right;"> 0.6484006 </td>
   <td style="text-align:right;"> 3.409335 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.1889265 </td>
   <td style="text-align:right;"> 0.0157439 </td>
   <td style="text-align:right;"> 0.4053027 </td>
   <td style="text-align:right;"> 0.9530223 </td>
   <td style="text-align:right;"> 3.927877 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 1.5537892 </td>
   <td style="text-align:right;"> 0.0388447 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 32.304063 </td>
  </tr>
</tbody>
</table>

### Communtiy size with canker count


```r
cat("Model formula: ", deparse(update(design_with_canker, log_copy_number ~ .)), "\n")
```

```
# Model formula:  log_copy_number ~ Site + Storage + Scion + Cankers + Site:Storage +      Site:Scion + Storage:Scion + Site:Cankers + Storage:Cankers +      Scion:Cankers + Site:Storage:Scion + Site:Storage:Cankers +      Site:Scion:Cankers + Storage:Scion:Cankers + Site:Storage:Scion:Cankers
```

```r
abundance_canker_anova <- aov(update(design_with_canker, log_copy_number ~ .), data = as.data.frame(colData(dds)))
summary(abundance_canker_anova)
```

```
#                            Df Sum Sq Mean Sq F value Pr(>F)  
# Site                        2 1.7693  0.8847  16.092 0.0585 .
# Storage                     1 0.0367  0.0367   0.667 0.4999  
# Scion                       6 0.5080  0.0847   1.540 0.4445  
# Cankers                     1 0.0000  0.0000   0.000 0.9986  
# Site:Storage                2 0.0996  0.0498   0.906 0.5246  
# Site:Scion                 12 0.2182  0.0182   0.331 0.9136  
# Storage:Scion               6 0.2595  0.0432   0.787 0.6535  
# Site:Cankers                2 0.0549  0.0275   0.500 0.6668  
# Storage:Cankers             1 0.0580  0.0580   1.056 0.4122  
# Scion:Cankers               6 0.2467  0.0411   0.748 0.6691  
# Site:Storage:Scion         11 0.2916  0.0265   0.482 0.8279  
# Site:Storage:Cankers        2 0.0546  0.0273   0.497 0.6682  
# Site:Scion:Cankers         11 0.2677  0.0243   0.443 0.8493  
# Storage:Scion:Cankers       5 0.3179  0.0636   1.156 0.5241  
# Site:Storage:Scion:Cankers  3 0.0725  0.0242   0.440 0.7494  
# Residuals                   2 0.1099  0.0550                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 8 observations deleted due to missingness
```

```r
abundance_canker_results <- abundance_canker_anova %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(abundance_canker_results$Sum.Sq)
abundance_canker_results$Perc.Var <- abundance_canker_results$Sum.Sq / total_variance * 100
abundance_canker_results %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> Sum.Sq </th>
   <th style="text-align:right;"> Mean.Sq </th>
   <th style="text-align:right;"> F.value </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1.7693200 </td>
   <td style="text-align:right;"> 0.8846600 </td>
   <td style="text-align:right;"> 16.0921702 </td>
   <td style="text-align:right;"> 0.0585063 </td>
   <td style="text-align:right;"> 40.5329702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0366603 </td>
   <td style="text-align:right;"> 0.0366603 </td>
   <td style="text-align:right;"> 0.6668591 </td>
   <td style="text-align:right;"> 0.4999459 </td>
   <td style="text-align:right;"> 0.8398426 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.5079574 </td>
   <td style="text-align:right;"> 0.0846596 </td>
   <td style="text-align:right;"> 1.5399771 </td>
   <td style="text-align:right;"> 0.4444623 </td>
   <td style="text-align:right;"> 11.6366865 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0000002 </td>
   <td style="text-align:right;"> 0.0000002 </td>
   <td style="text-align:right;"> 0.0000041 </td>
   <td style="text-align:right;"> 0.9985624 </td>
   <td style="text-align:right;"> 0.0000052 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0996205 </td>
   <td style="text-align:right;"> 0.0498102 </td>
   <td style="text-align:right;"> 0.9060598 </td>
   <td style="text-align:right;"> 0.5246425 </td>
   <td style="text-align:right;"> 2.2821840 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 0.2182080 </td>
   <td style="text-align:right;"> 0.0181840 </td>
   <td style="text-align:right;"> 0.3307712 </td>
   <td style="text-align:right;"> 0.9135563 </td>
   <td style="text-align:right;"> 4.9988806 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.2594567 </td>
   <td style="text-align:right;"> 0.0432428 </td>
   <td style="text-align:right;"> 0.7865963 </td>
   <td style="text-align:right;"> 0.6535158 </td>
   <td style="text-align:right;"> 5.9438380 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0549499 </td>
   <td style="text-align:right;"> 0.0274749 </td>
   <td style="text-align:right;"> 0.4997757 </td>
   <td style="text-align:right;"> 0.6667664 </td>
   <td style="text-align:right;"> 1.2588354 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0580426 </td>
   <td style="text-align:right;"> 0.0580426 </td>
   <td style="text-align:right;"> 1.0558078 </td>
   <td style="text-align:right;"> 0.4122002 </td>
   <td style="text-align:right;"> 1.3296847 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 0.2466772 </td>
   <td style="text-align:right;"> 0.0411129 </td>
   <td style="text-align:right;"> 0.7478526 </td>
   <td style="text-align:right;"> 0.6690618 </td>
   <td style="text-align:right;"> 5.6510750 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.2915615 </td>
   <td style="text-align:right;"> 0.0265056 </td>
   <td style="text-align:right;"> 0.4821428 </td>
   <td style="text-align:right;"> 0.8279391 </td>
   <td style="text-align:right;"> 6.6793190 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.0546069 </td>
   <td style="text-align:right;"> 0.0273034 </td>
   <td style="text-align:right;"> 0.4966559 </td>
   <td style="text-align:right;"> 0.6681562 </td>
   <td style="text-align:right;"> 1.2509773 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 0.2677051 </td>
   <td style="text-align:right;"> 0.0243368 </td>
   <td style="text-align:right;"> 0.4426926 </td>
   <td style="text-align:right;"> 0.8493071 </td>
   <td style="text-align:right;"> 6.1327994 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 0.3178814 </td>
   <td style="text-align:right;"> 0.0635763 </td>
   <td style="text-align:right;"> 1.1564672 </td>
   <td style="text-align:right;"> 0.5241354 </td>
   <td style="text-align:right;"> 7.2822761 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 0.0725410 </td>
   <td style="text-align:right;"> 0.0241803 </td>
   <td style="text-align:right;"> 0.4398457 </td>
   <td style="text-align:right;"> 0.7493798 </td>
   <td style="text-align:right;"> 1.6618254 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 0.1099491 </td>
   <td style="text-align:right;"> 0.0549746 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 2.5188007 </td>
  </tr>
</tbody>
</table>

## Alpha diversity analysis

### Alpha diversity plot


```r
# plot alpha diversity - plot_alpha will convert normalised abundances to integer values

bac_alpha_plot <- plot_alpha(
  counts(dds,normalize = F), colData(dds),
  design = "Scion", colour = "Site",
  measures = c("Shannon", "Simpson"),
  type="box"
) + 
  scale_colour_manual(values = cbPalette) + 
  theme(axis.title.x =  element_blank()) + 
  ggtitle("Bacterial α-diversity")

abundance_plot <- ggboxplot(
  data = as.data.frame(colData(dds)), x = "Site", y = "log_copy_number", 
  color = "Scion", add = "jitter", legend = "top", 
  title = "Bacterial abundance", xlab = "Site", ylab = "log10 copy number"
)

ggsave(
  filename = "bac_alpha.png", plot = bac_alpha_plot, path = "figures/", 
  height = 20, width = 40, units = "cm"
)

bac_alpha_plot
```

![](root_endophytes_files/figure-html/BAC α-diversity plot-1.png)<!-- -->

### Permutation based anova on α-diversity index ranks


```r
# get the diversity index data
all_alpha_ord <- plot_alpha(
  counts(dds, normalize = F), colData(dds), design = "Site", returnData = T
)

# join diversity indices and metadata
all_alpha_ord <- all_alpha_ord[
  as.data.table(colData(dds), keep.rownames = "Samples"), on = "Samples"
]

bac_alpha <- all_alpha_ord

formula <- FULL_DESIGN

# Chao1

setkey(all_alpha_ord, S.chao1)
all_alpha_ord[, measure := as.numeric(as.factor(S.chao1))]
result <- aovp(update(formula, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS "
```

```r
summary(result)
```

```
# Component 1 :
#                    Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# Site                2  20914.4   10457.2 5000   <2e-16 ***
# Storage             1    113.5     113.5  110   0.4818    
# Site:Storage        2   3319.8    1659.9 5000   0.0094 ** 
# Scion               6   1291.7     215.3  407   0.5037    
# Site:Scion         12   5425.8     452.1 2201   0.1577    
# Storage:Scion       6    867.4     144.6  699   0.8054    
# Site:Storage:Scion 12   1912.4     159.4  551   0.9456    
# Residuals          40  12095.5     302.4                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 20914.3677 </td>
   <td style="text-align:right;"> 10457.1839 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 45.5249023 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 113.5359 </td>
   <td style="text-align:right;"> 113.5359 </td>
   <td style="text-align:right;"> 110 </td>
   <td style="text-align:right;"> 0.4818182 </td>
   <td style="text-align:right;"> 0.2471368 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3319.8162 </td>
   <td style="text-align:right;"> 1659.9081 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0094000 </td>
   <td style="text-align:right;"> 7.2263388 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1291.7022 </td>
   <td style="text-align:right;"> 215.2837 </td>
   <td style="text-align:right;"> 407 </td>
   <td style="text-align:right;"> 0.5036855 </td>
   <td style="text-align:right;"> 2.8116852 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 5425.7655 </td>
   <td style="text-align:right;"> 452.1471 </td>
   <td style="text-align:right;"> 2201 </td>
   <td style="text-align:right;"> 0.1576556 </td>
   <td style="text-align:right;"> 11.8104189 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 867.3693 </td>
   <td style="text-align:right;"> 144.5615 </td>
   <td style="text-align:right;"> 699 </td>
   <td style="text-align:right;"> 0.8054363 </td>
   <td style="text-align:right;"> 1.8880275 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 1912.4432 </td>
   <td style="text-align:right;"> 159.3703 </td>
   <td style="text-align:right;"> 551 </td>
   <td style="text-align:right;"> 0.9455535 </td>
   <td style="text-align:right;"> 4.1628699 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 12095.5000 </td>
   <td style="text-align:right;"> 302.3875 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 26.3286207 </td>
  </tr>
</tbody>
</table>

```r
# Shannon

setkey(all_alpha_ord, shannon)
all_alpha_ord[, measure := as.numeric(as.factor(shannon))]
result <- aovp(update(formula, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS "
```

```r
summary(result)
```

```
# Component 1 :
#                    Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# Site                2  23583.2   11791.6 5000   <2e-16 ***
# Storage             1     89.3      89.3  120   0.4583    
# Site:Storage        2    898.4     449.2  714   0.2983    
# Scion               6    994.0     165.7  165   0.7939    
# Site:Scion         12   3435.7     286.3 5000   0.7100    
# Storage:Scion       6    481.5      80.3   84   1.0000    
# Site:Storage:Scion 12   2143.3     178.6 2635   0.9450    
# Residuals          40  14315.0     357.9                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 23583.23942 </td>
   <td style="text-align:right;"> 11791.61971 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 51.3343116 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 89.26865 </td>
   <td style="text-align:right;"> 89.26865 </td>
   <td style="text-align:right;"> 120 </td>
   <td style="text-align:right;"> 0.4583333 </td>
   <td style="text-align:right;"> 0.1943136 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 898.39852 </td>
   <td style="text-align:right;"> 449.19926 </td>
   <td style="text-align:right;"> 714 </td>
   <td style="text-align:right;"> 0.2983193 </td>
   <td style="text-align:right;"> 1.9555698 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 994.03480 </td>
   <td style="text-align:right;"> 165.67247 </td>
   <td style="text-align:right;"> 165 </td>
   <td style="text-align:right;"> 0.7939394 </td>
   <td style="text-align:right;"> 2.1637440 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3435.70147 </td>
   <td style="text-align:right;"> 286.30846 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.7100000 </td>
   <td style="text-align:right;"> 7.4785896 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 481.54603 </td>
   <td style="text-align:right;"> 80.25767 </td>
   <td style="text-align:right;"> 84 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 1.0481950 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2143.31111 </td>
   <td style="text-align:right;"> 178.60926 </td>
   <td style="text-align:right;"> 2635 </td>
   <td style="text-align:right;"> 0.9449715 </td>
   <td style="text-align:right;"> 4.6654066 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 14315.00000 </td>
   <td style="text-align:right;"> 357.87500 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 31.1598698 </td>
  </tr>
</tbody>
</table>

```r
# Simpson

setkey(all_alpha_ord, simpson)
all_alpha_ord[, measure := as.numeric(as.factor(simpson))]
result <- aovp(update(formula, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS "
```

```r
summary(result)
```

```
# Component 1 :
#                    Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# Site                2  20851.3   10425.7 5000   <2e-16 ***
# Storage             1     31.0      31.0  163   0.3804    
# Site:Storage        2   1422.3     711.2 1690   0.1118    
# Scion               6   1995.8     332.6  609   0.3924    
# Site:Scion         12   3933.8     327.8 1536   0.4915    
# Storage:Scion       6   1003.5     167.3  677   0.5968    
# Site:Storage:Scion 12   3012.3     251.0 3909   0.7135    
# Residuals          40  13690.5     342.3                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 20851.30423 </td>
   <td style="text-align:right;"> 10425.65212 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 45.3876302 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 30.96912 </td>
   <td style="text-align:right;"> 30.96912 </td>
   <td style="text-align:right;"> 163 </td>
   <td style="text-align:right;"> 0.3803681 </td>
   <td style="text-align:right;"> 0.0674114 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1422.32005 </td>
   <td style="text-align:right;"> 711.16002 </td>
   <td style="text-align:right;"> 1690 </td>
   <td style="text-align:right;"> 0.1118343 </td>
   <td style="text-align:right;"> 3.0960047 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1995.82345 </td>
   <td style="text-align:right;"> 332.63724 </td>
   <td style="text-align:right;"> 609 </td>
   <td style="text-align:right;"> 0.3924466 </td>
   <td style="text-align:right;"> 4.3443660 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3933.80635 </td>
   <td style="text-align:right;"> 327.81720 </td>
   <td style="text-align:right;"> 1536 </td>
   <td style="text-align:right;"> 0.4915365 </td>
   <td style="text-align:right;"> 8.5628288 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1003.51689 </td>
   <td style="text-align:right;"> 167.25281 </td>
   <td style="text-align:right;"> 677 </td>
   <td style="text-align:right;"> 0.5967504 </td>
   <td style="text-align:right;"> 2.1843839 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3012.25990 </td>
   <td style="text-align:right;"> 251.02166 </td>
   <td style="text-align:right;"> 3909 </td>
   <td style="text-align:right;"> 0.7134817 </td>
   <td style="text-align:right;"> 6.5568722 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 13690.50000 </td>
   <td style="text-align:right;"> 342.26250 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 29.8005028 </td>
  </tr>
</tbody>
</table>

#### Permutation based anova on α-diversity index ranks with canker lesion counts


```r
cat("Model formula: ", deparse(design_with_canker), "\n")
```

```
# Model formula:  y ~ Site * Storage * Scion * Cankers
```

```r
# Chao1

setkey(all_alpha_ord, S.chao1)
all_alpha_ord[, measure := as.numeric(as.factor(S.chao1))]
result <- aovp(update(design_with_canker, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS : numeric variables centered"
```

```r
summary(result)
```

```
# Component 1 :
#                            Df R Sum Sq R Mean Sq Iter Pr(Prob)    
# Site                        2  17286.5    8643.2 5000  < 2e-16 ***
# Storage                     1     25.7      25.7  159  0.38994    
# Site:Storage                2   2459.2    1229.6 5000  0.00860 ** 
# Scion                       6    959.7     160.0 2258  0.09655 .  
# Site:Scion                 12   5232.0     436.0 5000  0.03540 *  
# Storage:Scion               6   1440.7     240.1 5000  0.08520 .  
# Site:Storage:Scion         11   1947.9     177.1 3371  0.09641 .  
# Cankers                     1      9.1       9.1  146  0.41096    
# Site:Cankers                2   1346.0     673.0 5000  0.02440 *  
# Storage:Cankers             1    139.8     139.8 1182  0.07868 .  
# Site:Storage:Cankers        2    384.4     192.2 1589  0.07992 .  
# Scion:Cankers               6   2396.5     399.4 5000  0.03240 *  
# Site:Scion:Cankers         11   3673.2     333.9 5000  0.04820 *  
# Storage:Scion:Cankers       5   1630.3     326.1 5000  0.03500 *  
# Site:Storage:Scion:Cankers  3    783.7     261.2 3717  0.08475 .  
# Residuals                   2     36.5      18.3                  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 17286.462354 </td>
   <td style="text-align:right;"> 8643.231177 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 43.4867581 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 25.691923 </td>
   <td style="text-align:right;"> 25.691923 </td>
   <td style="text-align:right;"> 159 </td>
   <td style="text-align:right;"> 0.3899371 </td>
   <td style="text-align:right;"> 0.0646320 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2459.189669 </td>
   <td style="text-align:right;"> 1229.594834 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0086000 </td>
   <td style="text-align:right;"> 6.1864703 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 959.749728 </td>
   <td style="text-align:right;"> 159.958288 </td>
   <td style="text-align:right;"> 2258 </td>
   <td style="text-align:right;"> 0.0965456 </td>
   <td style="text-align:right;"> 2.4143982 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 5231.991372 </td>
   <td style="text-align:right;"> 435.999281 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0354000 </td>
   <td style="text-align:right;"> 13.1618800 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1440.657689 </td>
   <td style="text-align:right;"> 240.109615 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0852000 </td>
   <td style="text-align:right;"> 3.6241963 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 1947.851861 </td>
   <td style="text-align:right;"> 177.077442 </td>
   <td style="text-align:right;"> 3371 </td>
   <td style="text-align:right;"> 0.0964106 </td>
   <td style="text-align:right;"> 4.9001213 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9.058379 </td>
   <td style="text-align:right;"> 9.058379 </td>
   <td style="text-align:right;"> 146 </td>
   <td style="text-align:right;"> 0.4109589 </td>
   <td style="text-align:right;"> 0.0227877 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 1346.041139 </td>
   <td style="text-align:right;"> 673.020569 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0244000 </td>
   <td style="text-align:right;"> 3.3861738 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 139.829164 </td>
   <td style="text-align:right;"> 139.829164 </td>
   <td style="text-align:right;"> 1182 </td>
   <td style="text-align:right;"> 0.0786802 </td>
   <td style="text-align:right;"> 0.3517618 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 384.389932 </td>
   <td style="text-align:right;"> 192.194966 </td>
   <td style="text-align:right;"> 1589 </td>
   <td style="text-align:right;"> 0.0799245 </td>
   <td style="text-align:right;"> 0.9669921 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2396.539956 </td>
   <td style="text-align:right;"> 399.423326 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0324000 </td>
   <td style="text-align:right;"> 6.0288653 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 3673.152926 </td>
   <td style="text-align:right;"> 333.922993 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0482000 </td>
   <td style="text-align:right;"> 9.2403818 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1630.262608 </td>
   <td style="text-align:right;"> 326.052522 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0350000 </td>
   <td style="text-align:right;"> 4.1011766 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 783.725895 </td>
   <td style="text-align:right;"> 261.241965 </td>
   <td style="text-align:right;"> 3717 </td>
   <td style="text-align:right;"> 0.0847458 </td>
   <td style="text-align:right;"> 1.9715832 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 36.500000 </td>
   <td style="text-align:right;"> 18.250000 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.0918214 </td>
  </tr>
</tbody>
</table>

```r
# Shannon

setkey(all_alpha_ord, shannon)
all_alpha_ord[, measure := as.numeric(as.factor(shannon))]
result <- aovp(update(design_with_canker, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS : numeric variables centered"
```

```r
summary(result)
```

```
# Component 1 :
#                            Df R Sum Sq R Mean Sq Iter Pr(Prob)   
# Site                        2  21716.8   10858.4 5000  0.00500 **
# Storage                     1    144.4     144.4   51  1.00000   
# Site:Storage                2    536.9     268.5   51  0.84314   
# Scion                       6   1135.6     189.3 1279  0.39484   
# Site:Scion                 12   3509.4     292.4 1280  0.24453   
# Storage:Scion               6   1258.1     209.7  424  0.54009   
# Site:Storage:Scion         11   1737.3     157.9 1255  0.46056   
# Cankers                     1    240.9     240.9   51  1.00000   
# Site:Cankers                2   3168.1    1584.0 1314  0.09513 . 
# Storage:Cankers             1    222.3     222.3   51  0.78431   
# Site:Storage:Cankers        2    336.6     168.3   53  1.00000   
# Scion:Cankers               6   1314.3     219.1  363  0.49587   
# Site:Scion:Cankers         11   4569.1     415.4 1302  0.20430   
# Storage:Scion:Cankers       5   1819.4     363.9  759  0.32148   
# Site:Storage:Scion:Cankers  3    186.3      62.1  231  0.65368   
# Residuals                   2    153.0      76.5                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 21716.7551 </td>
   <td style="text-align:right;"> 10858.37753 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0050000 </td>
   <td style="text-align:right;"> 51.6469197 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 144.3662 </td>
   <td style="text-align:right;"> 144.36619 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.3433325 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 536.9468 </td>
   <td style="text-align:right;"> 268.47340 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 0.8431373 </td>
   <td style="text-align:right;"> 1.2769702 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1135.6281 </td>
   <td style="text-align:right;"> 189.27134 </td>
   <td style="text-align:right;"> 1279 </td>
   <td style="text-align:right;"> 0.3948397 </td>
   <td style="text-align:right;"> 2.7007576 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3509.3930 </td>
   <td style="text-align:right;"> 292.44942 </td>
   <td style="text-align:right;"> 1280 </td>
   <td style="text-align:right;"> 0.2445312 </td>
   <td style="text-align:right;"> 8.3460600 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1258.0868 </td>
   <td style="text-align:right;"> 209.68113 </td>
   <td style="text-align:right;"> 424 </td>
   <td style="text-align:right;"> 0.5400943 </td>
   <td style="text-align:right;"> 2.9919897 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 1737.3241 </td>
   <td style="text-align:right;"> 157.93855 </td>
   <td style="text-align:right;"> 1255 </td>
   <td style="text-align:right;"> 0.4605578 </td>
   <td style="text-align:right;"> 4.1317147 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 240.9417 </td>
   <td style="text-align:right;"> 240.94166 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.5730089 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3168.0618 </td>
   <td style="text-align:right;"> 1584.03091 </td>
   <td style="text-align:right;"> 1314 </td>
   <td style="text-align:right;"> 0.0951294 </td>
   <td style="text-align:right;"> 7.5343040 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 222.2887 </td>
   <td style="text-align:right;"> 222.28869 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 0.7843137 </td>
   <td style="text-align:right;"> 0.5286483 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 336.6340 </td>
   <td style="text-align:right;"> 168.31699 </td>
   <td style="text-align:right;"> 53 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> 0.8005850 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1314.3344 </td>
   <td style="text-align:right;"> 219.05573 </td>
   <td style="text-align:right;"> 363 </td>
   <td style="text-align:right;"> 0.4958678 </td>
   <td style="text-align:right;"> 3.1257581 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 4569.0691 </td>
   <td style="text-align:right;"> 415.36992 </td>
   <td style="text-align:right;"> 1302 </td>
   <td style="text-align:right;"> 0.2043011 </td>
   <td style="text-align:right;"> 10.8661882 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1819.4173 </td>
   <td style="text-align:right;"> 363.88346 </td>
   <td style="text-align:right;"> 759 </td>
   <td style="text-align:right;"> 0.3214756 </td>
   <td style="text-align:right;"> 4.3269493 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 186.2531 </td>
   <td style="text-align:right;"> 62.08436 </td>
   <td style="text-align:right;"> 231 </td>
   <td style="text-align:right;"> 0.6536797 </td>
   <td style="text-align:right;"> 0.4429482 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 153.0000 </td>
   <td style="text-align:right;"> 76.50000 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.3638655 </td>
  </tr>
</tbody>
</table>

```r
# Simpson

setkey(all_alpha_ord, simpson)
all_alpha_ord[, measure := as.numeric(as.factor(simpson))]
result <- aovp(update(design_with_canker, measure ~ .), all_alpha_ord, seqs = T)
```

```
# [1] "Settings:  sequential SS : numeric variables centered"
```

```r
summary(result)
```

```
# Component 1 :
#                            Df R Sum Sq R Mean Sq Iter Pr(Prob)   
# Site                        2  20046.3   10023.1 5000  0.00260 **
# Storage                     1     16.9      16.9   51  0.88235   
# Site:Storage                2    915.5     457.7 3201  0.06061 . 
# Scion                       6   2312.5     385.4 2714  0.10022   
# Site:Scion                 12   3655.7     304.6 3955  0.11757   
# Storage:Scion               6   1149.8     191.6 1224  0.19444   
# Site:Storage:Scion         11   3055.0     277.7 3576  0.11242   
# Cankers                     1    147.7     147.7   87  0.54023   
# Site:Cankers                2   3560.6    1780.3 5000  0.01140 * 
# Storage:Cankers             1    174.1     174.1 1030  0.08932 . 
# Site:Storage:Cankers        2    570.5     285.3 2377  0.10265   
# Scion:Cankers               6   1465.7     244.3 1170  0.13590   
# Site:Scion:Cankers         11   4620.8     420.1 3963  0.05905 . 
# Storage:Scion:Cankers       5   1226.8     245.4 2716  0.15869   
# Site:Storage:Scion:Cankers  3    116.3      38.8   77  0.83117   
# Residuals                   2    116.5      58.2                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% summary() %>% unclass() %>% data.frame()
total_variance <- sum(df$R.Sum.Sq)
df$Perc.Var <- df$R.Sum.Sq / total_variance * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> R.Sum.Sq </th>
   <th style="text-align:right;"> R.Mean.Sq </th>
   <th style="text-align:right;"> Iter </th>
   <th style="text-align:right;"> Pr.Prob. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 20046.25092 </td>
   <td style="text-align:right;"> 10023.12546 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0026000 </td>
   <td style="text-align:right;"> 46.4564155 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 16.91887 </td>
   <td style="text-align:right;"> 16.91887 </td>
   <td style="text-align:right;"> 51 </td>
   <td style="text-align:right;"> 0.8823529 </td>
   <td style="text-align:right;"> 0.0392088 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 915.49527 </td>
   <td style="text-align:right;"> 457.74763 </td>
   <td style="text-align:right;"> 3201 </td>
   <td style="text-align:right;"> 0.0606061 </td>
   <td style="text-align:right;"> 2.1216251 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2312.47791 </td>
   <td style="text-align:right;"> 385.41299 </td>
   <td style="text-align:right;"> 2714 </td>
   <td style="text-align:right;"> 0.1002211 </td>
   <td style="text-align:right;"> 5.3590786 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 3655.69909 </td>
   <td style="text-align:right;"> 304.64159 </td>
   <td style="text-align:right;"> 3955 </td>
   <td style="text-align:right;"> 0.1175727 </td>
   <td style="text-align:right;"> 8.4719420 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1149.79543 </td>
   <td style="text-align:right;"> 191.63257 </td>
   <td style="text-align:right;"> 1224 </td>
   <td style="text-align:right;"> 0.1944444 </td>
   <td style="text-align:right;"> 2.6646067 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 3055.02468 </td>
   <td style="text-align:right;"> 277.72952 </td>
   <td style="text-align:right;"> 3576 </td>
   <td style="text-align:right;"> 0.1124161 </td>
   <td style="text-align:right;"> 7.0799022 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 147.74792 </td>
   <td style="text-align:right;"> 147.74792 </td>
   <td style="text-align:right;"> 87 </td>
   <td style="text-align:right;"> 0.5402299 </td>
   <td style="text-align:right;"> 0.3424001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3560.56456 </td>
   <td style="text-align:right;"> 1780.28228 </td>
   <td style="text-align:right;"> 5000 </td>
   <td style="text-align:right;"> 0.0114000 </td>
   <td style="text-align:right;"> 8.2514714 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 174.11409 </td>
   <td style="text-align:right;"> 174.11409 </td>
   <td style="text-align:right;"> 1030 </td>
   <td style="text-align:right;"> 0.0893204 </td>
   <td style="text-align:right;"> 0.4035027 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 570.51731 </td>
   <td style="text-align:right;"> 285.25866 </td>
   <td style="text-align:right;"> 2377 </td>
   <td style="text-align:right;"> 0.1026504 </td>
   <td style="text-align:right;"> 1.3221519 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 1465.71745 </td>
   <td style="text-align:right;"> 244.28624 </td>
   <td style="text-align:right;"> 1170 </td>
   <td style="text-align:right;"> 0.1358974 </td>
   <td style="text-align:right;"> 3.3967438 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 4620.83118 </td>
   <td style="text-align:right;"> 420.07556 </td>
   <td style="text-align:right;"> 3963 </td>
   <td style="text-align:right;"> 0.0590462 </td>
   <td style="text-align:right;"> 10.7085986 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 1226.75147 </td>
   <td style="text-align:right;"> 245.35029 </td>
   <td style="text-align:right;"> 2716 </td>
   <td style="text-align:right;"> 0.1586892 </td>
   <td style="text-align:right;"> 2.8429493 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 116.25601 </td>
   <td style="text-align:right;"> 38.75200 </td>
   <td style="text-align:right;"> 77 </td>
   <td style="text-align:right;"> 0.8311688 </td>
   <td style="text-align:right;"> 0.2694188 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residuals </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 116.50000 </td>
   <td style="text-align:right;"> 58.25000 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 0.2699843 </td>
  </tr>
</tbody>
</table>

## Beta diversity PCA/NMDS

### PCA 

```r
# Number of PCs to include
n_pcs <- 10

# perform PC decomposition of DES object
mypca <- des_to_pca(dds)

# to get pca plot axis into the same scale create a dataframe of PC scores multiplied by their variance
bac_pca <- t(data.frame(t(mypca$x) * mypca$percentVar))

formula = FULL_DESIGN
```

#### Percent variation in first 10 PCs 

```r
# Cumulative percentage of variance explained
pca_cum_var <- data.frame(
  cumulative = cumsum(mypca$percentVar * 100),
  no = 1:length(mypca$percentVar)
)

# Plot cumulative percentage of variance explained
bac_cum_pca <- ggline(
  pca_cum_var, x = "no", y = "cumulative", plot_type = "l",
  xlab = "Number of PCs", ylab = "Cumulative % variance explained",
  title = "Bacteria: cumulative % variance explained by PCs"
)
ggsave(filename = "bac_cum_pca.png", plot = bac_cum_pca, path = "figures/",)
bac_cum_pca
```

![](root_endophytes_files/figure-html/BAC PCA var-1.png)<!-- -->

```r
pca_var <- data.frame(
  PC = paste0("PC", 1:n_pcs),
  perc_var = round(mypca$percentVar[1:n_pcs] * 100, 1)
)

pca_var %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:right;"> perc_var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:right;"> 18.6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:right;"> 12.7 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:right;"> 7.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:right;"> 4.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:right;"> 3.0 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:right;"> 2.1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:right;"> 1.8 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:right;"> 1.6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:right;"> 1.3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:right;"> 1.3 </td>
  </tr>
</tbody>
</table>

#### ANOVA of first 10 PCs

```r
pca_summary <- apply(
  mypca$x[, 1:n_pcs], 2, 
  function(x){
    summary(aov(update(formula, x ~ .), data = as.data.frame(cbind(x, colData(dds)))))
  }
)

pca_summary
```

```
# $PC1
#                    Df Sum Sq Mean Sq F value Pr(>F)    
# Site                2 190646   95323 224.072 <2e-16 ***
# Storage             1     25      25   0.058 0.8106    
# Scion               6   4227     704   1.656 0.1572    
# Site:Storage        2   2566    1283   3.015 0.0603 .  
# Site:Scion         12   8401     700   1.646 0.1176    
# Storage:Scion       6   1527     255   0.598 0.7298    
# Site:Storage:Scion 12   5124     427   1.004 0.4634    
# Residuals          40  17016     425                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC2
#                    Df Sum Sq Mean Sq F value Pr(>F)    
# Site                2 137668   68834 276.644 <2e-16 ***
# Storage             1    142     142   0.572  0.454    
# Scion               6   1474     246   0.987  0.447    
# Site:Storage        2    604     302   1.213  0.308    
# Site:Scion         12   4425     369   1.482  0.172    
# Storage:Scion       6    711     119   0.476  0.822    
# Site:Storage:Scion 12   2035     170   0.681  0.759    
# Residuals          40   9953     249                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC3
#                    Df Sum Sq Mean Sq F value   Pr(>F)    
# Site                2   2673    1337   1.562 0.222149    
# Storage             1     30      30   0.035 0.851920    
# Scion               6  13665    2278   2.663 0.028714 *  
# Site:Storage        2  14748    7374   8.621 0.000771 ***
# Site:Scion         12  10740     895   1.046 0.428331    
# Storage:Scion       6   8140    1357   1.586 0.176467    
# Site:Storage:Scion 12   2870     239   0.280 0.989466    
# Residuals          40  34217     855                     
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC4
#                    Df Sum Sq Mean Sq F value  Pr(>F)   
# Site                2   5281  2640.3   5.834 0.00598 **
# Storage             1    635   634.6   1.402 0.24336   
# Scion               6   3161   526.9   1.164 0.34456   
# Site:Storage        2    630   315.1   0.696 0.50442   
# Site:Scion         12  13692  1141.0   2.521 0.01418 * 
# Storage:Scion       6   1827   304.5   0.673 0.67221   
# Site:Storage:Scion 12   6439   536.6   1.186 0.32589   
# Residuals          40  18104   452.6                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC5
#                    Df Sum Sq Mean Sq F value  Pr(>F)   
# Site                2      5       2   0.006 0.99447   
# Storage             1   4033    4033   9.433 0.00382 **
# Scion               6   1458     243   0.568 0.75299   
# Site:Storage        2   5967    2984   6.979 0.00251 **
# Site:Scion         12   2582     215   0.503 0.90010   
# Storage:Scion       6   1984     331   0.773 0.59554   
# Site:Storage:Scion 12   3480     290   0.678 0.76149   
# Residuals          40  17102     428                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC6
#                    Df Sum Sq Mean Sq F value Pr(>F)  
# Site                2     32    16.2   0.055 0.9467  
# Storage             1   1498  1497.8   5.070 0.0299 *
# Scion               6   4754   792.3   2.682 0.0278 *
# Site:Storage        2    811   405.3   1.372 0.2653  
# Site:Scion         12   2040   170.0   0.575 0.8486  
# Storage:Scion       6   1724   287.4   0.973 0.4560  
# Site:Storage:Scion 12   3166   263.8   0.893 0.5611  
# Residuals          40  11817   295.4                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC7
#                    Df Sum Sq Mean Sq F value Pr(>F)  
# Site                2     85    42.7   0.157 0.8556  
# Storage             1    448   448.2   1.645 0.2070  
# Scion               6   2285   380.9   1.398 0.2393  
# Site:Storage        2   2000   999.9   3.670 0.0344 *
# Site:Scion         12   2462   205.1   0.753 0.6925  
# Storage:Scion       6    983   163.9   0.602 0.7273  
# Site:Storage:Scion 12   3081   256.7   0.942 0.5166  
# Residuals          40  10898   272.4                 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC8
#                    Df Sum Sq Mean Sq F value Pr(>F)
# Site                2    110    55.1   0.197  0.822
# Storage             1    549   548.5   1.963  0.169
# Scion               6    750   125.1   0.448  0.842
# Site:Storage        2    334   167.1   0.598  0.555
# Site:Scion         12   1153    96.1   0.344  0.975
# Storage:Scion       6   1856   309.4   1.108  0.375
# Site:Storage:Scion 12   3233   269.4   0.964  0.497
# Residuals          40  11175   279.4               
# 
# $PC9
#                    Df Sum Sq Mean Sq F value  Pr(>F)   
# Site                2    291   145.4   0.624 0.54101   
# Storage             1    120   120.0   0.515 0.47731   
# Scion               6    658   109.6   0.470 0.82619   
# Site:Storage        2   3474  1737.2   7.452 0.00177 **
# Site:Scion         12   1215   101.2   0.434 0.93963   
# Storage:Scion       6    606   101.0   0.433 0.85225   
# Site:Storage:Scion 12    970    80.9   0.347 0.97417   
# Residuals          40   9324   233.1                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# $PC10
#                    Df Sum Sq Mean Sq F value  Pr(>F)   
# Site                2    131    65.4   0.396 0.67572   
# Storage             1    365   365.5   2.213 0.14466   
# Scion               6   1060   176.7   1.070 0.39635   
# Site:Storage        2   2202  1100.8   6.666 0.00317 **
# Site:Scion         12   2043   170.2   1.031 0.44086   
# Storage:Scion       6   1947   324.4   1.965 0.09381 . 
# Site:Storage:Scion 12   1980   165.0   0.999 0.46711   
# Residuals          40   6605   165.1                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

#### Percent variation in first 10 PCs for each factor


```r
# Extract PC scores as a list of dataframes
pcas <- lapply(pca_summary, function(i) data.frame(unclass(i)))

# Merge into single dataframe
pcs_factors_tidy <- lapply(
  names(pcas),
  function(name) {
    pcas[[name]] %>%
    mutate(
      PC = name, #substring(name, 3),
      Factor = gsub(" ", "", rownames(pcas[[name]])),
      var = Sum.Sq / sum(pcas[[name]]$Sum.Sq) * 100,
      pc_var = subset(pca_var, PC == name)$"perc_var",
      total_var = var * pc_var / 100,
      sig = case_when(
        is.na(Pr..F.) ~ "",
        Pr..F. < 0.001 ~ "***",
        Pr..F. < 0.01 ~ "**",
        Pr..F. < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variance = ifelse(
        total_var < 0.01, paste0("<0.01", sig),
        paste0(round(total_var, 2), sig)
      )
    )
  }
) %>% bind_rows() %>% data.table()

# Order PCs and factors
pcs_factors_tidy$PC <- factor(pcs_factors_tidy$PC, levels = paste0("PC", 1:n_pcs))
pcs_factors_tidy$Factor <- factor(pcs_factors_tidy$Factor, levels = unique(pcs_factors_tidy$Factor))

# Significant factors
pcs_factors_tidy[
  Pr..F. < 0.05, 
  c("PC", "Factor", "Df", "F.value", "Pr..F.", "var", "pc_var", "total_var")
] %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:left;"> Factor </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> F.value </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> var </th>
   <th style="text-align:right;"> pc_var </th>
   <th style="text-align:right;"> total_var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 224.072105 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 83.058689 </td>
   <td style="text-align:right;"> 18.6 </td>
   <td style="text-align:right;"> 15.4489162 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 276.643717 </td>
   <td style="text-align:right;"> 0.0000000 </td>
   <td style="text-align:right;"> 87.680340 </td>
   <td style="text-align:right;"> 12.7 </td>
   <td style="text-align:right;"> 11.1354032 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2.662540 </td>
   <td style="text-align:right;"> 0.0287140 </td>
   <td style="text-align:right;"> 15.692236 </td>
   <td style="text-align:right;"> 7.1 </td>
   <td style="text-align:right;"> 1.1141487 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 8.620505 </td>
   <td style="text-align:right;"> 0.0007710 </td>
   <td style="text-align:right;"> 16.935585 </td>
   <td style="text-align:right;"> 7.1 </td>
   <td style="text-align:right;"> 1.2024266 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 5.833732 </td>
   <td style="text-align:right;"> 0.0059821 </td>
   <td style="text-align:right;"> 10.610344 </td>
   <td style="text-align:right;"> 4.0 </td>
   <td style="text-align:right;"> 0.4244137 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 2.521003 </td>
   <td style="text-align:right;"> 0.0141754 </td>
   <td style="text-align:right;"> 27.511079 </td>
   <td style="text-align:right;"> 4.0 </td>
   <td style="text-align:right;"> 1.1004432 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 9.432657 </td>
   <td style="text-align:right;"> 0.0038232 </td>
   <td style="text-align:right;"> 11.015538 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.3304662 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 6.978586 </td>
   <td style="text-align:right;"> 0.0025131 </td>
   <td style="text-align:right;"> 16.299308 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.4889792 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 5.069838 </td>
   <td style="text-align:right;"> 0.0299066 </td>
   <td style="text-align:right;"> 5.795920 </td>
   <td style="text-align:right;"> 2.1 </td>
   <td style="text-align:right;"> 0.1217143 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 2.681842 </td>
   <td style="text-align:right;"> 0.0277905 </td>
   <td style="text-align:right;"> 18.395545 </td>
   <td style="text-align:right;"> 2.1 </td>
   <td style="text-align:right;"> 0.3863064 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 3.670184 </td>
   <td style="text-align:right;"> 0.0344007 </td>
   <td style="text-align:right;"> 8.991247 </td>
   <td style="text-align:right;"> 1.8 </td>
   <td style="text-align:right;"> 0.1618424 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 7.452304 </td>
   <td style="text-align:right;"> 0.0017743 </td>
   <td style="text-align:right;"> 20.856556 </td>
   <td style="text-align:right;"> 1.3 </td>
   <td style="text-align:right;"> 0.2711352 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 6.666103 </td>
   <td style="text-align:right;"> 0.0031726 </td>
   <td style="text-align:right;"> 13.479337 </td>
   <td style="text-align:right;"> 1.3 </td>
   <td style="text-align:right;"> 0.1752314 </td>
  </tr>
</tbody>
</table>

```r
# Table with factors as columns and PCs as rows
# pcs_factors <- dcast(pcs_factors_tidy, PC ~ Factor, value.var = "variance")
pcs_factors <- pcs_factors_tidy %>%
  select(PC, pc_var, Factor, variance) %>%
  spread(key = Factor, value = variance)

pcs_factors %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:right;"> pc_var </th>
   <th style="text-align:left;"> Site </th>
   <th style="text-align:left;"> Storage </th>
   <th style="text-align:left;"> Scion </th>
   <th style="text-align:left;"> Site:Storage </th>
   <th style="text-align:left;"> Site:Scion </th>
   <th style="text-align:left;"> Storage:Scion </th>
   <th style="text-align:left;"> Site:Storage:Scion </th>
   <th style="text-align:left;"> Residuals </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:right;"> 18.6 </td>
   <td style="text-align:left;"> 15.45*** </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.34 </td>
   <td style="text-align:left;"> 0.21 </td>
   <td style="text-align:left;"> 0.68 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.42 </td>
   <td style="text-align:left;"> 1.38 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:right;"> 12.7 </td>
   <td style="text-align:left;"> 11.14*** </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.36 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.81 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:right;"> 7.1 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 1.11* </td>
   <td style="text-align:left;"> 1.2*** </td>
   <td style="text-align:left;"> 0.88 </td>
   <td style="text-align:left;"> 0.66 </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> 2.79 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:right;"> 4.0 </td>
   <td style="text-align:left;"> 0.42** </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 1.1* </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.52 </td>
   <td style="text-align:left;"> 1.46 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.33** </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.49** </td>
   <td style="text-align:left;"> 0.21 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.29 </td>
   <td style="text-align:left;"> 1.4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:right;"> 2.1 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.12* </td>
   <td style="text-align:left;"> 0.39* </td>
   <td style="text-align:left;"> 0.07 </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.26 </td>
   <td style="text-align:left;"> 0.96 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:right;"> 1.8 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> 0.16* </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> 0.88 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:right;"> 1.6 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.27 </td>
   <td style="text-align:left;"> 0.93 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:right;"> 1.3 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.27** </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.73 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:right;"> 1.3 </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.18** </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.53 </td>
  </tr>
</tbody>
</table>

#### PCA plot

```r
bac_pca_plot <- plotOrd(
  bac_pca,
  colData(dds),
  design = "Site",
  shape = "Storage",
  axes = c(1, 2),
  # facet = "Storage", 
  cbPalette = T,
  alpha = 0.75,
) #+ facet_wrap(~facet)

ggsave(filename = "bac_pca_plot.png", plot = bac_pca_plot, path = "figures/")

bac_pca_plot
```

![](root_endophytes_files/figure-html/BAC PCA plot-1.png)<!-- -->

```r
bac_pca_3_6_plot <- plotOrd(
  bac_pca,
  colData(dds),
  design = "Scion",
  shape = "Storage",
  axes = c(3, 6), 
  cbPalette = T,
  alpha = 0.75,
)

ggsave(filename = "bac_pca_3_6_plot.png", plot = bac_pca_3_6_plot, path = "figures/")

bac_pca_3_6_plot
```

![](root_endophytes_files/figure-html/BAC PCA plot-2.png)<!-- -->

#### PCA sum of squares (% var)


```r
sum_squares <- apply(mypca$x, 2 ,function(x) 
  summary(aov(update(formula, x ~ .), data = cbind(x, colData(dds))))[[1]][2]
)
sum_squares <- do.call(cbind, sum_squares)
x <- t(apply(sum_squares, 2, prop.table))
perVar <- x * mypca$percentVar
#colSums(perVar)
round(colSums(perVar) / sum(colSums(perVar)) * 100, 3)
```

```
# Site               Storage            Scion              Site:Storage       
#             27.377              1.163              6.201              3.608 
# Site:Scion         Storage:Scion      Site:Storage:Scion Residuals          
#             11.077              5.262              9.804             35.509
```

### PCA with canker lesion counts


```r
pca_summary <- apply(
  mypca$x[, 1:n_pcs], 2, 
  function(x){
    summary(aov(update(design_with_canker, x ~ .), data = as.data.frame(cbind(x, colData(dds)))))
  }
)

pca_summary
```

```
# $PC1
#                            Df Sum Sq Mean Sq F value  Pr(>F)   
# Site                        2 183473   91736 149.732 0.00663 **
# Storage                     1     29      29   0.047 0.84848   
# Scion                       6   2734     456   0.744 0.67074   
# Cankers                     1     83      83   0.136 0.74743   
# Site:Storage                2   2221    1111   1.813 0.35553   
# Site:Scion                 12   8158     680   1.110 0.56813   
# Storage:Scion               6   1710     285   0.465 0.80231   
# Site:Cankers                2    870     435   0.710 0.58480   
# Storage:Cankers             1    303     303   0.494 0.55477   
# Scion:Cankers               6   3085     514   0.839 0.63341   
# Site:Storage:Scion         11   6208     564   0.921 0.62869   
# Site:Storage:Cankers        2    115      58   0.094 0.91403   
# Site:Scion:Cankers         11   3731     339   0.554 0.79031   
# Storage:Scion:Cankers       5   1057     211   0.345 0.85403   
# Site:Storage:Scion:Cankers  3    468     156   0.254 0.85483   
# Residuals                   2   1225     613                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 8 observations deleted due to missingness
# 
# $PC2
#                            Df Sum Sq Mean Sq  F value   Pr(>F)    
# Site                        2 117349   58675 2292.220 0.000436 ***
# Storage                     1     63      63    2.466 0.256906    
# Scion                       6   1557     260   10.139 0.092486 .  
# Cankers                     1     16      16    0.616 0.514794    
# Site:Storage                2    427     213    8.335 0.107123    
# Site:Scion                 12   3860     322   12.566 0.076010 .  
# Storage:Scion               6    849     141    5.526 0.161130    
# Site:Cankers                2    292     146    5.699 0.149277    
# Storage:Cankers             1      1       1    0.040 0.859373    
# Scion:Cankers               6   2183     364   14.211 0.067193 .  
# Site:Storage:Scion         11   3230     294   11.473 0.082845 .  
# Site:Storage:Cankers        2    116      58    2.267 0.306088    
# Site:Scion:Cankers         11   3376     307   11.990 0.079446 .  
# Storage:Scion:Cankers       5   1521     304   11.885 0.079420 .  
# Site:Storage:Scion:Cankers  3    370     123    4.812 0.176843    
# Residuals                   2     51      26                      
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 8 observations deleted due to missingness
# 
# $PC3
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2   2729    1364   1.473  0.404
# Storage                     1     21      21   0.023  0.894
# Scion                       6  15566    2594   2.800  0.286
# Cankers                     1      0       0   0.000  0.987
# Site:Storage                2  12273    6136   6.624  0.131
# Site:Scion                 12  10972     914   0.987  0.608
# Storage:Scion               6  10110    1685   1.819  0.396
# Site:Cankers                2    322     161   0.174  0.852
# Storage:Cankers             1     40      40   0.043  0.855
# Scion:Cankers               6   2168     361   0.390  0.843
# Site:Storage:Scion         11  10272     934   1.008  0.598
# Site:Storage:Cankers        2   1233     617   0.666  0.600
# Site:Scion:Cankers         11   6836     621   0.671  0.733
# Storage:Scion:Cankers       5   2296     459   0.496  0.772
# Site:Storage:Scion:Cankers  3    414     138   0.149  0.922
# Residuals                   2   1853     926               
# 8 observations deleted due to missingness
# 
# $PC4
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2   5635  2817.7   4.432  0.184
# Storage                     1    448   447.6   0.704  0.490
# Scion                       6   2964   494.0   0.777  0.657
# Cankers                     1    273   272.7   0.429  0.580
# Site:Storage                2    439   219.3   0.345  0.744
# Site:Scion                 12  12850  1070.8   1.684  0.432
# Storage:Scion               6   2251   375.2   0.590  0.739
# Site:Cankers                2   2668  1333.8   2.098  0.323
# Storage:Cankers             1     11    11.4   0.018  0.906
# Scion:Cankers               6   2516   419.3   0.659  0.707
# Site:Storage:Scion         11   7600   690.9   1.087  0.573
# Site:Storage:Cankers        2    260   129.9   0.204  0.830
# Site:Scion:Cankers         11   2583   234.8   0.369  0.889
# Storage:Scion:Cankers       5    876   175.1   0.275  0.894
# Site:Storage:Scion:Cankers  3   1787   595.6   0.937  0.553
# Residuals                   2   1272   635.8               
# 8 observations deleted due to missingness
# 
# $PC5
#                            Df Sum Sq Mean Sq F value  Pr(>F)   
# Site                        2    303     152   33.68 0.02883 * 
# Storage                     1   3952    3952  877.50 0.00114 **
# Scion                       6    839     140   31.05 0.03153 * 
# Cankers                     1    552     552  122.63 0.00806 **
# Site:Storage                2   3907    1954  433.74 0.00230 **
# Site:Scion                 12   2729     227   50.49 0.01958 * 
# Storage:Scion               6   2906     484  107.54 0.00924 **
# Site:Cankers                2   1703     852  189.08 0.00526 **
# Storage:Cankers             1    381     381   84.69 0.01160 * 
# Scion:Cankers               6   1540     257   56.98 0.01735 * 
# Site:Storage:Scion         11   5162     469  104.20 0.00954 **
# Site:Storage:Cankers        2    306     153   33.94 0.02862 * 
# Site:Scion:Cankers         11   3860     351   77.91 0.01274 * 
# Storage:Scion:Cankers       5   1637     327   72.69 0.01363 * 
# Site:Storage:Scion:Cankers  3   3435    1145  254.20 0.00392 **
# Residuals                   2      9       5                   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 8 observations deleted due to missingness
# 
# $PC6
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2      8     3.9   0.013  0.987
# Storage                     1   1047  1047.0   3.523  0.201
# Scion                       6   4366   727.7   2.448  0.318
# Cankers                     1      3     3.2   0.011  0.927
# Site:Storage                2   1489   744.6   2.505  0.285
# Site:Scion                 12   2430   202.5   0.681  0.731
# Storage:Scion               6   2727   454.5   1.529  0.447
# Site:Cankers                2     14     7.1   0.024  0.977
# Storage:Cankers             1     68    68.1   0.229  0.679
# Scion:Cankers               6   1628   271.4   0.913  0.607
# Site:Storage:Scion         11   1737   157.9   0.531  0.802
# Site:Storage:Cankers        2   1292   646.1   2.174  0.315
# Site:Scion:Cankers         11   2653   241.2   0.811  0.671
# Storage:Scion:Cankers       5   1539   307.9   1.036  0.558
# Site:Storage:Scion:Cankers  3   1362   453.9   1.527  0.419
# Residuals                   2    594   297.2               
# 8 observations deleted due to missingness
# 
# $PC7
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2  562.5   281.2   0.777  0.563
# Storage                     1   34.6    34.6   0.096  0.786
# Scion                       6 1604.1   267.3   0.739  0.673
# Cankers                     1    0.2     0.2   0.000  0.986
# Site:Storage                2  795.5   397.8   1.100  0.476
# Site:Scion                 12 1494.7   124.6   0.344  0.906
# Storage:Scion               6  898.8   149.8   0.414  0.830
# Site:Cankers                2  313.5   156.7   0.433  0.698
# Storage:Cankers             1  295.4   295.4   0.817  0.462
# Scion:Cankers               6 2154.0   359.0   0.992  0.581
# Site:Storage:Scion         11 1858.0   168.9   0.467  0.836
# Site:Storage:Cankers        2 1131.1   565.6   1.563  0.390
# Site:Scion:Cankers         11 2681.2   243.7   0.674  0.731
# Storage:Scion:Cankers       5  420.6    84.1   0.233  0.918
# Site:Storage:Scion:Cankers  3  904.0   301.3   0.833  0.586
# Residuals                   2  723.6   361.8               
# 8 observations deleted due to missingness
# 
# $PC8
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2    283   141.6   0.630  0.613
# Storage                     1    182   182.4   0.812  0.463
# Scion                       6   1073   178.8   0.796  0.650
# Cankers                     1      3     2.5   0.011  0.925
# Site:Storage                2     92    45.9   0.204  0.830
# Site:Scion                 12   1138    94.9   0.422  0.864
# Storage:Scion               6   1475   245.9   1.095  0.549
# Site:Cankers                2   2198  1099.0   4.894  0.170
# Storage:Cankers             1     88    88.3   0.393  0.595
# Scion:Cankers               6   1277   212.8   0.947  0.595
# Site:Storage:Scion         11   3803   345.8   1.540  0.459
# Site:Storage:Cankers        2    102    50.8   0.226  0.815
# Site:Scion:Cankers         11   1894   172.2   0.767  0.690
# Storage:Scion:Cankers       5   1523   304.6   1.356  0.476
# Site:Storage:Scion:Cankers  3    902   300.5   1.338  0.455
# Residuals                   2    449   224.6               
# 8 observations deleted due to missingness
# 
# $PC9
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2    228   113.8   0.299  0.770
# Storage                     1    122   122.5   0.322  0.628
# Scion                       6    528    88.0   0.231  0.931
# Cankers                     1     49    49.3   0.130  0.753
# Site:Storage                2   3625  1812.3   4.764  0.173
# Site:Scion                 12   1940   161.7   0.425  0.863
# Storage:Scion               6   1399   233.2   0.613  0.728
# Site:Cankers                2    168    83.9   0.221  0.819
# Storage:Cankers             1     77    76.6   0.201  0.698
# Scion:Cankers               6    281    46.9   0.123  0.980
# Site:Storage:Scion         11   1970   179.1   0.471  0.834
# Site:Storage:Cankers        2    982   490.9   1.290  0.437
# Site:Scion:Cankers         11   2719   247.2   0.650  0.743
# Storage:Scion:Cankers       5     88    17.6   0.046  0.997
# Site:Storage:Scion:Cankers  3    660   220.0   0.578  0.683
# Residuals                   2    761   380.4               
# 8 observations deleted due to missingness
# 
# $PC10
#                            Df Sum Sq Mean Sq F value Pr(>F)
# Site                        2   62.8    31.4   0.236  0.809
# Storage                     1  290.4   290.4   2.185  0.277
# Scion                       6  732.1   122.0   0.918  0.605
# Cankers                     1  244.9   244.9   1.842  0.308
# Site:Storage                2 2115.5  1057.8   7.957  0.112
# Site:Scion                 12 2601.3   216.8   1.631  0.442
# Storage:Scion               6 1657.3   276.2   2.078  0.360
# Site:Cankers                2  128.0    64.0   0.482  0.675
# Storage:Cankers             1   13.9    13.9   0.104  0.777
# Scion:Cankers               6 1599.4   266.6   2.005  0.370
# Site:Storage:Scion         11 2198.2   199.8   1.503  0.466
# Site:Storage:Cankers        2  622.4   311.2   2.341  0.299
# Site:Scion:Cankers         11 1480.1   134.6   1.012  0.597
# Storage:Scion:Cankers       5  463.9    92.8   0.698  0.678
# Site:Storage:Scion:Cankers  3  694.3   231.4   1.741  0.385
# Residuals                   2  265.9   132.9               
# 8 observations deleted due to missingness
```

```r
# Extract PC scores as a list of dataframes
pcas <- lapply(pca_summary, function(i) data.frame(unclass(i)))

# Merge into single dataframe
pcs_factors_tidy <- lapply(
  names(pcas),
  function(name) {
    pcas[[name]] %>%
    mutate(
      PC = name, #substring(name, 3),
      Factor = gsub(" ", "", rownames(pcas[[name]])),
      var = Sum.Sq / sum(pcas[[name]]$Sum.Sq) * 100,
      pc_var = subset(pca_var, PC == name)$"perc_var",
      total_var = var * pc_var / 100,
      sig = case_when(
        is.na(Pr..F.) ~ "",
        Pr..F. < 0.001 ~ "***",
        Pr..F. < 0.01 ~ "**",
        Pr..F. < 0.05 ~ "*",
        TRUE ~ ""
      ),
      variance = ifelse(
        total_var < 0.01, paste0("<0.01", sig),
        paste0(round(total_var, 2), sig)
      )
    )
  }
) %>% bind_rows() %>% data.table()

# Order PCs and factors
pcs_factors_tidy$PC <- factor(pcs_factors_tidy$PC, levels = paste0("PC", 1:n_pcs))
pcs_factors_tidy$Factor <- factor(pcs_factors_tidy$Factor, levels = unique(pcs_factors_tidy$Factor))

# Significant factors
pcs_factors_tidy[
  Pr..F. < 0.05, 
  c("PC", "Factor", "Df", "F.value", "Pr..F.", "var", "pc_var", "total_var")
] %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:left;"> Factor </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> F.value </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> var </th>
   <th style="text-align:right;"> pc_var </th>
   <th style="text-align:right;"> total_var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 149.73167 </td>
   <td style="text-align:right;"> 0.0066343 </td>
   <td style="text-align:right;"> 85.1501398 </td>
   <td style="text-align:right;"> 18.6 </td>
   <td style="text-align:right;"> 15.8379260 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 2292.22027 </td>
   <td style="text-align:right;"> 0.0004361 </td>
   <td style="text-align:right;"> 86.7580635 </td>
   <td style="text-align:right;"> 12.7 </td>
   <td style="text-align:right;"> 11.0182741 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 33.68218 </td>
   <td style="text-align:right;"> 0.0288333 </td>
   <td style="text-align:right;"> 0.9132743 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.0273982 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Storage </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 877.50144 </td>
   <td style="text-align:right;"> 0.0011377 </td>
   <td style="text-align:right;"> 11.8964924 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.3568948 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 31.04713 </td>
   <td style="text-align:right;"> 0.0315297 </td>
   <td style="text-align:right;"> 2.5254793 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.0757644 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 122.62687 </td>
   <td style="text-align:right;"> 0.0080564 </td>
   <td style="text-align:right;"> 1.6624812 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.0498744 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site:Storage </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 433.73563 </td>
   <td style="text-align:right;"> 0.0023002 </td>
   <td style="text-align:right;"> 11.7605109 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.3528153 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site:Scion </td>
   <td style="text-align:right;"> 12 </td>
   <td style="text-align:right;"> 50.48597 </td>
   <td style="text-align:right;"> 0.0195806 </td>
   <td style="text-align:right;"> 8.2134005 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.2464020 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Storage:Scion </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 107.54378 </td>
   <td style="text-align:right;"> 0.0092412 </td>
   <td style="text-align:right;"> 8.7479775 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.2624393 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 189.08185 </td>
   <td style="text-align:right;"> 0.0052609 </td>
   <td style="text-align:right;"> 5.1268537 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.1538056 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Storage:Cankers </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 84.69040 </td>
   <td style="text-align:right;"> 0.0116026 </td>
   <td style="text-align:right;"> 1.1481676 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.0344450 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Scion:Cankers </td>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:right;"> 56.97860 </td>
   <td style="text-align:right;"> 0.0173471 </td>
   <td style="text-align:right;"> 4.6348334 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.1390450 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site:Storage:Scion </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 104.19848 </td>
   <td style="text-align:right;"> 0.0095429 </td>
   <td style="text-align:right;"> 15.5390755 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.4661723 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site:Storage:Cankers </td>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:right;"> 33.93941 </td>
   <td style="text-align:right;"> 0.0286210 </td>
   <td style="text-align:right;"> 0.9202490 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.0276075 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site:Scion:Cankers </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;"> 77.90771 </td>
   <td style="text-align:right;"> 0.0127389 </td>
   <td style="text-align:right;"> 11.6183444 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.3485503 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:right;"> 72.68466 </td>
   <td style="text-align:right;"> 0.0136266 </td>
   <td style="text-align:right;"> 4.9270149 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.1478104 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:left;"> Site:Storage:Scion:Cankers </td>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:right;"> 254.19961 </td>
   <td style="text-align:right;"> 0.0039211 </td>
   <td style="text-align:right;"> 10.3387310 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:right;"> 0.3101619 </td>
  </tr>
</tbody>
</table>

```r
# Table with factors as columns and PCs as rows
# pcs_factors <- dcast(pcs_factors_tidy, PC ~ Factor, value.var = "variance")
pcs_factors <- pcs_factors_tidy %>%
  select(PC, pc_var, Factor, variance) %>%
  spread(key = Factor, value = variance)

# Add a row with sum of each column by removing '<' and '*' and converting to numeric
pc_sums <- pcs_factors %>%
  mutate_all(~ gsub("[<|*]", "", .)) %>%
  mutate_all(as.numeric) %>%
  colSums()

pcs_factors %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:right;"> pc_var </th>
   <th style="text-align:left;"> Site </th>
   <th style="text-align:left;"> Storage </th>
   <th style="text-align:left;"> Scion </th>
   <th style="text-align:left;"> Cankers </th>
   <th style="text-align:left;"> Site:Storage </th>
   <th style="text-align:left;"> Site:Scion </th>
   <th style="text-align:left;"> Storage:Scion </th>
   <th style="text-align:left;"> Site:Cankers </th>
   <th style="text-align:left;"> Storage:Cankers </th>
   <th style="text-align:left;"> Scion:Cankers </th>
   <th style="text-align:left;"> Site:Storage:Scion </th>
   <th style="text-align:left;"> Site:Storage:Cankers </th>
   <th style="text-align:left;"> Site:Scion:Cankers </th>
   <th style="text-align:left;"> Storage:Scion:Cankers </th>
   <th style="text-align:left;"> Site:Storage:Scion:Cankers </th>
   <th style="text-align:left;"> Residuals </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:right;"> 18.6 </td>
   <td style="text-align:left;"> 15.84** </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.24 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> 0.7 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> 0.27 </td>
   <td style="text-align:left;"> 0.54 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.32 </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:right;"> 12.7 </td>
   <td style="text-align:left;"> 11.02*** </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.36 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> 0.3 </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> 0.32 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:right;"> 7.1 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 1.43 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 1.13 </td>
   <td style="text-align:left;"> 1.01 </td>
   <td style="text-align:left;"> 0.93 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> 0.95 </td>
   <td style="text-align:left;"> 0.11 </td>
   <td style="text-align:left;"> 0.63 </td>
   <td style="text-align:left;"> 0.21 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.17 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:right;"> 4.0 </td>
   <td style="text-align:left;"> 0.51 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.27 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 1.16 </td>
   <td style="text-align:left;"> 0.2 </td>
   <td style="text-align:left;"> 0.24 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> 0.68 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.11 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:right;"> 3.0 </td>
   <td style="text-align:left;"> 0.03* </td>
   <td style="text-align:left;"> 0.36** </td>
   <td style="text-align:left;"> 0.08* </td>
   <td style="text-align:left;"> 0.05** </td>
   <td style="text-align:left;"> 0.35** </td>
   <td style="text-align:left;"> 0.25* </td>
   <td style="text-align:left;"> 0.26** </td>
   <td style="text-align:left;"> 0.15** </td>
   <td style="text-align:left;"> 0.03* </td>
   <td style="text-align:left;"> 0.14* </td>
   <td style="text-align:left;"> 0.47** </td>
   <td style="text-align:left;"> 0.03* </td>
   <td style="text-align:left;"> 0.35* </td>
   <td style="text-align:left;"> 0.15* </td>
   <td style="text-align:left;"> 0.31** </td>
   <td style="text-align:left;"> &lt;0.01 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:right;"> 2.1 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.4 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> 0.25 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.24 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.05 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:right;"> 1.8 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> 0.17 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> 0.24 </td>
   <td style="text-align:left;"> 0.21 </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> 0.3 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0.08 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:right;"> 1.6 </td>
   <td style="text-align:left;"> 0.03 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.11 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.21 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.37 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> 0.15 </td>
   <td style="text-align:left;"> 0.09 </td>
   <td style="text-align:left;"> 0.04 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:right;"> 1.3 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.3 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.12 </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.16 </td>
   <td style="text-align:left;"> 0.08 </td>
   <td style="text-align:left;"> 0.23 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.06 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:right;"> 1.3 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.02 </td>
   <td style="text-align:left;"> 0.18 </td>
   <td style="text-align:left;"> 0.22 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.01 </td>
   <td style="text-align:left;"> &lt;0.01 </td>
   <td style="text-align:left;"> 0.14 </td>
   <td style="text-align:left;"> 0.19 </td>
   <td style="text-align:left;"> 0.05 </td>
   <td style="text-align:left;"> 0.13 </td>
   <td style="text-align:left;"> 0.04 </td>
   <td style="text-align:left;"> 0.06 </td>
   <td style="text-align:left;"> 0.02 </td>
  </tr>
</tbody>
</table>

### ADONIS


```r
# Calculate Bray-Curtis distance matrix
vg <- vegdist(t(counts(dds, normalize = T)), method = "bray")

formula <- update(FULL_DESIGN, vg ~ .)

set.seed(sum(utf8ToInt("Hamish McLean")))
result <- adonis2(formula, colData(dds), permutations = 1000)
result
```

```
# Permutation test for adonis under reduced model
# Permutation: free
# Number of permutations: 1000
# 
# adonis2(formula = formula, data = colData(dds), permutations = 1000)
#          Df SumOfSqs      R2      F   Pr(>F)    
# Model    41   9.9754 0.67431 2.0199 0.000999 ***
# Residual 40   4.8180 0.32569                    
# Total    81  14.7934 1.00000                    
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
df <- result %>% data.frame()
df$Perc.Var <- df$SumOfSqs / df["Total", "SumOfSqs"] * 100
df %>%
  kbl() %>%
  kable_styling("striped", full_width = F)
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> Df </th>
   <th style="text-align:right;"> SumOfSqs </th>
   <th style="text-align:right;"> R2 </th>
   <th style="text-align:right;"> F </th>
   <th style="text-align:right;"> Pr..F. </th>
   <th style="text-align:right;"> Perc.Var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Model </td>
   <td style="text-align:right;"> 41 </td>
   <td style="text-align:right;"> 9.975362 </td>
   <td style="text-align:right;"> 0.6743129 </td>
   <td style="text-align:right;"> 2.019934 </td>
   <td style="text-align:right;"> 0.000999 </td>
   <td style="text-align:right;"> 67.43129 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Residual </td>
   <td style="text-align:right;"> 40 </td>
   <td style="text-align:right;"> 4.818010 </td>
   <td style="text-align:right;"> 0.3256871 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 32.56871 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Total </td>
   <td style="text-align:right;"> 81 </td>
   <td style="text-align:right;"> 14.793372 </td>
   <td style="text-align:right;"> 1.0000000 </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> 100.00000 </td>
  </tr>
</tbody>
</table>

### NMDS ordination


```r
set.seed(sum(utf8ToInt("Hamish McLean")))
ord <- metaMDS(vg,trace=0) 
#sratmax=20000,maxit=20000,try = 177, trymax = 177

bac_nmds <- scores(ord)

bac_nmds_plot <- plotOrd(
  bac_nmds, colData(dds), 
  design = "Site", 
  shape = "Storage", 
  alpha = 0.75, cbPalette = T
) #+ theme(text = element_text(size = 14))

ggsave(filename = "fun_nmds_plot.png", plot = bac_nmds_plot, path = "figures/")

bac_nmds_plot
```

![](root_endophytes_files/figure-html/BAC NMDS-1.png)<!-- -->

## ASV abundance

### Explore distribution of ASV counts


```r
# Extract normalised counts from DESeq object
asv_counts <- counts(dds, normalize = T) %>% as.data.frame()

# Sum ASV counts across samples
total_asv_counts <- rowSums(asv_counts)

# Sort ASVs by abundance
total_asv_counts <- total_asv_counts[order(total_asv_counts, decreasing = T)]

# Caculate cumulative percentage
cumulative <- data.frame(
  cumulative = cumsum(total_asv_counts) / sum(total_asv_counts) * 100,
  no = seq_along(total_asv_counts)
)

# Plot cumulative percentage of ASVs
bac_cum_asv <- ggline(
  data = cumulative, x = "no", y = "cumulative", 
  plot_type = "l", palette = cbPalette,
  title = "Cumulative percentage of bacterial ASVs", xlab = "Number of ASVs", 
  ylab = "Cumulative percentage of reads"
)
ggsave(filename = "bac_cum_asv.png", plot = bac_cum_asv, path = "figures/")
bac_cum_asv
```

![](root_endophytes_files/figure-html/BAC top ASVs-1.png)<!-- -->

```r
# Find the number of ASVs that account for 50%, 80%, and 99% of total reads
cat(
  "Number of ASVs that account for 50%, 80%, 90%, and 99% of total reads", "\n\n",
  "50%:", sum(cumulative <= 50), "\n",
  "80%:", sum(cumulative <= 80), "\n",
  "90%:", sum(cumulative <= 90), "\n",
  "99%:", sum(cumulative <= 99), "\n"
)
```

```
# Number of ASVs that account for 50%, 80%, 90%, and 99% of total reads 
# 
#  50%: 205 
#  80%: 1055 
#  90%: 2036 
#  99%: 5037
```

```r
# Find the cumulative percentage accounted for by top x ASVs
cat(
  "Percentage of total reads accounted for by the top 100, 200,and 500 ASVs:", "\n\n",
  "100:", round(cumulative[cumulative$no == 100, "cumulative"], 1) , "\n",
  "200:", round(cumulative[cumulative$no == 200, "cumulative"], 1) , "\n",
  "500:", round(cumulative[cumulative$no == 500, "cumulative"], 1) , "\n"
)
```

```
# Percentage of total reads accounted for by the top 100, 200,and 500 ASVs: 
# 
#  100: 43.7 
#  200: 53.9 
#  500: 69.1
```

```r
# Average ASV counts in order
mean_asv_counts <- rowMeans(asv_counts)
mean_asv_counts <- mean_asv_counts[order(mean_asv_counts, decreasing = T)]

# Plot read count distribution
bac_asv_counts <- ggline(
  data = data.frame(ASV = seq_along(mean_asv_counts), counts = mean_asv_counts),
  x = "ASV", y = "counts", plot_type = "l",
  title = "Bacterial ASV read count distribution", xlab = "ASV", ylab = "Mean read count"
)
ggsave(filename = "bac_asv_counts.png", plot = bac_asv_counts, path = "figures/")
bac_asv_counts
```

![](root_endophytes_files/figure-html/BAC top ASVs-2.png)<!-- -->

```r
# Number of ASVs with mean read count > 100, 200, and 500
cat(
  "Number of ASVs with mean read count > 100, 200, and 500", "\n\n",
  "100:", sum(rowMeans(asv_counts) > 100), "\n",
  "200:", sum(rowMeans(asv_counts) > 200), "\n",
  "500:", sum(rowMeans(asv_counts) > 500), "\n"
)
```

```
# Number of ASVs with mean read count > 100, 200, and 500 
# 
#  100: 71 
#  200: 32 
#  500: 8
```

### Filter top ASVs with 100 % of reads


```r
# Filter the top x abundant ASVs by the sum of their normalised counts
# top_asvs <- asv_counts[order(rowSums(asv_counts), decreasing = T)[1:DIFFOTU], ]

# Filter ASVs with mean read count > 100
# top_asvs <- asv_counts[rowMeans(asv_counts) > 100, ]

# Filter top ASVs with 90% of reads
# top_asvs <- asv_counts[filter_otus(asv_counts, BACASVFILTER), ]

# No filtering
top_asvs <- asv_counts

# Check that sample names match
identical(names(top_asvs), rownames(colData))
```

```
# [1] TRUE
```

```r
# Extract taxonomic data for top ASVs
top_taxa <- taxData[rownames(top_asvs), ]

# Log transform normalised counts
# top_asvs <- log10(top_asvs + 1) # Log transform in models instead

top_asv_data <- data.frame(t(top_asvs))
top_asv_ids <- rownames(top_asvs)
identical(rownames(top_asv_data), rownames(colData))
```

```
# [1] TRUE
```

```r
top_asv_data <- merge(top_asv_data, colData, by = 0) %>% column_to_rownames("Row.names")
```

### Effect of design factors on abundance of top ASVs

Effect of Site, Scion, and Storage on abundance of top 200 ASVs


```r
# Perform ANOVA on list of top ASVs
top_asvs_anova_results <- lapply(
  top_asv_ids, 
  function(asv) {
    asv_lm_anova(asv, formula, top_asv_data) %>%
    extend_asv_anova(asv)
  }
) %>% bind_rows() %>% data.table()

# Group by factor and adjust p-values
top_asvs_anova_results <- top_asvs_anova_results %>% 
  group_by(Factor) %>% 
  mutate(p.adj = p.adjust(`Pr..F.`, method = "BH")) %>% 
  data.table()

# Order factors by original order
top_asvs_anova_results$Factor <- factor(top_asvs_anova_results$Factor, levels = unique(top_asvs_anova_results$Factor))

# Summary of top ASV ANOVA results
top_asvs_anova_summary <- top_asvs_anova_results %>% 
  select(ASV, Taxonomy, Abundance, Factor, var, p.adj) %>% 
  pivot_wider(names_from = Factor, values_from = c(var, p.adj), names_glue = "{Factor}_{.value}") %>%
  select(
    ASV, Taxonomy, Abundance, Site_var, Site_p.adj, Storage_var, Storage_p.adj, 
    Scion_var, Scion_p.adj, "Site:Storage_var", "Site:Storage_p.adj", "Site:Scion_var", 
    "Site:Scion_p.adj", "Storage:Scion_var", "Storage:Scion_p.adj", "Site:Storage:Scion_var", 
    "Site:Storage:Scion_p.adj"
  ) %>%
  data.table()

# top_asvs_anova_summary %>%
#   kbl() %>%
#   kable_styling("striped", full_width = F)

cat(
  "Number of ASVs with statistically significant (*P* < 0.05) adjusted p-values", "\n\n",
  "Site:", nrow(top_asvs_anova_results[Factor == "Site" & p.adj < 0.05, ]), "\n",
  "Storage:", nrow(top_asvs_anova_results[Factor == "Storage" & p.adj < 0.05, ]), "\n",
  "Scion:", nrow(top_asvs_anova_results[Factor == "Scion" & p.adj < 0.05, ]), "\n",
  "Site:Storage:", nrow(top_asvs_anova_results[Factor == "Site:Storage" & p.adj < 0.05, ]), "\n",
  "Site:Scion:", nrow(top_asvs_anova_results[Factor == "Site:Scion" & p.adj < 0.05, ]), "\n",
  "Storage:Scion:", nrow(top_asvs_anova_results[Factor == "Storage:Scion" & p.adj < 0.05, ]), "\n",
  "Site:Storage:Scion:", nrow(top_asvs_anova_results[Factor == "Site:Storage:Scion" & p.adj < 0.05, ]), "\n\n",
  "Total ASVs:", length(unique(top_asvs_anova_results$ASV)), "\n\n"
)

# Filter by significant effect of scion and its interactions
scion_asvs <- top_asvs_anova_results[grepl("Scion", Factor) & p.adj < 0.05, ]

# scion_asvs %>%
#   kbl() %>%
#   kable_styling("striped", full_width = F)

cat(
  length(scion_asvs$ASV), 
  "ASVs with significant (*P* < 0.05) adjusted p-values for the effect of Scion and its interactions.", "\n\n"
)

# Summary of ASVs with significant Scion effect
# bac_significant_scion <- top_asvs_anova_summary[ASV %in% scion_asvs$ASV, ] %>%
#   arrange(desc(Abundance))

# kbl(bac_significant_scion) %>%
#   kable_styling("striped")

# kbl(bac_significant_scion) %>%
#   kable_styling("striped") %>%
#   save_kable("tables/BAC_scion_asvs.html")

top_asvs_anova_summary %>%
  arrange(desc(Abundance)) %>%
  kbl() %>%
  kable_styling("striped") %>%
  save_kable("tables/BAC_all_asvs.html")
```


```r
# Export significant ASVs as fasta

# Read BAC ASVs
BAC_asvs <- read.fasta("data/BAC.zotus.fa")
# Replace 'OTU' with 'ASV' in sequence names
names(BAC_asvs) <- gsub("OTU", "ASV", names(BAC_asvs))

# Write significant ASVs to fasta
write.fasta(
  sequences = BAC_asvs[scion_asvs$ASV],
  names = paste(scion_asvs$ASV, taxData[scion_asvs$ASV, "rank"]),
  file = "fasta/BAC_scion_asvs.fasta"
)
```

## Canker counts

Testing the effects of of ASV abundance, α-diversity, and β-diversity on canker counts.

This uses a nested negative binomial regression model.

The base model for canker counts uses the formula: Cankers ~ Site * Storage * Scion.


```r
# Filter out samples with missing canker count
canker_abundance_data <- colData[complete.cases(colData$Cankers), ]

# Base model
canker_design = "Cankers ~ Site * Storage * Scion"
base_model <- glm.nb(canker_design, data = canker_abundance_data)

# Abundance model
abundance_design = paste(canker_design, "+ log(copy_number)")
abundance_model <- glm.nb(abundance_design, data = canker_abundance_data)

# ANOVA of abundance with canker count
anova(base_model, abundance_model) %>%
  kbl() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> Model </th>
   <th style="text-align:right;"> theta </th>
   <th style="text-align:right;"> Resid. df </th>
   <th style="text-align:right;">    2 x log-lik. </th>
   <th style="text-align:left;"> Test </th>
   <th style="text-align:right;">    df </th>
   <th style="text-align:right;"> LR stat. </th>
   <th style="text-align:right;"> Pr(Chi) </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Site * Storage * Scion </td>
   <td style="text-align:right;"> 2.882191 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> -554.6118 </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
   <td style="text-align:right;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Site * Storage * Scion + log(copy_number) </td>
   <td style="text-align:right;"> 2.980708 </td>
   <td style="text-align:right;"> 32 </td>
   <td style="text-align:right;"> -551.9342 </td>
   <td style="text-align:left;"> 1 vs 2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 2.677645 </td>
   <td style="text-align:right;"> 0.1017661 </td>
  </tr>
</tbody>
</table>

### Effect of ASV abundance on canker count


```r
# Filter out samples with missing canker count
top_asv_data <- top_asv_data[complete.cases(top_asv_data$Cankers), ]

# Base model design
canker_design = "Cankers ~ Site * Storage * Scion"

# Base model with ASV abundance data
base_model <- glm.nb(canker_design, data = top_asv_data)

# Effect of ASV abundance on canker count for top ASVs
asv_canker_results <- lapply(
  top_asv_ids, 
  function(x) asv_canker_anova(x, canker_design, base_model, top_asv_data)
) %>% bind_rows() %>% data.table()

# Adjust p-values for multiple testing
asv_canker_results$p_adjusted <- p.adjust(asv_canker_results$p, method = "BH")

# Summary of ASVs with statistically significant (*P* < 0.05) adjusted p-values
cat(
  nrow(asv_canker_results[p_adjusted < 0.05, ]), 
  "ASVs have statistically significant (*P* < 0.05) adjusted p-values"
)

bac_sig_canker_asvs <- asv_canker_results[p_adjusted < 0.05, ] %>%
  arrange(desc(Abundance))

# kbl(bac_sig_canker_asvs) %>%
#   kable_styling("striped")

kbl(bac_sig_canker_asvs) %>%
  kable_styling("striped") %>%
  save_kable("tables/BAC_canker_asvs.html")

cat(
  "\nSignificant canker bacterial ASVs\n",
  "\nTotal ASVs: ", nrow(asv_canker_results),
  "\nWarnings: ", nrow(asv_canker_results[!is.na(warning)]),
  "\nTotal ASVs without warnings: ", nrow(asv_canker_results[is.na(warning)]),
  "\nTotal significant ASVs: ", nrow(asv_canker_results[p_adjusted < 0.05]),
  "\nASVs with negative effect: ", nrow(asv_canker_results[p_adjusted < 0.05 & coef < 0]),
  "\nASVs with positive effect: ", nrow(asv_canker_results[p_adjusted < 0.05 & coef > 0]), "\n"
)
```

#### Effect of ASV abundance on canker count per site

Filter top ASVs with 100 % of reads per site
and test the effect of ASV abundance on canker count per site.


```r
# For each site, select top ASVs with 50% of reads
top_asvs_per_site <- lapply(
  unique(colData$Site),
  function(site) {
    samples <- filter(colData, Site == site)
    top_asv_data <- select(asv_counts, rownames(samples))
    # top_asvs <- filter(top_asv_data, rowMeans(top_asv_data) > 100)
    # top_asvs <- top_asv_data[filter_otus(top_asv_data, BACASVFILTER), ]
    top_asvs <- top_asv_data[rowSums(top_asv_data) > 0, ]
    # top_asvs <- top_asvs[ , colSums(is.na(top_asvs)) == 0] # Remove ASVs with missing values
    top_asv_ids <- rownames(top_asvs)
    top_asvs <- data.frame(t(top_asvs)) %>% merge(samples, by = 0) %>% column_to_rownames("Row.names")
    top_asvs <- top_asvs[complete.cases(top_asvs$Cankers), ]
    return(list(asvs = top_asv_ids, data = top_asvs))
  }
)

# Summary of top ASVs per site
data.table(
  Site = 1:3,
  topASVs = sapply(top_asvs_per_site, function(x) length(x$asvs)),
  totalASVs = sapply(top_asvs_per_site, function(x) sum(rowSums(asv_counts[, rownames(x$data)]) > 0, na.rm = TRUE)),
  topreads = sapply(top_asvs_per_site, function(x) sum(select(x$data, starts_with("ASV")))),
  totalreads = sapply(top_asvs_per_site, function(x) sum(asv_counts[, rownames(x$data)]))
) %>% mutate(
  "ASV%" = round(topASVs / totalASVs * 100),
  "Read%" = round(topreads / totalreads * 100)
) %>% kbl() %>% kable_styling("striped")

canker_site_design <- "Cankers ~ Storage * Scion"

# ANOVA of ASV abundance with canker count per ASV
asv_canker_site_anova <- function(asvs, data) {
  base_model <- glm.nb(canker_site_design, data = data)
  results <- lapply(
    asvs, 
    function(asv) asv_canker_anova(asv, canker_site_design, base_model, data)
  ) %>% bind_rows() %>% data.table()
  results$p_adjusted <- p.adjust(results$p, method = "BH")
  return(results)
}

# Run ANOVA per site
asv_canker_site_results <- lapply(
  top_asvs_per_site,
  function(x) asv_canker_site_anova(x$asvs, x$data)
)

# Add site to each result as new column and merge into single data.table
asv_canker_site_results <- lapply(
  1:3, 
  function(site) {
    result <- asv_canker_site_results[[site]]
    result$Site <- site
    return(result)
  }
) %>% bind_rows()

# Significant ASVs
significant_asvs <- asv_canker_site_results[p_adjusted < 0.05 & is.na(warning), ]

bac_canker_site_asvs <- significant_asvs[, c("Site", "ASV", "Taxonomy", "Abundance", "coef", "var", "p_adjusted", "warning")] %>%
  arrange(Site, desc(Abundance))

# kbl(bac_canker_site_asvs) %>%
#   kable_styling("striped")
  
kbl(bac_canker_site_asvs) %>%
  kable_styling("striped") %>%
  save_kable("tables/BAC_significant_site_canker_asvs.html")

for(site in 1:3){
  cat(
    "\nSignificant canker fungal ASVs for site ", site, "\n",
    "\nTotal ASVs: ", nrow(asv_canker_site_results[Site == site]),
    "\nWarnings: ", nrow(asv_canker_site_results[Site == site & !is.na(warning)]),
    "\nTotal ASVs without warnings: ", nrow(asv_canker_site_results[Site == site & is.na(warning)]),
    "\nTotal significant ASVs: ", nrow(asv_canker_site_results[Site == site & p_adjusted < 0.05]),
    "\nASVs with negative effect: ", nrow(asv_canker_site_results[Site == site & p_adjusted < 0.05 & coef < 0]),
    "\nASVs with positive effect: ", nrow(asv_canker_site_results[Site == site & p_adjusted < 0.05 & coef > 0]), "\n"
  )
}
```


```r
# Export significant ASVs as FASTA
write.fasta(
  sequences = BAC_asvs[as.character(significant_asvs$ASV)],
  names = paste(significant_asvs$ASV, taxData[as.character(significant_asvs$ASV), "rank"]),
  file = "fasta/BAC_canker_asvs.fasta"
)
```

##### Plot of ASV abundance against canker count


```r
# List of significant ASVs
significant_asv_list <- significant_asvs$ASV %>% unlist()

significant_asv_data <- asv_counts[significant_asv_list, ] %>% 
  t() %>% 
  data.frame() %>% 
  merge(colData, by = 0) %>% 
  column_to_rownames("Row.names") %>%
  select(c(significant_asv_list, "Site", "Storage", "Scion", "Cankers"))

# Melt data for ggplot
significant_asv_long_data <- significant_asv_data %>% reshape2::melt(
  id.vars = c("Site", "Storage", "Scion", "Cankers"), variable.name = "ASV", value.name = "Abundance"
)

# Log trasnform abundance
significant_asv_long_data$log10_abundance <- log10(significant_asv_long_data$Abundance + 1)

bac_asv_canker_plot <- ggscatter(
  data = significant_asv_long_data, x = "log10_abundance", y = "Cankers", 
  color = "Storage", facet.by = c("ASV", "Site"),
  xlab = "ASV abundance (log10)", ylab = "Canker count",
  palette = cbPalette, legend = "bottom"
)

ggsave(
  filename = "bac_asv_canker_plot.png", plot = bac_asv_canker_plot, path = "figures/",
  height = 40, width = 20, units = "cm"
)

bac_asv_canker_plot
```

#### Effect of aggregated genera on canker counts


```r
# Add genus from taxData to countData
bac_genus_data <- counts(dds, normalize = T) %>% as.data.frame() %>% mutate(
  genus = taxData[rownames(countData), "genus"]
)

# Group by genus and set genus to rownames
bac_genus_data <- bac_genus_data %>% group_by(genus) %>% summarise_all(sum) %>% as.data.frame()

# Set rownames as genus
rownames(bac_genus_data) <- bac_genus_data$genus
bac_genus_data <- dplyr::select(bac_genus_data, -genus)

# Filter genera with mean abundance < 100
bac_genus_data <- bac_genus_data[rowMeans(bac_genus_data) > 10, ]

# Rank not genus
not_genus <- rownames(bac_genus_data)[grep("\\([a-z]\\)", rownames(bac_genus_data))]
# Remove rows with genus in not_genus
bac_genus_data <- bac_genus_data[!rownames(bac_genus_data) %in% not_genus, ]
cat(
  length(not_genus), "non-genus ranks removed:\n\n",
  not_genus, "\n"
)

# Remove any with slashes or underscores
removals <- rownames(bac_genus_data)[grep("[/_]", rownames(bac_genus_data))]
# Remove rows with symbols in the name
bac_genus_data <- bac_genus_data[!rownames(bac_genus_data) %in% removals, ]
cat(
  length(removals), "ranks with slashes or underscores removed:\n\n",
  removals, "\n"
)

# Set rownames as genus
bac_genera <- rownames(bac_genus_data)

# Transpose and add metadata from colData
bac_genus_data <- t(bac_genus_data) %>% as.data.frame() %>% mutate(
  Site = colData[rownames(.), "Site"],
  Storage = colData[rownames(.), "Storage"],
  Scion = colData[rownames(.), "Scion"],
  Cankers = colData[rownames(.), "Cankers"]
)

# Filter out samples with missing canker count
bac_genus_data <- bac_genus_data[complete.cases(bac_genus_data$Cankers), ]
```


```r
# Base model design
canker_design = "Cankers ~ Site * Storage * Scion"

# Base model with genus abundance data
base_model <- glm.nb(canker_design, data = bac_genus_data)

# ANOVA of genus abundance with canker count
genus_canker_anova <- function(genus, design, base_model, data) {
  log_genus = paste0("log(", genus, " + 1)")
  f = paste(design, "+", log_genus)#, "+", log_genus, ":Site")
  m = glm.nb(f, data = data)
  a = anova(base_model, m) %>% data.frame()
  b = suppressWarnings(anova(m)) %>% data.frame()
  total_deviance = sum(b$Deviance, na.rm = T) + tail(b$Resid..Dev, 1)
  d = data.table(
    Genus = genus,
    Abundance = mean(data[[genus]]),
    coef = m$coefficients[log_genus],
    var = b[log_genus, 'Deviance'] / total_deviance * 100,
    p = a[2, 'Pr.Chi.']
  )
  return(d)
}

# genus_canker_anova(bac_genera[1], canker_design, base_model, bac_genus_data)

# Effect of genera abundance on canker counts
genus_canker_results <- lapply(
  bac_genera, 
  function(x) genus_canker_anova(x, canker_design, base_model, bac_genus_data)
) %>% bind_rows() %>% data.table()

# Adjust p-values for multiple testing
genus_canker_results$p_adjusted <- p.adjust(genus_canker_results$p, method = "BH")

# Summary of ASVs with statistically significant (*P* < 0.05) adjusted p-values
cat(
  nrow(genus_canker_results[p_adjusted < 0.05, ]), "of", nrow(genus_canker_results),
  "genera have statistically significant (*P* < 0.05) adjusted p-values\n\n"
)
if(nrow(genus_canker_results[p_adjusted < 0.05, ]) > 0) {
  genus_canker_results[p_adjusted < 0.05, ] %>%
    arrange(desc(Abundance)) %>%
    kbl() %>%
    kable_styling("striped")
}

sig_genus <- genus_canker_results[p_adjusted < 0.05]$Genus %>% unlist()

for (genus in sig_genus) {
  log_genus = paste0("log(", genus, " + 1)")
  f = paste(canker_design, "+", log_genus)
  m = glm.nb(f, data = bac_genus_data)
  print(anova(base_model, m))
}
```

##### Plot of genus abundance against canker count


```r
significant_genera <- genus_canker_results[p_adjusted < 0.05]$Genus %>% unlist()

significant_genera_data <- bac_genus_data[, c(significant_genera, FACTORS, "Cankers")]

# Melt data for ggplot
significant_genera_data <- significant_genera_data %>% reshape2::melt(
  id.vars = c(FACTORS, "Cankers"), variable.name = "Genus", value.name = "Abundance"
)

# Log transform abundance
significant_genera_data$log10_abundance <- log10(significant_genera_data$Abundance + 1)

# Plot of genus abundance against canker count
bac_genus_canker_plot <- ggscatter(
  data = significant_genera_data, x = "log10_abundance", y = "Cankers", 
  color = "Site", shape = "Storage",
  xlab = "Genus abundance (log10)", ylab = "Canker count",
  free_x = T, free_y = T, palette = cbPalette
) + facet_wrap(~ Genus, scales = "free")

ggsave(
  filename = "bac_genus_canker_plot.png", plot = bac_genus_canker_plot, path = "figures/",
  height = 10, width = 20, units = "cm"
)
```

### Effect of α-diversity on canker count


```r
# ANOVA of α-diversity with canker count

# Base model with α-diversity data
base_model <- glm.nb(canker_design, data = all_alpha_ord)

measures <- c("S.chao1", "shannon", "simpson")

# ANOVA of α-diversity with canker count
alpha_canker_anova <- function(measure, data) {
  f = paste(canker_design, "+", measure)
  m = glm.nb(f, data = data)
  a = anova(base_model, m) %>% data.frame()
  b = anova(m) %>% data.frame()
  total_deviance = sum(b$Deviance, na.rm = T) + tail(b$Resid..Dev, 1)
  d = data.frame(
    measure = measure,
    df = a[2, 'X...df'],
    coef = m$coefficients[measure],
    lr = a[2, 'LR.stat.'],
    p = a[2, 'Pr.Chi.'],
    var = b[measure, 'Deviance'] / total_deviance * 100
  )
  return(d)
}

# Effect of α-diversity on canker count for each measure
alpha_canker_results <- data.table(bind_rows(lapply(measures, function(x) alpha_canker_anova(x, all_alpha_ord))))

alpha_canker_results %>%
  kbl() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> measure </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> coef </th>
   <th style="text-align:right;"> lr </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:right;"> var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> S.chao1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.000326 </td>
   <td style="text-align:right;"> 0.8269849 </td>
   <td style="text-align:right;"> 0.3631458 </td>
   <td style="text-align:right;"> 0.1949192 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> shannon </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 1.047427 </td>
   <td style="text-align:right;"> 4.4866793 </td>
   <td style="text-align:right;"> 0.0341600 </td>
   <td style="text-align:right;"> 3.0434692 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> simpson </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 74.201152 </td>
   <td style="text-align:right;"> 6.3333138 </td>
   <td style="text-align:right;"> 0.0118491 </td>
   <td style="text-align:right;"> 2.5540565 </td>
  </tr>
</tbody>
</table>

```r
# ANOVA results
for (measure in measures) {
  f = paste(canker_design, "+", measure)
  m = glm.nb(f, data = all_alpha_ord)
  print(anova(base_model, m))
}
```

```
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: Cankers
#                              Model    theta Resid. df    2 x log-lik.   Test
# 1           Site * Storage * Scion 2.882191        33       -554.6118       
# 2 Site * Storage * Scion + S.chao1 2.904996        32       -553.7848 1 vs 2
#      df  LR stat.   Pr(Chi)
# 1                          
# 2     1 0.8269849 0.3631458
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: Cankers
#                              Model    theta Resid. df    2 x log-lik.   Test
# 1           Site * Storage * Scion 2.882191        33       -554.6118       
# 2 Site * Storage * Scion + shannon 3.049289        32       -550.1251 1 vs 2
#      df LR stat.    Pr(Chi)
# 1                          
# 2     1 4.486679 0.03415997
# Likelihood ratio tests of Negative Binomial Models
# 
# Response: Cankers
#                              Model    theta Resid. df    2 x log-lik.   Test
# 1           Site * Storage * Scion 2.882191        33       -554.6118       
# 2 Site * Storage * Scion + simpson 3.136392        32       -548.2785 1 vs 2
#      df LR stat.    Pr(Chi)
# 1                          
# 2     1 6.333314 0.01184907
```

Plot of Simpson against canker counts


```r
ggscatter(all_alpha_ord, y = "Cankers", x = "simpson", add = "reg.line", conf.int = T, cor.coef = T)
```

![](root_endophytes_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
ggsave("figures/bac_alpha_test.png")

m <- glm.nb("Cankers ~ Site * Storage * Scion + simpson", data = all_alpha_ord)
```

Spearman correlation of α-diversity indices with canker lesion



```r
cor.test(all_alpha_ord$S.chao1, all_alpha_ord$Cankers, method = "spearman", use = "complete.obs")
```

```
# 
# 	Spearman's rank correlation rho
# 
# data:  all_alpha_ord$S.chao1 and all_alpha_ord$Cankers
# S = 45495, p-value = 0.004557
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.3262464
```

```r
cor.test(all_alpha_ord$shannon, all_alpha_ord$Cankers, method = "spearman", use = "complete.obs")
```

```
# 
# 	Spearman's rank correlation rho
# 
# data:  all_alpha_ord$shannon and all_alpha_ord$Cankers
# S = 33051, p-value = 3.366e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5105344
```

```r
cor.test(all_alpha_ord$simpson, all_alpha_ord$Cankers, method = "spearman", use = "complete.obs")
```

```
# 
# 	Spearman's rank correlation rho
# 
# data:  all_alpha_ord$simpson and all_alpha_ord$Cankers
# S = 32444, p-value = 2.104e-06
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho 
# 0.5195309
```

### Effect of β-diversity on canker count


```r
no_pcs <- 10

# Merge PC scores with canker data
pc_scores <- merge(colData, data.frame(mypca$x[, 1:no_pcs]), by = "row.names") %>% 
  column_to_rownames("Row.names")

pcs <- tail(colnames(pc_scores), no_pcs)

# Base model with β-diversity data
base_model <- glm.nb(canker_design, data = pc_scores)

# ANOVA of β-diversity with canker count
beta_canker_anova <- function(pc, data) {
  f = paste0(canker_design, "+", pc)
  m = glm.nb(f, data = data)
  a = anova(base_model, m) %>% data.frame()
  b = anova(m) %>% data.frame()
  total_deviance = sum(b$Deviance, na.rm = T) + tail(b$Resid..Dev, 1)
  d = data.frame(
    PC = pc,
    df = a[2, 'X...df'],
    coef = m$coefficients[pc],
    lr = a[2, 'LR.stat.'],
    p = a[2, 'Pr.Chi.'],
    var = b[pc, 'Deviance'] / total_deviance * 100
  )
  return(d)
}

# Effect of β-diversity on canker count for each PC
beta_canker_results <- data.table(bind_rows(lapply(pcs, function(x) beta_canker_anova(x, pc_scores))))

beta_canker_results %>%
  kbl() %>%
  kable_styling("striped")
```

<table class="table table-striped" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> PC </th>
   <th style="text-align:right;"> df </th>
   <th style="text-align:right;"> coef </th>
   <th style="text-align:right;"> lr </th>
   <th style="text-align:right;"> p </th>
   <th style="text-align:right;"> var </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> PC1 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0247556 </td>
   <td style="text-align:right;"> 14.5763203 </td>
   <td style="text-align:right;"> 0.0001346 </td>
   <td style="text-align:right;"> 0.2245606 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC2 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0071783 </td>
   <td style="text-align:right;"> 0.9833974 </td>
   <td style="text-align:right;"> 0.3213615 </td>
   <td style="text-align:right;"> 0.0276326 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC3 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0036667 </td>
   <td style="text-align:right;"> 0.4855541 </td>
   <td style="text-align:right;"> 0.4859172 </td>
   <td style="text-align:right;"> 0.3145477 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC4 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0157086 </td>
   <td style="text-align:right;"> 6.8155264 </td>
   <td style="text-align:right;"> 0.0090369 </td>
   <td style="text-align:right;"> 0.2983212 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC5 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0040380 </td>
   <td style="text-align:right;"> 0.4772993 </td>
   <td style="text-align:right;"> 0.4896482 </td>
   <td style="text-align:right;"> 3.4298480 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC6 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0092406 </td>
   <td style="text-align:right;"> 1.0251900 </td>
   <td style="text-align:right;"> 0.3112911 </td>
   <td style="text-align:right;"> 0.9320288 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC7 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> -0.0017239 </td>
   <td style="text-align:right;"> 0.0318831 </td>
   <td style="text-align:right;"> 0.8582844 </td>
   <td style="text-align:right;"> 0.0292505 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC8 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0127517 </td>
   <td style="text-align:right;"> 2.2531119 </td>
   <td style="text-align:right;"> 0.1333460 </td>
   <td style="text-align:right;"> 0.4595563 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC9 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0061285 </td>
   <td style="text-align:right;"> 0.3383622 </td>
   <td style="text-align:right;"> 0.5607761 </td>
   <td style="text-align:right;"> 0.0412642 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> PC10 </td>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:right;"> 0.0012980 </td>
   <td style="text-align:right;"> 0.0097534 </td>
   <td style="text-align:right;"> 0.9213293 </td>
   <td style="text-align:right;"> 3.5283651 </td>
  </tr>
</tbody>
</table>

# Extra figures

## Abundance


```r
abundance_combined <- rbind(
  as.data.frame(colData(ubiome_FUN$dds)), as.data.frame(colData(ubiome_BAC$dds))
) %>% mutate(kingdom = ifelse(Target == "ITS", "Fungi", "Bacteria")) %>%
  mutate(kingdom = factor(kingdom, levels = c("Fungi", "Bacteria"))) %>%
  mutate(Storage = factor(Storage, levels = c("no", "yes")))

# abundance_bar <- ggbarplot(
#   data = abundance_combined, x = "Storage", y = "log_copy_number", 
#   fill = "Site", add = "mean_se", facet.by = "kingdom",
#   palette = cbPalette, position = position_dodge(0.8),
#   ylab = "Mean copy number (log10)", xlab = F, legend = "right"
# ) + guides(fill = guide_legend(title = "Site"))

# ggsave(
#   filename = "abundance_bar.png", plot = abundance_bar, path = "figures/", 
#   height = 12, width = 24, units = "cm"
# )

# abundance_bar

abundance_box <- ggboxplot(
  data = abundance_combined, x = "Site", y = "log_copy_number", 
  color = "Storage", add = "jitter", facet.by = "kingdom",
  palette = cbPalette, legend = "bottom",
  ylab = "Copy number (log10)", xlab = "Site"
) 

ggsave(
  filename = "abundance_box.png", plot = abundance_box, path = "figures/", 
  height = 12, width = 24, units = "cm"
)

abundance_box
```

![](root_endophytes_files/figure-html/abundance figures-1.png)<!-- -->

## Alpha diversity


```r
alpha_combined <- rbind(fun_alpha, bac_alpha) %>% 
  subset(select = c("Site", "Storage", "Scion", "Target", "S.chao1", "shannon", "simpson")) %>%
  mutate(kingdom = ifelse(Target == "ITS", "Fungi", "Bacteria")) %>%
  mutate(kingdom = factor(kingdom, levels = c("Fungi", "Bacteria"))) %>%
  mutate(Storage = factor(Storage, levels = c("no", "yes"))) %>%
  rename(Chao1 = S.chao1, Shannon = shannon, Simpson = simpson) %>%
  pivot_longer(
    cols = c("Chao1", "Shannon", "Simpson"), names_to = "measure", values_to = "value"
  )
```

```
# Error in rename(., Chao1 = S.chao1, Shannon = shannon, Simpson = simpson): unused arguments (Chao1 = S.chao1, Shannon = shannon, Simpson = simpson)
```

```r
# alpha_bar <- ggbarplot(
#   data = alpha_combined, x = "Storage", y = "value", 
#   fill = "Site", add = "mean_se", facet.by = c("measure", "kingdom"), 
#   palette = cbPalette, position = position_dodge(0.8), scales = "free_y",
#   ylab = "Mean diversity index", xlab = F, legend = "right"
# ) + guides(fill = guide_legend(title = "Site"))

# ggsave(
#   filename = "alpha_bar.png", plot = alpha_bar, path = "figures/", 
#   height = 12, width = 24, units = "cm"
# )

# alpha_bar

alpha_box <- ggboxplot(
  data = alpha_combined, x = "Site", y = "value", 
  color = "Storage", add = "jitter", facet.by = c("measure", "kingdom"),
  palette = cbPalette, scales = "free_y", legend = "bottom",
  ylab = "Mean diversity index", xlab = "Site"
) #+ guides(color = guide_legend(position = "right"))
```

```
# Error in eval(expr, envir, enclos): object 'alpha_combined' not found
```

```r
ggsave(
  filename = "alpha_box.png", plot = alpha_box, path = "figures/", 
  height = 24, width = 24, units = "cm"
)
```

```
# Error in eval(expr, envir, enclos): object 'alpha_box' not found
```

```r
alpha_box
```

```
# Error in eval(expr, envir, enclos): object 'alpha_box' not found
```

```r
alpha_box_fungi <- ggboxplot(
  data = alpha_combined[alpha_combined$kingdom == "Fungi", ], x = "Site", y = "value", 
  color = "Storage", add = "jitter", facet.by = "measure",
  palette = cbPalette, scales = "free_y", legend = "bottom",
  ylab = "Diversity index", xlab = "Site"
)
```

```
# Error in eval(expr, envir, enclos): object 'alpha_combined' not found
```

```r
ggsave(
  filename = "alpha_box_fungi.png", plot = alpha_box_fungi, path = "figures/", 
  height = 12, width = 24, units = "cm"
)
```

```
# Error in eval(expr, envir, enclos): object 'alpha_box_fungi' not found
```

```r
alpha_box_bacteria <- ggboxplot(
  data = alpha_combined[alpha_combined$kingdom == "Bacteria", ], x = "Site", y = "value", 
  color = "Storage", add = "jitter", facet.by = "measure",
  palette = cbPalette, scales = "free_y", legend = "bottom",
  ylab = "Diversity index", xlab = "Site"
)
```

```
# Error in eval(expr, envir, enclos): object 'alpha_combined' not found
```

```r
ggsave(
  filename = "alpha_box_bacteria.png", plot = alpha_box_bacteria, path = "figures/", 
  height = 12, width = 24, units = "cm"
)
```

```
# Error in eval(expr, envir, enclos): object 'alpha_box_bacteria' not found
```

```r
alpha_box_combined <- ggarrange(
  alpha_box_fungi, alpha_box_bacteria,
  ncol = 1, nrow = 2, labels = c("A", "B"),
  common.legend = T, legend = "bottom"
)
```

```
# Error in eval(expr, envir, enclos): object 'alpha_box_fungi' not found
```

```r
ggsave(
  filename = "alpha_box_combined.png", plot = alpha_box_combined, path = "figures/", 
  height = 20, width = 24, units = "cm"
)
```

```
# Error in eval(expr, envir, enclos): object 'alpha_box_combined' not found
```

## PCA
  

```r
pca_combo_plot <- ggarrange(
  fun_pca_plot, bac_pca_plot,
  ncol = 2, nrow = 1, labels = c("A", "B"),
  common.legend = T, legend = "bottom"
)

ggsave(
  filename = "pca_combo_plot.png", plot = pca_combo_plot, path = "figures/", 
  height = 10, width = 24, units = "cm"
)

pca_combo_plot
```

![](root_endophytes_files/figure-html/PCA figures-1.png)<!-- -->

```r
nmds_combo_plot <- ggarrange(
  fun_nmds_plot, bac_nmds_plot,
  ncol = 2, nrow = 1, widths = c(1.19, 1),
  common.legend = T, legend = "bottom"
)

ggsave(
  filename = "nmds_combo_plot.png", plot = nmds_combo_plot, path = "figures/", 
  height = 13, width = 24, units = "cm"
)

nmds_combo_plot
```

![](root_endophytes_files/figure-html/PCA figures-2.png)<!-- -->

```r
mega_combo_plot <- ggarrange(
  fun_pca_plot, bac_pca_plot,
  fun_nmds_plot, bac_nmds_plot,
  ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"),
  common.legend = T, legend = "bottom"
) + labs(shape = "Storage")

ggsave(
  filename = "mega_combo_plot.png", plot = mega_combo_plot, path = "figures/", 
  height = 25, width = 30, units = "cm"
)

mega_combo_plot
```

![](root_endophytes_files/figure-html/PCA figures-3.png)<!-- -->


```r
# Save environment
save.image("BAC.RData")
```
