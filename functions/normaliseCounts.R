#' Normalise count data
#' Normalise count data using qPCR theoretical copy number data
#' by multiplying the count data by the ratio of copy number:total OTUs
#' using the formula adjustedOTUs = OTUcount * copyNumber / totalOTUs.
#' `countData` is a matrix of OTU counts, `abundance` is a data frame of
#' qPCR copy number data, `ID` is the column name of the sample ID in both dataframes.
#' 'countData' is overwritten with the normalised data.
normaliseCountData <- function(data) {
  # Calculate total OTUs per sample
  totalOTUs <- apply(data$countData, 2, sum)
  totalOTUs <- data.frame("ID" = names(totalOTUs), "TotalOTUs" = totalOTUs, row.names = NULL)
  
  # Merge total OTUs with abundance data
  data[["abundance"]] <- merge(data$abundance, totalOTUs, by = "ID")

  # Calculate adjusted total OTUs
  data[["abundance"]] <- data$abundance %>%
    mutate(AdjustedTotalOTUs = MeanAdjustedTCN / TotalOTUs)

  # Save raw count data
  data$rawCountData <- data$countData

  # Multiply count data by adjusted total OTUs for each sample
  abundance = t(data.frame("AdjustedTotalOTUs" = data$abundance$AdjustedTotalOTUs, row.names = data$abundance$ID))
  data$countData<- data$rawCountData %>%
    mutate(across(everything(), ~ . * abundance[, cur_column()]))

  return(data)
}

# ubiome_BAC <- normaliseCountData(ubiome_BAC)
# ubiome_FUN <- normaliseCountData(ubiome_FUN)

# identical(ubiome_FUN$rawCountData, ubiome_BAC$countData)