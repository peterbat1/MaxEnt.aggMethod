# Collate information on covariates used in ENMs
#
# First developed progressively between January and November 2021 for NSW Dept
# Planning and Environment PLP SoS Project and Natural Resources Council FMIP
# project
#
# Peter D. Wilson
# Adjunct Fellow
# Dept. of Biological Sciences
# Faculty of Science and Engineering
# Macquarie University, Sydney, Australia
#
# 2021-03-11
# 2022-08-16: Code clean-up

#library(fitMaxnet)

#' Collate information on covariates used in ENMs
#'
#' @param numThins
#' @param numFolds
#' @param basePath Character. Path to the base folder below which model files for each taxon will be found
#' @param envCovarPath Character. Path to the folder in which the covariate raster used for model fitting are found
#' @param taxonTable_fileName
#'
#' @return
#' @export
#'
#' @examples
Collate_covar_importance <- function(numThins = 5,
                                     numFolds = 5,
                                     basePath = "",
                                     envCovarPath = "",
                                     taxonTable_fileName = "")
{
  taxonTable <- read.csv(taxonTable_fileName, stringsAsFactors = FALSE)
  rownames(taxonTable) <- taxonTable$taxon
  taxonList <- taxonTable$taxon

  # numThins  <- 5
  # numFolds <- 5

  # Path to the base folder below which model files for each taxon will be found
  #basePath <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"

  # Path to the folder in which the covariate raster used for model fitting are found
  #envCovarPath <- "/home/peterw/Data_and_Projects/Climate/CHELSA/Current/eastOZ/bioclim"

  # Gather a list of covariates available for modelling
  covFiles <- gsub(".tif", "", list.files(envCovarPath), fixed = TRUE)
  if (length(covFiles) == 0) stop("Environmental covariate files not found")

  # Comment-out the line below if no covars are excluded, or adjust search pattern to suit
  # Drop bioclim_08 and bioclim_09
  # covFiles <- covFiles[-grep("_08_|_09_", covFiles)]

  #####################################################################
  meanResults <- data.frame(taxon = taxonList, matrix(0, length(taxonList), length(covFiles), dimnames = list(taxonList, covFiles)))
  minResults <- data.frame(taxon = taxonList, matrix(0, length(taxonList), length(covFiles), dimnames = list(taxonList, covFiles)))
  maxResults <- data.frame(taxon = taxonList, matrix(0, length(taxonList), length(covFiles), dimnames = list(taxonList, covFiles)))

  cat("MaxEnt models fitted using Aggregation Method: Gather variable importance information\n=====================================================================================\n")

  for (thisTaxon in taxonList)
  {
    cat("  ", thisTaxon, ": ")
    this_Taxon <- gsub(" ", "_", thisTaxon, fixed = TRUE)

    bestReg <- taxonTable[thisTaxon, "bestReg"]

    aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "best_AggLevel"])

    taxonFolder <- paste0(basePath, thisTaxon, "/", aggStr)

    if (bestReg > 0)
    {
      varContribSums <- rep(0, length(covFiles))
      varContribMins <- rep(999, length(covFiles))
      varContribMaxs <- rep(-999, length(covFiles))

      for (thisThin in 1:numThins)
      {
        occSWD <- read.csv(paste0(taxonFolder, "/thin_", thisThin, "/", this_Taxon, "_SWD.csv"), stringsAsFactors = FALSE)
        bkgSWD <- read.csv(paste0(taxonFolder, "/thin_", thisThin, "/", this_Taxon, "_SWD_BKG.csv"), stringsAsFactors = FALSE)

        for (thisFold in 1:numFolds)
        {
          load(paste0(taxonFolder, "/thin_", thisThin, "/Fold_", thisFold, "/reg_", bestReg, "/", this_Taxon, "_Fold_", thisFold, "_reg_", bestReg, ".Rd"))

          modelVarContribs <- fitMaxnet::varImportance(maxnet_model, occSWD, bkgSWD, responseType = "cloglog")
          #print(modelVarContribs)

          varContribSums <- varContribSums + modelVarContribs
          #print(varContribSums)

          newMaxInd <- which(modelVarContribs > varContribMaxs)
          if (length(newMaxInd) > 0) varContribMaxs[newMaxInd] <- modelVarContribs[newMaxInd]

          newMinInd <- which(modelVarContribs < varContribMins)
          if (length(newMinInd) > 0) varContribMins[newMinInd] <- modelVarContribs[newMinInd]
        }
      }

      # Store as mean contribution across all replicate model fits
      varContribSums <- varContribSums/(numThins * numFolds)

      meanResults[thisTaxon, 2:ncol(meanResults)] <- varContribSums
      maxResults[thisTaxon, 2:ncol(maxResults)] <- varContribMaxs
      minResults[thisTaxon, 2:ncol(minResults)] <- varContribMins
      cat("done\n")
    } else
    {
      cat("No MaxEnt model present - skipped\n")
    }
  }

  write.csv(meanResults, "/home/peterw/Nyctimene/PLP/Metadata/MEAN_variable_importance_summary_refit.csv", row.names = FALSE)
  write.csv(minResults, "/home/peterw/Nyctimene/PLP/Metadata/MINIMUM_variable_importance_summary_refit.csv", row.names = FALSE)
  write.csv(maxResults, "/home/peterw/Nyctimene/PLP/Metadata/MAXIMUM_variable_importance_summary_refit.csv", row.names = FALSE)

  cat("*** End of processing\n")
}
