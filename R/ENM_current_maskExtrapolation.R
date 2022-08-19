# SoS PLP Project
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
# 2021-09-02; 2021-10-13: Adapted for a full run through the taxon list
# 2022-08-16: Code clean-up and conversion to using package 'terra'

# library(fitMaxnet)
# library(terra)
# library(Rfast)

#' ENM mask extrapolation areas under current environment
#'
#' @param numThins Integer
#' @param numFolds Integer
#' @param basePath Character
#' @param env_folder Character
#' @param taxonTable_filename Character
#'
#' @return
#' @export
#'
#' @examples
ENM_current_maskExtrapolation <- function(numThins = 5,
                                          numFolds = 5,
                                          basePath = "",
                                          env_folder = "",
                                          taxonTable_filename = "")
{
  # basePath <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"
  # env_folder <- "/home/peterw/Data_and_Projects/Climate/CHELSA/Current/eastOZ/bioclim/"

  cat("MaxEnt ENMs fitted by Aggregation Method: Mask extrapolation regions from projection onto current climate data:\n")
  cat("===============================================================================================================\n")


  cat("  Preparing projection data set...")
  fitMaxnet::prepProjData(env_folder)
  cat("done.\n")

  taxonTable <- read.csv(taxonTable_filename, stringsAsFactors = FALSE)
  rownames(taxonTable) <- taxonTable$taxon
  taxonList <- taxonTable$taxon

  # numThins <- 5
  # numFolds <- 5

  cat("  Start of masking process:\n")

  for (thisTaxon in taxonList)
  {
    cat("   ", thisTaxon, ":\n")
    this_Taxon <- gsub(" ", "_", thisTaxon)

    thisRegVal <- taxonTable[thisTaxon, "bestReg"]
    thisRegStr <- paste0("reg_", gsub(".", "_", thisRegVal, fixed = TRUE))

    aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "best_AggLevel"])

    if (thisRegVal > 0)
    {
      for (thisThin in 1:numThins)
      {
        cat("     Thin", thisThin, ": ")

        for (thisFold in 1:numFolds)
        {
          cat("Fold_", thisFold, "  ", sep = "")

          modelFile <- paste0(basePath, thisTaxon, "/", aggStr, "/thin_", thisThin, "/Fold_", thisFold, "/", thisRegStr, "/", this_Taxon, "_Fold_", thisFold, "_", thisRegStr, ".Rd")

          modelRas <- paste0(basePath, thisTaxon, "/", aggStr, "/thin_", thisThin, "/Fold_", thisFold, "/", thisRegStr, "/", this_Taxon, "_projection.tif")

          fitMaxnet::maskExtrapolation(modelFile,
                                       modelRas,
                                       maskOutpath = paste0(basePath, thisTaxon, "/", aggStr, "/thin_", thisThin, "/Fold_", thisFold, "/", thisRegStr),
                                       makePlots = FALSE)
        }

        cat("\n")
      }
    }
    else
    {
      cat("Skipping this taxon because regularisation = 0...nothing can be done.\n")
    }

  }

  cat("  \nMaking mean and std dev rasters:\n")

  if (exists("projData")) rm(projData)
  gc()

  rasTemplate <- terra::rast(list.files(env_folder, "*.tif", full.names = TRUE)[1])

  for (thisTaxon in taxonList)
  {
    cat("      ", thisTaxon, ": ")
    this_Taxon <- gsub(" ", "_", thisTaxon)
    thisRegVal <- taxonTable[thisTaxon, "bestReg"]
    thisRegStr <- paste0("reg_", gsub(".", "_", thisRegVal, fixed = TRUE))
    aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "best_AggLevel"])

    if (thisRegVal > 0)
    {
      outFolder <- paste0(basePath, thisTaxon, "/ENM_Results/Current_masked/")
      if (!dir.exists(outFolder)) dir.create(outFolder, recursive = TRUE)

      theRasFiles <- paste0(paste0(basePath, thisTaxon, "/", aggStr, "/thin_", 1:numThins), "/Fold_", rep(1:numFolds, each = numThins), "/", thisRegStr, "/", this_Taxon, "_projection_masked.tif")

      rasData <- as.matrix(terra::rast(theRasFiles))
      meanVec <- Rfast::rowmeans(rasData)
      meanRas <- rasTemplate
      meanRas[] <- meanVec
      terra::writeRaster(meanRas, filename = paste0(outFolder,  this_Taxon, "_mean_masked.tif"), overwrite = TRUE)

      sdVec <- Rfast::rowVars(rasData, std = TRUE)
      sdRas <- rasTemplate
      sdRas[] <- sdVec
      terra::writeRaster(sdRas, filename = paste0(outFolder, this_Taxon, "_stdDev_masked.tif"), overwrite = TRUE)
      taxonTable[thisTaxon, "meanProjection_current_masked"] <- as.character(Sys.Date())
      cat("done\n")
    }
    else
    {
      cat("Skipping this taxon because regularisation = 0...nothing can be done.\n")
    }

  }
}
