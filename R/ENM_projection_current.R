# Make MEAN FULL model projections for current data
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
# 2021-02-23; 2021-03-14: Code tidy-up
# 2022-08-16: Further cleaning of code

# library(fitMaxnet)
# library(terra)
# library(stringr)
# library(Rfast)

#' ENM Projection onto Current Environment
#'
#' @param numThins Integer
#' @param numFolds Integer
#' @param env_folder Character
#' @param taxonTable_filename Character
#'
#' @return
#' @export
#'
#' @examples
ENM_projection_current <- function(numThins = 5,
                                   numFolds = 5,
                                   env_folder = "",
                                   taxonTable_filename = "")
{
  ###############################################################################
  # Set critical parameters to drive the process

  # Number of cross-validation folds
  numFolds <- 5

  # Number of replicate thinning runs
  numThins <- 5

  env_folder <- "/home/peterw/Data_and_Projects/Climate/CHELSA/Current/eastOZ/bioclim/"

  taxonTable_filename <- "/home/peterw/Downloads/workflow_tests/taxonTable_demonstration.csv"
  taxonTable <- read.csv(taxonTable_filename, stringsAsFactors = FALSE)
  rownames(taxonTable) <- taxonTable$taxon
  taxonList <- taxonTable$taxon

  cat("MaxEnt ENMs fitted by Aggregation Method: Projection onto current climate data:\n")
  cat("===============================================================================\n")

  basePath <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"

  cat("  Preparing projection data set")
  fitMaxnet::prepProjData(env_folder)
  cat("done.\n")

  cat("  Start of projection runs:\n")
  cat("  -------------------------\n")

  for (thisTaxon in taxonList)
  {
    this_Taxon <- gsub(" ", "_", thisTaxon)
    thisRegVal <- taxonTable[thisTaxon, "bestReg"]
    thisRegStr <- paste0("reg_", gsub(".", "_", thisRegVal, fixed = TRUE))
    aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "best_AggLevel"])

    cat("  ", thisTaxon, ":\n", sep = "")
    cat("    Aggregation level:", taxonTable[thisTaxon, "best_AggLevel"], "\n")
    cat("    Regularisation value:", thisRegVal, "\n")

    if (thisRegVal > 0)
    {
      cat("    Projecting models\n")
      # Generate vector of model file names
      theModelFiles <- paste0(paste0(basePath, thisTaxon, "/", aggStr, "/thin_", 1:numThins), "/Fold_", rep(1:numFolds, each = numThins), "/", thisRegStr, "/", this_Taxon,  "_Fold_", rep(1:numFolds, each = numThins), "_", thisRegStr, ".Rd")

      for (thisFile in theModelFiles)
      {
        fitMaxnet::projectMaxnet(thisTaxon,
                                 thisFile,
                                 type = "cloglog",
                                 baseOutputPath = dirname(thisFile))
      }

      cat("\n----------------------------------------------------------------\n")
    }
    else
    {
      cat("Skipping this taxon because regularisation = 0...nothing can be done.\n")
      cat("\n----------------------------------------------------------------\n")
    }
  }


  cat("  Making mean and std dev rasters:\n")

  if (exists("projData")) rm(projData)
  gc()

  rasTemplate <- terra::rast(list.files(env_folder, "*.tif", full.names = TRUE)[1])

  for (thisTaxon in taxonList)
  {
    cat("  ", thisTaxon, ": ")
    this_Taxon <- gsub(" ", "_", thisTaxon)
    thisRegVal <- taxonTable[thisTaxon, "bestReg"]
    thisRegStr <- paste0("reg_", gsub(".", "_", thisRegVal, fixed = TRUE))
    aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "best_AggLevel"])

    if (thisRegVal > 0)
    {
      outFolder <- paste0(basePath, thisTaxon, "/ENM_Results/Current/")
      if (!dir.exists(outFolder)) dir.create(outFolder, recursive = TRUE)

      theRasFiles <- paste0(paste0(basePath, thisTaxon, "/", aggStr, "/thin_", 1:numThins), "/Fold_", rep(1:numFolds, each = numThins), "/", thisRegStr, "/", this_Taxon, "_projection.tif")

      rasData <- as.matrix(terra::rast(theRasFiles))
      meanVec <- Rfast::rowmeans(rasData)
      meanRas <- rasTemplate
      meanRas[] <- meanVec
      terra::writeRaster(meanRas, filename = paste0(outFolder,  this_Taxon, "_mean.tif"), overwrite = TRUE)

      sdVec <- Rfast::rowVars(rasData, std = TRUE)
      sdRas <- rasTemplate
      sdRas[] <- sdVec
      terra::writeRaster(sdRas, filename = paste0(outFolder, this_Taxon, "_stdDev.tif"), overwrite = TRUE)
      taxonTable[thisTaxon, "meanProjection_current"] <- as.character(Sys.Date())
      cat("done\n")
    }
    else
    {
      cat("Skipping this taxon because regularisation = 0...nothing can be done.\n")
    }

  }

  write.csv(taxonTable, taxonTable_filename, row.names = FALSE)

  cat("\n*** End of processing\n")
}
