# Model projections onto FUTURE data
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
# 2022-08-16: Code clean-up

# library(fitMaxnet)
# library(stringr)
# library(parallel)

ENM_future_projections <- function(numThins = 5,
                                   numFolds = 5,
                                   gcm_list = "",
                                   rcp_ssp_list = "",
                                   time_steps = "",
                                   base_path = "",
                                   base_env_path = "",
                                   taxonTable_filename = "")
{
  ###############################################################################
  # Set critical parameters to drive the process

  # Number of cross-validation folds
  #numFolds <- 5

  # Number of replicate thinning runs
  # numThins <- 5
  #
  # gcm_list <- c("CCCMA31", "CSIRO-MK30", "ECHAM5", "MIROC32")
  #
  # rcp_ssp_list <- c("rcp45", "rcp85")
  #
  # time_steps <- c("2050", "2070")

  ###############################################################################

  #base_env_path <- "/home/peterw/Data_and_Projects/Climate/CHELSA/AR5/"

  taxonTable <- read.csv(taxonTable_filename, stringsAsFactors = FALSE)
  rownames(taxonTable) <- taxonTable$taxon
  taxonList <- taxonTable$taxon

  #base_path <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"

  cat("MaxEnt ENMs by Aggregation Method: Mask extrapolation and project onto FUTURE climate data:\n")
  cat("-------------------------------------------------------------------------------------------\n")

  for (thisGCM in gcmList)
  {
    cat(thisGCM, ":\n")
    for (thisTime in timeSteps)
    {
      cat("    Time:", thisTime, "\n")
      cat("       Preparing projection dataset...")
      fitMaxnet::prepProjData(paste0(base_env_path, thisGCM, "/", thisTime))
      cat("done\n")

      for (thisTaxon in taxonList)
      {
        cat("          Processing", thisTaxon, ": ")

        this_Taxon <- gsub(" ", "_", thisTaxon)
        thisRegVal <- taxonTable[thisTaxon, "bestReg"]
        thisRegStr <- paste0("reg_", gsub(".", "_", thisRegVal, fixed = TRUE))
        aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "best_AggLevel"])

        if (thisRegVal > 0)
        {
          taxonFolder <- paste0(base_path, thisTaxon, "/", aggStr, "/")

          for (thisThin in 1:numThins)
          {
            theseModels <- paste0(taxonFolder, "thin_", thisThin, "/Fold_", 1:numFolds, "/reg_", thisRegVal, "/", this_Taxon,  "_Fold_", 1:numFolds, "_reg_", thisRegVal, ".Rd")
            if (!all(file.exists(theseModels))) stop("Missing model files")

            rasFiles <- paste0(taxonFolder, "thin_", thisThin,"/future/", this_Taxon,"_projection_", thisGCM,"_", thisTime,"_thin_", thisThin, "_fold_", 1:numFolds, ".tif")
            if (!all(file.exists(rasFiles))) stop("Missing raster files")

            ans <- parallel::mclapply(1:numFolds,
                                      function(i) {
                                        fitMaxnet::maskExtrapolation(theseModels[i],
                                                                     projRas = rasFiles[i],
                                                                     makePlots = FALSE,
                                                                     saveMask = FALSE)
                                      },
                                      mc.cores = 4)
          }
        }

        cat("done\n")
      }
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
      outFolder <- paste0(base_path, thisTaxon, "/ENM_Results/Future_masked/")
      if (!dir.exists(outFolder)) dir.create(outFolder, recursive = TRUE)

      theRasFiles <- paste0(paste0(base_path, thisTaxon, "/", aggStr, "/thin_", 1:numThins), "/Fold_", rep(1:numFolds, each = numThins), "/", thisRegStr, "/", this_Taxon, "_projection_masked.tif")

      rasData <- as.matrix(terra::rast(theRasFiles))
      meanVec <- Rfast::rowmeans(rasData)
      meanRas <- rasTemplate
      meanRas[] <- meanVec
      terra::writeRaster(meanRas, filename = paste0(outFolder,  this_Taxon, "_mean_masked.tif"), overwrite = TRUE)

      sdVec <- Rfast::rowVars(rasData, std = TRUE)
      sdRas <- rasTemplate
      sdRas[] <- sdVec
      terra::writeRaster(sdRas, filename = paste0(outFolder, this_Taxon, "_stdDev_masked.tif"), overwrite = TRUE)
      taxonTable[thisTaxon, "meanProjection_future_masked"] <- as.character(Sys.Date())
      cat("done\n")
    }
    else
    {
      cat("Skipping this taxon because regularisation = 0...nothing can be done.\n")
    }
  }



  cat("*** End of processing\n")
}
