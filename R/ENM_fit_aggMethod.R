# Fit MaxEnt models using aggregation method
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
# 2022-08-16: Prepare a generic version of MaxEnt fitting using the aggregation
# method to adjust for sampling bias + switch to using R-package terra in place
# of pkg raster.

# library(fitMaxnet)
# library(terra)
# library(sf)
# library(gdata)
# library(doParallel)
# library(fasterize)

#' ENM fit aggMethod
#'
#' @param numThin Number of replicate applications of the aggregation thinning method
#' @param numFolds Number of cross-validation folds to be used in each replicate
#' @param occ_CRS Integer. EPSG integer code for the CRS for the occurrence data. Default is 4326, lat/long coordinates in WGS84 datum
#' @param env_CRS Integer. EPSG integer code for the CRS for the environmental raster data. Default is 4326, lat/long coordinates in WGS84 datum
#' @param occ_X Characater. Name of the field in occurrence data for the X-coordinate relevant to the CRS specified in \emph{occ_CRS}. Default is "longitude"
#' @param occ_Y Characater. Name of the field in occurrence data for the Y-coordinate relevant to the CRS specified in \emph{occ_CRS}. Default is "latitude"
#' @param aggFactorSet Integer vector. Set of aggregation values to be applied to the occurrence data
#' @param regSequence Numeric vector. Sequence of MaxEnt regularisation values to be used in fitting MaxEnt models
#' @param taxonTable_filename Character. Path the the file holding the metadata table which will supply necessary informaiton to drive the model fitting process
#' @param outputFolder Character. Path the the base folder below which output folders will be systematically created to store output
#' @param occBaseFolder Character. Path the the base folder below which occurrence data files will be found
#' @param envDataPath Character. Path to the folder storing the raster layers chosen to be the covariates
#'
#' @return
#' @export
#'
#' @examples
#'
ENM_fit_aggMethod <- function(numThin = 5, numFolds = 5, occ_CRS = 4326, env_CRS = 4326,
                              occ_X = "longitude",
                              occ_Y = "latitude",
                              aggFactorSet = c(2, 4, 8, 16),
                              regSequence = 1:10,
                              taxonTable_filename = "",
                              outputFolder = "",
                              occBaseFolder = "",
                              envDataPath = "")
{
  ####################################################################################
  # Set the number of CPU cores to use in parallel computations
  totalCores <- parallel::detectCores()
  numCores <- ifelse(totalCores > 2, totalCores - 2, 1)
  ####################################################################################

  # numThin <- 5
  #
  # numFolds <- 5
  #
  # occ_CRS <- 4326 # WGS84 ellipsiodial lat/long
  # #env_CRS <- 3577 # Australian Albers Equal Area projection
  # env_CRS <- 4326
  #
  # occ_X <- "longitude"
  # occ_Y <- "latitude"
  #
  # aggFactorSet <- c(2, 4, 8, 16) #, 32) #, 64)
  #
  # regSequence <- 1:10

  #source("/home/peterw/Nyctimene/MU Projects/MaxEnt_agg_method_workflow/list_files.R")

  #taxonTable_filename <- "/home/peterw/Downloads/workflow_tests/taxonTable_demonstration.csv"
  taxonTable <- read.csv(taxonTable_filename, stringsAsFactors = FALSE)
  rownames(taxonTable) <- taxonTable$taxon
  taxonList <- taxonTable$taxon


  # outputFolder <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"
  # occBaseFolder <- "/home/peterw/Downloads/workflow_tests/occ_data/"
  #
  # envDataPath <- "/home/peterw/Data_and_Projects/Climate/CHELSA/Current/eastOZ/bioclim/"

  # Provide the variable names of any environmental covariates you wish to remove
  # from the process. Comment out the following line to retain all env covars
  # found in the folder at envDataPath
  #excludedVars <-  c("bioclim_08_albers", "bioclim_09_albers")

  dataFileSet <- list.files(envDataPath, "*.tif", full.names = TRUE)

  # if (length(excludedVars) > 0)
  # {
  #   dropInd <- unlist(lapply(1:length(excludedVars), function(i) {grep(excludedVars[i], dataFileSet)}))
  #   dataFileSet <- dataFileSet[-dropInd]
  # }

  # Load env covariates into a "stack"
  envStack <- terra::rast(dataFileSet)

  # Take a copy of an environmental covariate raster to be used later as a
  # template to generate output rasters
  rasTemplate <- envStack[[1]]

  ###############################################################################
  fitModel <- function(regIndex)
  {
    # NOTE: the variables 'newFolder', 'thisTaxon', 'thisFold', 'trainSWD',
    # 'trainBKG' are a global variables generated in the main script

    regFolder <- paste0(newFolder, "/reg_", gsub(".", "_", regSequence[regIndex], fixed = TRUE))
    if (!dir.exists(regFolder)) dir.create(regFolder)

    #cat("         start model fit...")
    thisModel <- fitMaxnet::fit_maxnet(taxonName = thisTaxon,
                                       replTag = paste0("Fold_", thisFold),
                                       baseOutputPath = regFolder,
                                       predVar = c(rep(1, nrow(trainSWD)), rep(0, nrow(trainBKG))),
                                       envData = rbind(trainSWD[, -c(1:3)], trainBKG[, -c(1:3)]),
                                       regMult = regSequence[regIndex])
    #cat("done\n")
    return(TRUE)
  }


  ###############################################################################
  #thisTaxon <- taxonList[1]
  for (thisTaxon in taxonList)
  {
    cat(thisTaxon, "\n")
    this_Taxon <- gsub(" ", "_", thisTaxon, fixed = TRUE)
    cat("    Preliminary work...")
    taxonOutputFolder <- paste0(outputFolder, thisTaxon)
    if (!dir.exists(taxonOutputFolder)) dir.create(taxonOutputFolder, recursive = TRUE)

    bufferDist <- taxonTable[thisTaxon, "bufferDist_km"]

    numBkgPts <- taxonTable[thisTaxon, "numBkgPts"]

    if ((bufferDist == 0) | (numBkgPts == 0))
      cat("Skipping taxon: bufferDist or numBkgPts == 0\n\n")
    else
    {
      ########## Change this line to suit your file naming convention
      occData <- read.csv(paste0(occBaseFolder, thisTaxon, "/", this_Taxon, "_filtered_merged_data.csv"))

      badCoordInd <- which(is.na(occData[, occ_X]))
      if (length(badCoordInd) > 0) occData <- occData[-badCoordInd, ]

      occ_sf <- sf::st_as_sf(occData, coords = c(occ_X, occ_Y), crs = occ_CRS)

      if (occ_CRS != env_CRS)
        occ_sf <- sf::st_transform(occ_sf, crs = env_CRS)

      occInd <- terra::cellFromXY(rasTemplate, sf::st_coordinates(occ_sf))

      allEnvData <- terra::extract(envStack, occInd)
      naInd <- which(is.na(rowSums(allEnvData)))
      if (length(naInd) > 0)
      {
        allEnvData <- allEnvData[-naInd, ]
        occ_sf <- occ_sf[-naInd, ]
        occInd <- occInd[-naInd]
        occData <- occData[-naInd, ]
      }
      cat("done\n")

      for (aggFactor in aggFactorSet)
      {
        aggFolder <- paste0(taxonOutputFolder, "/aggFactor_", aggFactor)
        if (!dir.exists(aggFolder)) dir.create(aggFolder, recursive = TRUE)

        newRas <- terra::aggregate(rasTemplate, aggFactor)
        newOccCells <- terra::cellFromXY(newRas, sf::st_coordinates(occ_sf))
        newTab <- rle(sort(newOccCells))

        if (length(newTab$values) > 25)
        {
          # Thin using the tabulated data:
          stuff <- data.frame(rowInd = 1:nrow(occData), cellInd = newOccCells, selected = rep(0, nrow(occData)))
          targetCellInd <- newTab$values[newTab$lengths > 1]

          # Mark singleton cells as selected
          rowsToMark <- match(newTab$values[newTab$lengths == 1], stuff$cellInd)
          stuff$selected[rowsToMark] <- 1

          for (thisThin in 1:numThin)
          {
            # Reset cell selection
            stuff$selected <- 0
            stuff$selected[rowsToMark] <- 1

            thinFolder <- paste0(aggFolder, "/thin_", thisThin)
            if (!dir.exists(thinFolder)) dir.create(thinFolder)

            for (thisTarget in targetCellInd)
            {
              theseRows <- stuff$rowInd[stuff$cellInd == thisTarget]
              stuff$selected[sample(theseRows, 1)] <- 1
            }

            cat("Agg factor =", aggFactor, " : Number of occupied cells =", length(newTab$values), ":   Number of retained points =", sum(stuff$selected), "\n")
            envData_thinned <- data.frame(X = occData[stuff$rowInd[stuff$selected == 1], occ_X],
                                          Y = occData[stuff$rowInd[stuff$selected == 1], occ_Y],
                                          allEnvData[stuff$rowInd[stuff$selected == 1], ])
            occSWD <- data.frame(species = rep(thisTaxon, length(newTab$values)),
                                 envData_thinned,
                                 stringsAsFactors = FALSE)

            badRows <- unique(which(is.na(occSWD), arr.ind = TRUE)[ , 1])
            if (length(badRows) > 0) occSWD <- occSWD[-badRows, ]
            outputSWDname <- paste0(thinFolder, "/", this_Taxon, "_SWD.csv")
            write.csv(occSWD, outputSWDname, row.names = FALSE)

            cat("      Computing buffer...")
            # Make a bounds polygon buffering occurrence points by bufferDist km
            envData_thinned_sf <- sf::st_as_sf(envData_thinned, coords = c(1, 2), crs = 4283)

            buffPoly <- fitMaxnet::bufferPoints(envData_thinned_sf, bufferDist)
            cat("done\n")

            cat("      Determining available cells for background data extraction...")
            bkgRas <- fasterize::fasterize(sf::st_transform(buffPoly, env_CRS), raster::raster(rasTemplate))
            bkgPool <- which(!is.na(bkgRas[]))
            bkg_cells <- sample(bkgPool, numBkgPts)
            bkg_xy <- terra::xyFromCell(rasTemplate, bkg_cells)
            bkg_envData <- matrix(0, length(bkg_cells), length(dataFileSet))
            cat("done\n")

            cat("      Extracting background environmental data...")
            bkg_envData <- terra::extract(envStack, bkg_cells)

            colnames(bkg_envData) <- gsub(".tif$", "", basename(dataFileSet))

            bkg_envData <- cbind(bkg_xy, bkg_envData)
            cat("done\n")

            cat("      Saving Background SWD file...")

            bkgData <- data.frame(species = rep("background", length(bkg_cells)),
                                  bkg_envData,
                                  stringsAsFactors = FALSE)

            # Filter any cells which unfortunately fall into an NA cells in at least one raster layer
            badRows <- unique(which(is.na(bkgData), arr.ind = TRUE)[ , 1])
            if (length(badRows) > 0) bkgData <- bkgData[-badRows, ]
            outputSWDBKGname <- paste0(thinFolder, "/", this_Taxon, "_SWD_BKG.csv")
            write.csv(bkgData, outputSWDBKGname, row.names = FALSE)
            cat("done\n")
          }

          cat("\n   Make k-fold sample sets: k =", numFolds, "for", numThin, "thinned occurrence data sets: ")
          for (thisThin in 1:numThin)
          {
            thinFolder <- paste0(aggFolder, "/thin_", thisThin)

            swd <- read.csv(paste0(thinFolder, "/", this_Taxon, "_SWD.csv"), stringsAsFactors = FALSE)
            ans <- fitMaxnet::makeFolds(nrow(swd), numFolds)
            saveRDS(ans, file = paste0(thinFolder, "/SWD_folds.rds"))

            bkg <- read.csv(paste0(thinFolder, "/", this_Taxon, "_SWD_BKG.csv"), stringsAsFactors = FALSE)
            ans <- fitMaxnet::makeFolds(nrow(bkg), numFolds)
            saveRDS(ans, file = paste0(thinFolder, "/BKG_folds.rds"))
          }

          cat("completed\nStart of model fitting:\n")

          for (thisThin in 1:numThin)
          {
            cat("    Thinned dataset ", thisThin, ":\n", sep ="")
            swdData <- read.csv(paste0(aggFolder, "/thin_", thisThin, "/", this_Taxon, "_SWD.csv"), stringsAsFactors = FALSE)
            swdFolds <- readRDS(paste0(aggFolder, "/thin_", thisThin, "/SWD_folds.rds"))

            # For background data
            bkgData <- read.csv(paste0(aggFolder, "/thin_", thisThin, "/", this_Taxon, "_SWD_BKG.csv"), stringsAsFactors = FALSE)
            bkgFolds <- readRDS(paste0(aggFolder, "/thin_", thisThin, "/BKG_folds.rds"))

            for (thisFold in 1:numFolds)
            {
              cat("      Fold", thisFold, ":\n")

              # Assemble test and training split of data for this fold
              trainSWD <- swdData[-swdFolds[[thisFold]], ]
              testSWD <- swdData[swdFolds[[thisFold]], ]
              trainBKG <- bkgData[-bkgFolds[[thisFold]], ]
              testBKG <- bkgData[bkgFolds[[thisFold]], ]

              newFolder <- paste0(aggFolder, "/thin_", thisThin, "/Fold_", thisFold)
              if (!dir.exists(newFolder)) dir.create(newFolder)

              if (nrow(testSWD) > 0)
              {
                # Fit a sequence of models
                cat("       Start maxnet model fitting along regularisation set\n")

                cl <- parallel::makeCluster(numCores)
                doParallel::registerDoParallel(cl)

                foreach::foreach (j = 1:length(regSequence), .packages = "fitMaxnet") %dopar%
                  {
                    fitModel(j)
                  }

                parallel::stopCluster(cl)

                taxonTable[thisTaxon, "model_lastUpdated"] <- as.character(Sys.Date())
              }
              else
                cat("       Skipping fold", thisFold, ": no data\n")
            }
          }
        }
        else
        {
          cat("**** Too few occupied cells to continue\n")
        }

      }
    }

    write.csv(taxonTable, taxonTable_filename, row.names = FALSE)
  }

  cat("*** End of processing\n")
}
