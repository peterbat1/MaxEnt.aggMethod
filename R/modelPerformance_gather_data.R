# Fit MaxEnt models using aggregation method
#
# First developed progressively between January and November 2021 for NSW Dept
# Planning and Environment PLP SoS Project and Natural Resources Council FMIP
# project
#
# Peter D. Wilson
#
# 2020-10-09; 2021-03-14: Code tidy-up & generalisation
# 2022-08-16: Code clean-up

# library(fitMaxnet)
# library(hmeasure)
# library(enmSdm)

#' Title
#'
#' @param numThins Integer
#' @param numFolds Integer
#' @param aggregation_levels Integer vector
#' @param regularisation_set Numeric vector
#' @param taxonTable_filename Characater
#'
#' @return
#' @export
#'
#' @examples
ENM_gather_performance_data <- function(numThins = 5,
                                        numFolds = 5,
                                        aggregation_levels = c(2, 4, 8, 16),
                                        regularisation_set = 1:10,
                                        taxonTable_filename)
{
  cat("MaxEnt ENMs fitted using the Agregation Method: Compute model performance statistics\n
    -------------------------------------------------------------------------------------\n")

  ###############################################################################
  # Set critical parameters to drive the process

  # Number of cross-validation folds
  #numFolds <- 5

  # Number of replicate thinning runs
  #numThins <- 5

  # Sequence of regularisation values to be used in the search for the best
  # performing model
  #regularisation_set <- 1:10

  # Sequence of aggregation levels or factors
  #aggregation_levels <- c(2, 4, 8, 16)

  ###############################################################################

  # Make a helpful look-up lists
  regLabels <- paste0("reg_", gsub(".", "_", as.character(regularisation_set), fixed = TRUE))
  names(regularisation_set) <- regLabels

  aggFactorSet <- paste0("aggFactor_", aggregation_levels)

  taxonTable <- read.csv(taxonTable_filename, stringsAsFactors = FALSE)
  rownames(taxonTable) <- taxonTable$taxon
  taxonList <- taxonTable$taxon

  basePath <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"

  for (thisTaxon in taxonList)
  {
    cat("\n#####", thisTaxon, "#####\n")

    for (aggFactor in aggFactorSet)
    {
      cat("\n  ", thisTaxon, "at aggFactor = ", aggFactor, ":\n")
      this_Taxon <- gsub(" ", "_", thisTaxon, fixed = TRUE)

      taxonPath <- paste0(basePath, thisTaxon, "/")

      extraMeasures_test <- NULL
      extraMeasures_train <- NULL

      if (dir.exists(paste0(taxonPath, aggFactor, "/thin_1")))
      {
        for (thisThin in 1:numThins)
        {
          cat("thin", thisThin, "\n")
          thinPath <- paste0(taxonPath, aggFactor, "/thin_", thisThin)

          SWD_occ <- read.csv(paste0(thinPath, "/", this_Taxon, "_SWD.csv"))
          SWD_bkg <- read.csv(paste0(thinPath, "/", this_Taxon, "_SWD_BKG.csv"))

          #Load occ SWD folds data object. #### NOTE: We do not load the bkg SWD folds
          #object since, according to best model performance evaluation practice, we
          #evaluate model performance against the full background
          SWD_folds <- readRDS(paste0(thinPath, "/SWD_folds.rds"))

          #thisFold <- 1
          for (thisFold in 1:numFolds)
          {
            #cat("  fold", thisFold, ": ")
            foldPath <- paste0(thinPath, "/Fold_", thisFold)

            #thisReg = regLabels[1]
            for (thisReg in 1:length(regLabels))
            {
              thisRegLabel <- regLabels[thisReg]
              #cat("  ", thisRegLabel)
              #print(paste0(foldPath, "/", thisReg, "/", this_Taxon, "_projection_masked.tif"))

              if (file.exists(paste0(foldPath, "/", thisRegLabel, "/", this_Taxon, "_Fold_", thisFold, "_", thisRegLabel,".Rd")))
              {
                load(paste0(foldPath, "/", thisRegLabel, "/", this_Taxon, "_Fold_", thisFold, "_", thisRegLabel,".Rd"))

                bkg_vals <- predict(maxnet_model, SWD_bkg[, -c(1:3)], clamp = FALSE, type = "cloglog")[, 1]

                occ_vals <- predict(maxnet_model, SWD_occ[SWD_folds[[thisFold]], -c(1:3)], clamp = FALSE, type = "cloglog")[, 1]
                occ_vals_exp <- predict(maxnet_model, SWD_occ[SWD_folds[[thisFold]], -c(1:3)], clamp = FALSE, type = "exponential")[, 1]

                ans <- hmeasure::HMeasure(true.class = c(rep(1, length(occ_vals)), rep(0, length(bkg_vals))), scores = c(occ_vals, bkg_vals))

                k <- length(maxnet_model$betas)
                AIC <- 2*(k - log(prod(occ_vals_exp)))
                AIC_c <- AIC + 2*k*(k+1)/(length(occ_vals_exp))
                cBoyce <- enmSdm::contBoyce(occ_vals, bkg_vals)
                extraMeasures_test <- rbind(extraMeasures_test,
                                            unlist(c(thin = thisThin,
                                                     fold = thisFold,
                                                     reg = thisReg,
                                                     ans$metrics[1, ],
                                                     cBoyce = cBoyce,
                                                     AIC = AIC,
                                                     AIC_c = AIC_c)))

                occ_vals <- predict(maxnet_model, SWD_occ[-SWD_folds[[thisFold]], -c(1:3)], clamp = FALSE, type = "cloglog")[, 1]
                occ_vals_exp <- predict(maxnet_model, SWD_occ[-SWD_folds[[thisFold]], -c(1:3)], clamp = FALSE, type = "exponential")[, 1]
                ans <- hmeasure::HMeasure(true.class = c(rep(1, length(occ_vals)), rep(0, length(bkg_vals))), scores = c(occ_vals, bkg_vals))
                AIC <- 2*(k - log(prod(occ_vals_exp)))
                AIC_c <- AIC + 2*k*(k+1)/(length(occ_vals_exp))
                cBoyce <- enmSdm::contBoyce(occ_vals, bkg_vals)

                extraMeasures_train <- rbind(extraMeasures_train,
                                             unlist(c(thin = thisThin,
                                                      fold = thisFold,
                                                      reg = thisReg,
                                                      ans$metrics[1, ],
                                                      cBoyce = cBoyce,
                                                      AIC = AIC,
                                                      AIC_c = AIC_c)))
              }
            }
            #cat("\n")
          }
        }

        # Save tables...
        write.csv(extraMeasures_test, paste0(taxonPath, aggFactor , "/modelPerformance_test_extraMeasures_", this_Taxon, ".csv"), row.names = FALSE)
        write.csv(extraMeasures_train, paste0(taxonPath, aggFactor , "/modelPerformance_train_extraMeasures_", this_Taxon, ".csv"), row.names = FALSE)
      }
      else
        cat("    Failed fit...skipping this agg level\n")
    }
  }
  cat("*** End of processing\n")
}
