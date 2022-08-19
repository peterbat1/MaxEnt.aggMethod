# Collate table of features and coefficients
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
# 2021-03-10
# 2022-08-16: Tweak code to make it more generalised

# library(fitMaxnet)
# library(dplyr)
# library(stringr)
# library(Rfast)

#' ENM_features_vars_used
#'
#' @param numThis
#' @param numFolds
#' @param baseFolder
#' @param taxonTable_filename
#'
#' @return
#' @export
#'
#' @examples
ENM_features_vars_used <- function(numThis = 5, numFolds = 5,baseFolder = "", taxonTable_filename = "")
{
  ###############################################################################
  # Set critical parameters to drive the process

  # Number of cross-validation folds
  #numFolds <- 5

  # Number of replicate thinning runs
  #numThins <- 5
  ###############################################################################

  #baseFolder <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"

  taxonTable <- read.csv(taxonTable_filename, stringsAsFactors = FALSE)
  rownames(taxonTable) <- taxonTable$taxon
  taxonList <- taxonTable$taxon

  cat("MaxEnt ENMs fitted using Agregation Method: Tabulate features used and prevalence in replicate models\n=====================================================================================================\n")

  for (thisTaxon in taxonList)
  {
    cat(thisTaxon, ":")

    if (taxonTable[thisTaxon, "bestReg"] > 0)
    {
      this_Taxon <- gsub(" ", "_", thisTaxon, fixed = TRUE)

      bestReg <- taxonTable[thisTaxon, "bestReg"]

      aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "best_AggLevel"])

      taxonFolder <- paste0(baseFolder, thisTaxon, "/", aggStr, "/")

      resultsTable <- NULL
      theModelFiles <- paste0(paste0(taxonFolder, "thin_", 1:numThins), "/Fold_", rep(1:numFolds, each = numThins), "/reg_", bestReg, "/", this_Taxon,  "_Fold_", rep(1:numFolds, each = numThins), "_reg_", bestReg, ".Rd")

      for (thisModelFile in theModelFiles)
      {
        #cat("Fold ", thisFold, "\n")
        load(thisModelFile)

        varsUsed <- sort(names(maxnet_model$varmax))

        featureBetas <- maxnet_model$betas

        tmpNames <- gsub("I(", "", names(featureBetas), fixed = TRUE)
        ii <- grep("^2)", tmpNames, fixed = TRUE)
        tmpNames[ii] <- gsub("^2)", "", tmpNames[ii], fixed = TRUE)
        tmpNames[ii] <- paste0(tmpNames[ii],":", tmpNames[ii])
        names(featureBetas) <- tmpNames

        betaValues <- data.frame(betaNames = tmpNames, beta = featureBetas)

        featuresUsed <- gsub("I(", "", sort(names(maxnet_model$featuremaxs)), fixed = TRUE)
        ii <- grep("^2)", featuresUsed, fixed = TRUE)
        featuresUsed[ii] <- gsub("^2)", "", featuresUsed[ii], fixed = TRUE)
        featuresUsed[ii] <- paste0(featuresUsed[ii],":",featuresUsed[ii])

        featuresUsed_df <- data.frame(betaNames = featuresUsed)

        ans <- dplyr::left_join(featuresUsed_df, betaValues, "betaNames")
        ans <- ans[order(ans$betaNames), ]

        if (is.null(resultsTable))
        {
          resultsTable <- data.frame(ans[, 2])
          rownames(resultsTable) <- ans[, 1]
        }
        else
          resultsTable <- cbind(resultsTable,  ans[, 2])
      }

      #colnames(resultsTable) <- paste0(this_Taxon, "_", 1:length(theModelFiles))
      colnames(resultsTable) <- paste0("R_", stringr::str_pad(as.character(1:length(theModelFiles)),
                                                              side = "left",
                                                              width = 2,
                                                              pad = "0"))

      resultsTable_zero <- resultsTable
      resultsTable_zero[which(is.na(resultsTable_zero), arr.ind = TRUE)] <- 0

      mandm <- t(Rfast::rowMinsMaxs(as.matrix(resultsTable_zero)))
      rownames(mandm) <- rownames(resultsTable)
      minVals <- apply(resultsTable_zero, 1, min)

      resultsTable_full <- data.frame(resultsTable,
                                      votes = apply(resultsTable, 1, function(x) {sum(!is.na(x))}))

      resultsTable_summary <- data.frame(mandm,
                                         votes = apply(resultsTable, 1, function(x) {sum(!is.na(x))}))

      write.csv(resultsTable_full, paste0(taxonFolder, this_Taxon, "_feature_table.csv"))
      write.csv(resultsTable_summary, paste0(taxonFolder, this_Taxon, "_feature_summary_table.csv"))
      cat(" done\n")
    }
    else
      cat(" No model found...skipping\n")
  }
}
