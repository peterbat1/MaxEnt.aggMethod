# Collate information to create model multi-run response plots
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
# 2021-03-22
# 2022-08-16: Code clean-up

# library(maxnet)
# library(ggplot2)
# library(ggpubr)

#' ENM_response_plots
#'
#' @param numThins Integer
#' @param numFolds Integer
#' @param responseType Character
#' @param taxonTable_filename Character
#' @param envFolder Character
#' @param basePath Character
#'
#' @return
#' @export
#'
#' @examples
ENM_response_plots <- function(numThins = 5,
                               numFolds = 5,
                               responseType = "cloglog",
                               taxonTable_filename = "",
                               envFolder = "",
                               basePath = "")
{
###############################################################################
# Set critical parameters to drive the process

# Number of cross-validation folds
#numFolds <- 5

# Number of replicate thinning runs
#numThins <- 5

#responseType <- "cloglog"
###############################################################################

taxonTable <- read.csv(taxonTable_filename, stringsAsFactors = FALSE)
rownames(taxonTable) <- taxonTable$taxon
taxonList <- taxonTable$taxon

#envFolder <- "/home/peterw/Myotis/narclim/nsw_cropped"

#basePath <- "/home/peterw/Myotis/PLP/MaxEnt_models/"

allVarNames <- gsub(".tif", "", list.files(envFolder, "*.tif"), fixed = TRUE)

######### Change this line to list the names of covariates which are to be
######### excluded; comment the line out entirely if there are no excluded env
######### covariates
allVarNames <- allVarNames[-match(c("bioclim_08_albers", "bioclim_09_albers"), allVarNames)]

varImp <- read.csv("/home/peterw/Nyctimene/PLP/Metadata/MEAN_variable_importance_summary_refit.csv", stringsAsFactors = FALSE)
rownames(varImp) <- varImp$taxon

#thisTaxon <- "Petauroides volans"
for (thisTaxon in taxonList)
{
  cat(thisTaxon, "\n")
  this_Taxon <- gsub(" ", "_", thisTaxon, fixed = TRUE)

  bestReg <- taxonTable[thisTaxon, "bestReg"]

  aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "bestAggLevel"])

  taxonResults <- vector("list", length(allVarNames))
  names(taxonResults) <- allVarNames

  if (bestReg > 0)
  {
    for (thisThin in 1:numThins)
    {
      for (thisFold in 1:numFolds)
      {
        load(paste0(basePath, thisTaxon, "/", aggStr, "/thin_", thisThin, "/Fold_", thisFold, "/reg_", bestReg, "/", this_Taxon, "_Fold_", thisFold, "_reg_", bestReg, ".Rd"))

        modelVarNames <- names(maxnet_model$samplemeans)

        for (thisVar in modelVarNames)
        {
          # For the current variable, make a mean value matrix and then replace the
          # variable column with a sequence of values spanning the range of values for
          # the variable encountered during model fitting
          hasLevels <- !is.null(unlist(maxnet_model$levels[thisVar]))

          if (hasLevels)
            numRows <- unlist(maxnet_model$levels[thisVar])
          else
            numRows <- 100

          d <- data.frame(matrix(unlist(maxnet_model$samplemeans), numRows, length(maxnet_model$samplemeans), byrow = TRUE))
          colnames(d) <- names(maxnet_model$samplemeans)

          varMax <- maxnet_model$varmax[thisVar]
          varMin <- maxnet_model$varmin[thisVar]

          if (hasLevels) d[, thisVar] <- unlist(maxnet_model$levels[thisVar])     else
            d[, thisVar] <- seq(varMin - 0.1 * (varMax - varMin), varMax + 0.1 * (varMax - varMin), length = 100)

          # if (is.null(ylabel))
          #   ylabel <- paste0("Response (", responseType, ")")
          predVals <- predict(maxnet_model, d, type = responseType)

          modelResults <- data.frame(prediction = predVals[,1], d[, thisVar])
          colnames(modelResults) <- c("prediction", thisVar)
          #taxonResults[[thisVar]] <- rbind(taxonResults[[thisVar]], data.frame(prediction = predVals[,1], d[, thisVar]))
          taxonResults[[thisVar]] <- rbind(taxonResults[[thisVar]], modelResults)
        }
      }
    }

    #### saveRDS
    saveRDS(taxonResults, file = paste0(basePath, thisTaxon, "/", aggStr, "/", this_Taxon, "_response_plot_data.rds"))

    #### Make plots
    plotyBits <- vector("list", length(allVarNames))
    names(plotyBits) <- allVarNames
    #thisVar <- allVarNames[3]
    for (thisVar in allVarNames)
    {
      thisVarImp <- varImp[thisTaxon, ]

      plotData <- taxonResults[[thisVar]]

      plotyBits[[thisVar]] <- ggplot2::ggplot(plotData, aes(x = plotData[, thisVar], y = prediction)) +
        geom_smooth(colour = "blue", size = 1) + ylab("prediction") + xlab(thisVar) + ylim(c(0, 1)) +
        theme(axis.title.x = element_text(size = 8),
              axis.title.y = element_text(size = 8),
              axis.text.x = element_text(size = 6),
              axis.text.y = element_text(size = 6),
              plot.title = element_text(size = 8)) +
        ggtitle(paste0("Importance: ", round(thisVarImp[thisVar], 2), "%"))
    }

    ggpubr::ggarrange(plotlist = plotyBits, ncol = 4, nrow = 3) %>%
      ggpubr::ggexport(filename = paste0(basePath, thisTaxon, "/", aggStr, "/", this_Taxon, "_response_plots.png"))


  } else cat("    *** No models found...skipped\n")
}

}
