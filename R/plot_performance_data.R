# Plot model performance data for MaxEnt ENMs fitted using aggregation method
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
# 2020-10-16; 2021-03-14: Code tidy-up and generalisation
# 2022-08-16: Further clean-up of code

# library(ggplot2)
# library(ggpubr)


#' ENM_plot_model_performance
#'
#' @param aggregation_levels
#' @param regularisation_sequence
#' @param basePath
#' @param taxonTable_filename
#'
#' @return
#' @export
#'
#' @examples
ENM_plot_model_performance <- function(aggregation_levels = c(2, 4, 8, 16),
                                       regularisation_sequence = 1:10,
                                       basePath = "",
                                       taxonTable_filename = "")
{
  #regSequence <- 1:10

  taxonTable <- read.csv(taxonTable_filename, stringsAsFactors = FALSE)
  rownames(taxonTable) <- taxonTable$taxon
  taxonList <- taxonTable$taxon

  basePath <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"

  cat("MaxEnt ENMs fitted using Agregation Method: Produce model performance plots\n===========================================================================\n")

  # Sequence of aggregation levels or factors
  #aggregation_levels <- c(2, 4, 8, 16)
  aggFactorSet <- paste0("aggFactor_", aggregation_levels)

  plotMapping <- matrix(1:(2*length(aggregation_levels)), length(aggregation_levels), 2, byrow = TRUE)
  colnames(plotMapping) <- c("cBoyce", "AUC")
  rownames(plotMapping) <- aggFactorSet

  for (thisTaxon in taxonList)
  {
    cat("  ", thisTaxon, "...", sep = "")
    this_Taxon <- gsub(" ", "_", thisTaxon)

    taxonPath <- paste0(basePath, thisTaxon, "/")

    plotyBits <- vector("list", 2*length(aggregation_levels))

    for (aggFactor in aggFactorSet)
    {
      if (dir.exists(paste0(taxonPath, aggFactor, "/thin_1")))
      {
        d <- read.csv(paste0(taxonPath, aggFactor, "/modelPerformance_test_extraMeasures_", this_Taxon, ".csv"), stringsAsFactors = FALSE)

        #### Continuous Boyce
        plotyBits[[plotMapping[aggFactor, "cBoyce"]]] <- ggplot(d, aes(x = reg, y = cBoyce)) +
          geom_point(colour = d$thin) +
          geom_smooth() +
          ylim(0.5, 1) +
          xlab("Regularisation") +
          ggtitle(substitute(paste(bold(aggStr), "\t", bolditalic(spp)), list(spp = thisTaxon, aggStr = aggFactor))) +
          theme(plot.title = element_text(size = 9))

        #### AUC
        plotyBits[[plotMapping[aggFactor, "AUC"]]] <- ggplot(d, aes(x = reg, y = AUC)) +
          geom_point(colour = d$thin) +
          geom_smooth() +
          ylim(0.5, 1) +
          xlab("Regularisation") +
          ggtitle(substitute(paste(bold(aggStr), "\t", bolditalic(spp)), list(spp = thisTaxon, aggStr = aggFactor))) +
          theme(plot.title = element_text(size = 9))
      }
      else
      {
        # Make a "dummy" d data.frame which gives coords for a "No data" text
        # annotation thus avoiding the complexities of NULL ggplot objects causing
        # heart failure in ggarrange
        d <- data.frame(reg = 5, cBoyce = 0.5, AUC = 0.5)

        #### Continuous Boyce
        plotyBits[[plotMapping[aggFactor, "cBoyce"]]] <- ggplot(d, aes(x = reg, y = cBoyce)) +
          geom_text(x = 5, y = 0.5, label = "No data") +
          ylim(0.5, 1) +
          xlab("Regularisation") +
          ggtitle(substitute(paste(bold(aggStr), "\t", bolditalic(spp)), list(spp = thisTaxon, aggStr = aggFactor))) +
          theme(plot.title = element_text(size = 9))

        #### AUC
        plotyBits[[plotMapping[aggFactor, "AUC"]]] <- ggplot(d, aes(x = reg, y = AUC)) +
          geom_text(x = 5, y = 0.5, label = "No data") +
          ylim(0.5, 1) +
          xlab("Regularisation") +
          ggtitle(substitute(paste(bold(aggStr), "\t", bolditalic(spp)), list(spp = thisTaxon, aggStr = aggFactor))) +
          theme(plot.title = element_text(size = 9))
      }

      #### AIC
      # p <- ggplot(d, aes(x = reg, y = AIC)) +
      #   geom_point(colour = d$thin) +
      #   geom_smooth() +
      #   ggtitle(thisTaxon) +
      #   theme(plot.title=element_text(size = 16, face = "bold.italic"))
      #
      # png(paste0(basePath, thisTaxon, "/", aggFactor, "/modelPerformance_test_extraMeasures_", this_Taxon, "_AIC.png"),
      #     width = 1024, height = 768)
      # plot(p)
      # dev.off()
      #
      # #### TPR/Recall
      # p <- ggplot(d, aes(x = reg, y = TPR)) +
      #   geom_point(colour = d$thin) +
      #   geom_smooth() +
      #   ylim(0, 1) +
      #   ggtitle(thisTaxon) +
      #   theme(plot.title = element_text(size = 16, face = "bold.italic"))
      #
      # png(paste0(basePath, thisTaxon, "/", aggFactor, "/modelPerformance_test_extraMeasures_", this_Taxon, "_TPR.png"),
      #     width = 1024, height = 768)
      # plot(p)
      # dev.off()
      #
      # #### Youden
      # p <- ggplot(d, aes(x = reg, y = Youden)) +
      #   geom_point(colour = d$thin) +
      #   geom_smooth() +
      #   ylim(0, 1) +
      #   ggtitle(thisTaxon) +
      #   theme(plot.title = element_text(size = 16, face = "bold.italic"))
      #
      # png(paste0(basePath, thisTaxon, "/", aggFactor, "/modelPerformance_test_extraMeasures_", this_Taxon, "_YoudenJ.png"),
      #     width = 1024, height = 768)
      # plot(p)
      # dev.off()

      # #### Continuous Boyce
      # plotyBits[[plotMapping[aggFactor, "cBoyce"]]] <- ggplot(d, aes(x = reg, y = cBoyce)) +
      #   geom_point(colour = d$thin) +
      #   geom_smooth() +
      #   ylim(0, 1) +
      #   ggtitle(thisTaxon) +
      #   theme(plot.title = element_text(size = 16, face = "bold.italic"))
      #
      # # png(paste0(basePath, thisTaxon, "/", aggFactor, "/modelPerformance_test_extraMeasures_", this_Taxon, "_Boyce.png"),
      # #     width = 1024, height = 768)
      # # plot(p)
      # # dev.off()
      #
      # #### AUC
      # plotyBits[[plotMapping[aggFactor, "AUC"]]] <- ggplot(d, aes(x = reg, y = AUC)) +
      #   geom_point(colour = d$thin) +
      #   geom_smooth() +
      #   ylim(0, 1) +
      #   ggtitle(thisTaxon) +
      #   theme(plot.title = element_text(size = 16, face = "bold.italic"))

      # png(paste0(basePath, thisTaxon, "/", aggFactor, "/modelPerformance_test_extraMeasures_", this_Taxon, "_AUC.png"),
      #     width = 1024, height = 768)
      # plot(p)
      # dev.off()

      # #### H
      # p <- ggplot(d, aes(x = reg, y = H)) +
      #   geom_point(colour = d$thin) +
      #   geom_smooth() +
      #   ylim(0, 1) +
      #   ggtitle(thisTaxon) +
      #   theme(plot.title = element_text(size = 16, face = "bold.italic"))
      #
      # png(paste0(basePath, thisTaxon, "/", aggFactor, "/modelPerformance_test_extraMeasures_", this_Taxon, "_H.png"),
      #     width = 1024, height = 768)
      # plot(p)
      # dev.off()


    }

    cat("done\n")
    outFilename <- paste0(taxonPath, thisTaxon,  "_performance_plots.png")

    aggLabels <- rep(aggFactorSet, each = 2)

    dummy <- suppressMessages(ggpubr::ggarrange(plotlist = plotyBits, ncol = 2, nrow = 3) %>% #, labels = aggLabels, hjust = -3.5, font.label = list(size = 8))
                                ggpubr::ggexport(filename = outFilename,
                                                 width = 2100,
                                                 height = 2400,
                                                 res = 300,
                                                 verbose = FALSE))
  }
  cat("*** End of processing\n")
}
