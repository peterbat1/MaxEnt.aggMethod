# Calculate summary performance statistics
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
# 2021-03-23 
# 2022-08-16: Code clean-up

library(Rfast)

baseFolder <- "/home/peterw/Downloads/workflow_tests/ENM_Results/"

taskList <- c("test", "train")

taxonTable <- read.csv("/home/peterw/Downloads/workflow_tests/taxonTable_demonstration.csv", stringsAsFactors = FALSE)
rownames(taxonTable) <- taxonTable$taxon
taxonList <- taxonTable$taxon

for (thisTaxon in taxonList)
{
  cat(thisTaxon, ": ")
  
  this_Taxon <- gsub(" ", "_", thisTaxon, fixed = TRUE)
  
  bestReg <- taxonTable[thisTaxon, "bestReg"]
  
  aggStr <- paste0("aggFactor_", taxonTable[thisTaxon, "best_AggLevel"])
  
  if (bestReg > 0)
  {
    if ("test" %in% taskList)
    {
      # Test scores
      perfTable <- read.csv(paste0(baseFolder, thisTaxon, "/", aggStr, "/modelPerformance_test_extraMeasures_", this_Taxon, ".csv"), stringsAsFactors = FALSE)
      
      bestPerf <- perfTable[perfTable$reg == bestReg, ]
      
      meanPerf <- colmeans(as.matrix(bestPerf))
      names(meanPerf) <- colnames(bestPerf)
      
      sdPerf <- colVars(as.matrix(bestPerf), std = TRUE)
      names(sdPerf) <- colnames(bestPerf)
      
      ans <- data.frame(mean = meanPerf, stdDev = sdPerf)
      write.csv(ans, paste0(baseFolder, thisTaxon, "/", aggStr, "/modelPerformance_test_summary_", this_Taxon, ".csv"))
    }
    
    if ("train" %in% taskList)
    {
      #####################################################################################
      # Training scores
      perfTable <- read.csv(paste0(baseFolder, thisTaxon, "/", aggStr, "/modelPerformance_train_extraMeasures_", this_Taxon, ".csv"), stringsAsFactors = FALSE)
      
      bestPerf <- perfTable[perfTable$reg == bestReg, ]
      
      meanPerf <- colmeans(as.matrix(bestPerf))
      names(meanPerf) <- colnames(bestPerf)
      
      sdPerf <- colVars(as.matrix(bestPerf), std = TRUE)
      names(sdPerf) <- colnames(bestPerf)
      
      ans <- data.frame(mean = meanPerf, stdDev = sdPerf)
      write.csv(ans, paste0(baseFolder, thisTaxon, "/", aggStr, "/modelPerformance_train_summary_", this_Taxon, "_refit_aggMethod.csv"))
    }
    
    cat("done\n")
  }
  else
    cat("No mmodel found...skipping\n")
}

