insulinReader <- function(dateArg){
  require(tidyverse)
  # Assumptions
  # Assumes duplicate of duplicate (ie SOP, 4 well per cond)
  # Assumes condition layout is similar or identical
  # Reads the file, skipping the (uneeded) header, and packs it in a tibble
  path <- file.path("C:", "Users", "Adam", "Desktop") %>%
    setwd()
  insFileName <- paste0(dateArg, "_rat.txt")
  insFile <- readLines(insFileName, encoding = "latin1")
  endLocations <- grep("~End", unlist(insFile))
  workingTable <- read.table(insFileName, fill = TRUE, skip = endLocations[3], sep = "\t", col.names = 1:9, na.strings = "")
  endLocations2 <- grep("~End", unlist(workingTable))
  allYouNeed <- workingTable[2:(endLocations2[2]-1),]
  endLocations3 <- grep("~End", unlist(allYouNeed))
  standardData <- allYouNeed[1:(endLocations3[1]-1),]
  standardData <- standardData[,-ncol(standardData)]
  conditionData <- allYouNeed[(endLocations3[1]+2):length(allYouNeed[,1]),] %>%
    fill(c(X1,X5,X6,X7,X8,X9)) %>%
    as.tibble()
  colnames(conditionData) <- as.character(unlist(conditionData[1,]))
  conditionData <- conditionData[-1,]
  
  # ConditionData done. Standard Data needs some work
  
  # Add condition category:
  
  conditionNames <- paste("Condition", 1:(nrow(conditionData)/4)) %>%
    rep(each = 4)
  conditionData[,10] <- conditionNames
  
  groupCond <- as.tibble(conditionData) %>%
    group_by(V10) %>%
    mutate("numericMeanConc" = as.numeric(as.character(Mean_Conc))) %>%
    mutate("conditionMean" = mean(numericMeanConc)) %>%
    mutate("adjCondMean" = mean(conditionMean/4)) %>%
    mutate("condSD" = sd(c(numericMeanConc[1]/4, numericMeanConc[3]/4))) #Choosing all values would result in overcounting and lower SD - not cool.
  
  View(groupCond)
  
  # I'll still need the SEM, % above glucose, and labelling capabilities, but for now I'll try plotting this mess.
  
  ggplot(groupCond, mapping = aes(x = V10, y = adjCondMean, fill = V10)) +
    geom_col()
  
  # After all that drama, standardData and conditionData now are (semi) normal tables that can be worked with

#  conditionsTable <- workingTable[(endLocations2[1]+2):(endLocations2[2] - 1),]

  #insFile <- read.table(insFileName, fill = TRUE, skip = endLocations[3], sep = "\t")
  
  #print(endLocations)
  #insFile <- read.table(insFileName, sep = "\t", fill = TRUE, skip = endLocations[2])
    
    #read.table(sep = "\t", fill = TRUE) %>%
    #as.tibble()
  # This is really all you need
}
