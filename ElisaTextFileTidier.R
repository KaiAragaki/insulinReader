# This is an R script that will take in an ELISA test file and turn it into a tidy, tabular file.
# This requires you to know which well has which condition, and you will have to label them based on a strict set of labels
# If a strict set of labels is not used, further analyses become very difficult.

# TODO: 
# Force this thing to stop making so many damn assumptions
# Query the user more
# Fallback if date is not used as file name
# Standard Data workup
# Deal with "Range?"

# Look through data and make a guess at the number of wells used
  # If the number of wells used is even, ask if it they were run in duplicate
    # If yes, ask if they were run vertically, horizontally, or neither
      # If vertically, ask for condition names keeping vertical grouping in mind
      # If horizontally, ask for condition names keeping horizontal grouping in mind
      # If neither, read off each well and ask for a condition name
  # If the number of wells is odd, read off each well and ask for a condition name

# Ask for date, passage number, and cell line

# Do calculations here

# Name "[Date] ELISA file (tidy)"

# Spit out file


insulinReader <- function(fileName){
  # Load necessary packages
  require(tidyverse)
  require(stringr)
  # Assumptions:
  # Assumes the numbers at the beginning of the file name are the date in the form of MMDDYY
  # Assumes duplicate of duplicate (ie SOP, 4 well per cond)
  # Assumes condition layout is similar or identical
  
  # Reads the file, skipping the (uneeded) header, and packs it in a tibble
  path <- file.path("C:", "Users", "Adam", "Desktop", "R for Lab Data", "Elisa Text Files") %>%
    setwd()
  
  # Extract numbers from fileName and use them as date
  fileDate <- regmatches(fileName, gregexpr("[[:digit:]]+", fileName))
  
  # Find the actual file
  fileName <- paste0(fileName, ".txt")
  
  # Reads file
  insFile <- readLines(fileName, encoding = "latin1")
  
  # Denotes locations of the phrase "End"  
  endLocations <- grep("~End", unlist(insFile))
  workingTable <- read.table(fileName, fill = TRUE, skip = endLocations[3], sep = "\t", col.names = 1:9, na.strings = "")
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
  
  # Honestly, it's probably okay just to spit the standardData out as is because it's used so little
  
  # Add condition category:
  # Makes sure there's an 'normal' amount of conditions. If there isn't, it assumes the last one was a singlet.
  # Should probably check, rather than by multiples of 4, if # of unique unknowns is even.
  
  rowResults <- c()
  
  uniqueConditions <- unique(conditionData[,1])
  uniqueConditionsNum <- nrow(uniqueConditions)
  everyotherUnique <- uniqueConditions[seq(1,length(unlist(uniqueConditions)),2),]
  everyotherUnique <- pull(everyotherUnique, 1)
  for (i in 1:length(everyotherUnique)) {
    thing <- min(grep(everyotherUnique[i], pull(conditionData, 1)))
    rowResults[i] <- as.integer(thing)
  }

  # However inefficiently, I now have the row number for every other Unk#
  lengthbetween <-c()
  for (i in 1:(length(rowResults)-1)){
    lengthbetween[i] <- rowResults[i+1] - rowResults[i]
  }
  lastResult <- rowResults[length(rowResults)]
  newLengthBetween <- append(lengthbetween, (nrow(conditionData)-lastResult + 1))

  modulusConditionNumber <- uniqueConditionsNum%/%2
  remainderConditionNumber <- uniqueConditionsNum%%2
  if (remainderConditionNumber == 0){
    normalConditionCheck <- function(){
      usrResponse <- readline(prompt="Is Unknown 1 a duplicate of Unknown 2, etc? (Y/N): ")
      return(toupper(toString(usrResponse)))
    }
    if(normalConditionCheck() == "Y"){
      conditionNames <- paste("Condition ", rep(1:length(everyotherUnique), newLengthBetween))
    }
    else{
      irregularConditionsList <- function(){
        # Make sure capital and lowercase are both accepted!
        writeLines("List Unk# Conditions manually: \n Same condition: Separate by comma (eg 1,2) \n Different condition: Separate by semicolon (eg 1,2; 3,4)")
        usrResponse <- readline(prompt= "Condition list: ")
        
        return(usrResponse)
      }
      testarray <- c()
      newtestarray <- c()
      
      testarray <- irregularConditionsList()
      testarray <- str_replace_all(testarray, " ","")
      testarray <- strsplit(testarray, ";")
      testarray <- unlist(testarray)
      for (i in 1:length(testarray)){
        newtestarray[i] <-strsplit(testarray[i], ",")
      }
      for (i in 1:length(newtestarray)){
        for (j in 1:length(newtestarray[[i]])){
          if (as.integer(newtestarray[[i]][j]) < 10){
            # Entry cannot contain leading 0s. Alternately, just do some regex magic and check for numeric equality.
            newtestarray[[i]][j] <- paste0("0",newtestarray[[i]][j])
          }
          newtestarray[[i]][j] <- paste0("Un",newtestarray[[i]][j])
        }
      }
      print(newtestarray)
      # Pick up here. Loop through each row and add condition number "N" to each Un associated with that row.
    }
  }
  else {
    conditionNames <- paste("Condition", 1:((nrow(conditionData)/4))) %>%
      rep(each = 4)
    conditionNamesRemainder <- paste("Condition", modulusConditionNumber+1) %>%
      rep(each = remainderConditionNumber)
    conditionNames <- c(conditionNames, conditionNamesRemainder)
  }
  
  conditionData[,10] <- conditionNames
  groupCond <- as.tibble(conditionData) %>%
    group_by(V10) %>%
    mutate("numericMeanConc" = as.numeric(as.character(Mean_Conc))) %>%
    mutate("conditionMean" = mean(numericMeanConc)) %>%
    mutate("adjCondMean" = mean(conditionMean/4)) %>%
    mutate("condSD" = sd(c(numericMeanConc[1]/4, numericMeanConc[3]/4))) %>%#Choosing all values would result in overcounting and lower SD - not cool.
    mutate("lower" = adjCondMean-condSD) %>%
    mutate("upper" = adjCondMean+condSD)

  
  View(groupCond)  
  # I'll still need the SEM, % above glucose, and labelling capabilities, but for now I'll try plotting this mess.
  theNumbers <- seq(from = 1, to = length(groupCond$V10), by = 4)
  
  ggplot(groupCond[theNumbers,], mapping = aes(x = V10, y = adjCondMean, fill = V10)) +
    geom_col() +
    geom_errorbar(ymin = groupCond[theNumbers, ]$lower, ymax = groupCond[theNumbers, ]$upper)
  
  # After all that drama, standardData and conditionData now are (semi) normal tables that can be worked with
  
  #  conditionsTable <- workingTable[(endLocations2[1]+2):(endLocations2[2] - 1),]
  
  #insFile <- read.table(insFileName, fill = TRUE, skip = endLocations[3], sep = "\t")

  #insFile <- read.table(insFileName, sep = "\t", fill = TRUE, skip = endLocations[2])
  
  #read.table(sep = "\t", fill = TRUE) %>%
  #as.tibble()
}