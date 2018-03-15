# This is an R script that will take in an ELISA test file and turn it into a tidy, tabular file.
# This requires you to know which well has which condition, and you will have to label them based on a strict set of labels
# If a strict set of labels is not used, further analyses become very difficult.

# TODO: 
# Force this thing to stop making so many damn assumptions
# Query the user more
# Fallback if date is not used as file name
# Standard Data workup
# Deal with "Range?"
# Database loading?
# Change static file system to dynamic
# Spit standard data out as is (don't forget naming conventions, ie DATE_Standards)
# Ask for date, passage number, and cell line
# Do calculations
# Name "[Date] ELISA file (tidy)"
# Spit out file

# Assumptions:
# Assumes the numbers at the beginning of the file name are the date in the form of MMDDYY

insulinReader <- function(fileName){
  
  # Attach required packages
  library(tidyverse) 
  library(stringr)

  # Static file path. Set as working directory. (Will need to change before deployment)
  path <- file.path("C:", "Users", "Adam", "Desktop", "R for Lab Data", "Elisa Text Files") %>%
    setwd()
  
  # Extract numbers from fileName and use them as date
  fileDate <- regmatches(fileName, gregexpr("[[:digit:]]+", fileName))
  
  # Find file name
  fileName <- paste0(fileName, ".txt")

  # Attempts to find and read file
  insFile <- readLines(fileName, encoding = "latin1")
  
  # Denotes locations of the phrase "End"  
  endLocations <- grep("~End", unlist(insFile))
  workingTable <- read.table(fileName, fill = TRUE, skip = endLocations[3], sep = "\t", col.names = 1:9, na.strings = "")
  endLocations2 <- grep("~End", unlist(workingTable))
  allYouNeed <- workingTable[2:(endLocations2[2]-1),]
  endLocations3 <- grep("~End", unlist(allYouNeed))
  
  # Takes the standard data and puts it in a tibble
  standardData <- allYouNeed[1:(endLocations3[1]-1),]
  standardData <- standardData[,-ncol(standardData)]
  
  # Takes the condition data and puts it in a tibble
  conditionData <- allYouNeed[(endLocations3[1]+2):length(allYouNeed[,1]),] %>%
    fill(c(X1,X5,X6,X7,X8,X9)) %>%
    as.tibble()
  colnames(conditionData) <- as.character(unlist(conditionData[1,]))
  conditionData <- conditionData[-1,]

  # Checks for number and position of Un
  uniqueConditions <- unique(conditionData$Sample)
  uniqueConditionsNum <- length(uniqueConditions)
  conditionLocations <- match(uniqueConditions, conditionData$Sample)
  
  # Finds the length between Un(N) and Un(N+1) - Or rather, if the Un was done in duplicate, triplicate, etc.
  lengthbetween <- diff(conditionLocations)

  # Tacks on the last result length in a kind of awkward, but working, way.
  lengthbetween <- append(lengthbetween, (nrow(conditionData)- tail(conditionLocations,1) + 1))

  # Checks to see if you did duplicates (within the assay, not the ELISA)
  modulusConditionNumber <- uniqueConditionsNum%/%2
  remainderConditionNumber <- uniqueConditionsNum%%2
  
  normalConditionCheck <- function(){
    usrResponse <- readline(prompt="Is Un01 a duplicate of Un02, Un03 of Un04, etc? (Y/N): ")
    return(toupper(toString(usrResponse)))
  }
  
  # Checks to see if there's an even amount of Un and if you did it the normal way
  if (remainderConditionNumber == 0 && normalConditionCheck() == "Y"){
    conditionNames <- paste("Condition ", rep(1:length(everyotherUnique), newLengthBetween))
    # If you did it the normal way, your work is done. (will still need condition naming)
}
  # If you didn't do it the normal way, you have to enter in which unknowns are associated with which others.
  else{
      # Takes a specially formatted input from the user (whitespace independent, due to later steps.) don't add "Un"! Just the numbers
      irregularConditionsList <- function(){
        writeLines(
          " List Unk# Conditions manually: 
        \n Same condition: Separate by comma (eg 1,2) 
        \n Different condition: Separate by semicolon (eg 1,2; 3,4)\n"
          )
        usrResponse <- readline(prompt= " Condition list: ")
        return(usrResponse)
      }
      testarray <- c()
      newtestarray <- c()
      
      
      testarray <- irregularConditionsList()
      # Removes whitespace
      testarray <- str_replace_all(testarray, " ","")
      # Splits (first) by semicolon
      testarray <- strsplit(testarray, ";")
      testarray <- unlist(testarray)
      # Then, splits by colon (Now a 2D array)
      for (i in 1:length(testarray)){
        newtestarray[i] <-strsplit(testarray[i], ",")
      }
      for (i in 1:length(newtestarray)){
        for (j in 1:length(newtestarray[[i]])){
          # Checks to see if Un number is below 10, adds a leading 0 if it needs it.
          if (as.integer(newtestarray[[i]][j]) < 10 && substring(newtestarray[[i]][j], 0, 1)!=0){
            newtestarray[[i]][j] <- paste0("0",newtestarray[[i]][j])
          }
          newtestarray[[i]][j] <- paste0("Un",newtestarray[[i]][j])
        }
      }
      
      # Finds the rows in the original file with that Un associated with it, adds "Condition N" to each one, in order.
      # Needs to be able to accept names later!
      for (i in 1:length(newtestarray)){
        associatedRows <- grep(paste(newtestarray[[i]], collapse ="|"), conditionData$Sample)
        conditionData[associatedRows,10] <- paste("Condition", i)
      }
      
    }
  
  # This next part should be replaced by the manual querying method
  #conditionData[,10] <- conditionNames
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
