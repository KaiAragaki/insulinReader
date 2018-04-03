# This is an R script that will take in an ELISA test file and turn it into a tidy, tabular file.
# This requires you to know which well has which condition, and you will have to label them based on a strict set of labels
# If a strict set of labels is not used, further analyses become very difficult

# NOTES:
# Input file must be straight off the plate reader
# Filename cannot include .txt - just filename
# High Glucose condition will be used for GSIS, and must be named "High Glucose" or some capitalization derivative thereof
# Assumes the numbers at the beginning of the file name are the date in the form of MMDDYY
# Assumes files are one level below the script, in a folder called "Elisa Text Files"

# TODO:
# Documentation
# Do something with Passage Number and Cell Line Data
# What if the date ISN'T in the filename?

# MAYBE:
# Invalid response catch for Y/N
# More invalid entry catching
# Catch "looks like there were some unnamed conditions" check for NAs
# Make a "easy Excel" file
# Default values
# Setwd back as soon as possible in case of errors

insulinReader <- function(fileName){
  {
  # Attach required packages
  library(tidyverse)
     
  path <- file.path(".", "Elisa Text Files") %>%
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
  elisaData <- allYouNeed[(endLocations3[1]+2):length(allYouNeed[,1]),] %>%
    fill(c(X1,X5,X6,X7,X8,X9)) %>%
    as.tibble()
  colnames(elisaData) <- as.character(unlist(elisaData[1,]))
  elisaData <- elisaData[-1,]
  elisaData[elisaData == "Range?"] <- NA
  } # File reading and wrangling code
  {
  # Checks for number and position of Un
  uniqueUnknowns <- unique(elisaData$Sample)
  uniqueUnknownsNum <- length(uniqueUnknowns)
  unknownLocations <- match(uniqueUnknowns, elisaData$Sample)
  
  # Finds the length between Un(N) and Un(N+1) - Or rather, if the Un was done in duplicate, triplicate, etc.
  samplesPerUnknown <- diff(unknownLocations)

  # Tacks on the last result length in a kind of awkward, but working, way.
  samplesPerUnknown <- append(samplesPerUnknown, (nrow(elisaData)- tail(unknownLocations,1) + 1))

  # Checks to see if you did duplicates (within the assay, not the ELISA)
  remainderUnknownNumber <- uniqueUnknownsNum%%2
  
  # Query user if duplicate check passes
  normalConditionCheck <- function(){
    usrResponse <- readline(prompt="Is Un01 a duplicate of Un02, Un03 of Un04, etc? (Y/N): ")
    return(toupper(toString(usrResponse)))
  }
  
  # Checks to see if there's an even amount of Un and if you did it the normal way
  if (remainderUnknownNumber == 0 && normalConditionCheck() == "Y"){
    # conditionNumberSeq <- c(1:uniqueUnknownsNum)
    # conditionNumberSeq <- str_pad(conditionNumberSeq, 2, "left", pad = "0")
    # conditionNumberSeq <- paste0("Un", conditionNumberSeq)
    # splitByUnknown <- matrix(conditionNumberSeq, ncol=2, nrow = uniqueUnknownsNum/2, byrow = TRUE)
    splitByUnknown <- c()
    for (i in 1:(uniqueUnknownsNum/2)){
      splitByUnknown[[i]] <- c(paste0("Un", str_pad((i*2 - 1), 2, "left", pad = "0")),paste0("Un", str_pad((i*2), 2, "left", pad = "0")))
    }
  }
  # If not done in 'regular' duplicate, you have to enter in which unknowns are associated with which others.
  else{
      # Takes a specially formatted input from the user (whitespace independent). Don't add "Un"! Just the numbers.
      conditionsList <- function(){
        writeLines(
          " List Unk# Conditions manually: 
        \n Same condition: Separate by comma (eg 1,2) 
        \n Different condition: Separate by semicolon (eg 1,2; 3,4)\n"
          )
        usrResponse <- readline(prompt= " Condition list: ")
        return(usrResponse)
      }
      splitByUnknown <- c()
      # Removes whitespace
      splitByCondition <- str_replace_all(conditionsList(), " ","")
      # Splits (first) by semicolon
      splitByCondition <- strsplit(splitByCondition, ";")
      splitByCondition <- matrix(unlist(splitByCondition))
      # Then, splits by comma (Now a 2D array)
      splitByUnknown <- lapply(splitByCondition, function(x) unlist(strsplit(x, ",")))
      # Checks to see if Un number is below 10, adds a leading 0 if it needs it.
      splitByUnknown <- lapply(splitByUnknown, function(x){
        x <- str_pad(x, 2, "left", pad = "0")
        x <- paste0("Un", x)
      })
  }
  print(splitByUnknown)
  # Finds the rows in the original file with that Un associated with it, adds "Condition N" to each one, in order.
  # List the condition names, sequentially from Condition 1.
  # Look to your records to make sure you are naming the same conditions the same thing. 
  # Because of this, make sure your naming is DESCRIPTIVE and CONSISTENT.
  # Make sure this thing actually is caps insensitive
  # Query user for condition names
  userConditionNames <- function(){
    usrResponse <- readline(prompt="Enter in condition names, from Condition 1 to Condition 'n'. Caps insensitive. Separate by semicolons: ")
    return((tolower(as.character(usrResponse))))
  }
  
  givenConditionNames <- userConditionNames()
  givenConditionNames <- unlist(strsplit(givenConditionNames, ";"))
  givenConditionNames <- trimws(givenConditionNames)
  
  # Finds the rows each condition is associated with
  # Look, I know apply functions are faster. But this works.
  for (i in 1:length(splitByUnknown)){
    associatedRows <- grep(paste(splitByUnknown[[i]], collapse ="|"), elisaData$Sample)
    elisaData[associatedRows,10] <- paste(givenConditionNames[i])
  }
  
  # Asks for passage number of cells
  passageNumber <- function(){
    usrResponse <- readline(prompt="What passage number are these cells? (If unknown, enter 'NA'): ")
    return(toupper(usrResponse))
  }
  passageNumber <- passageNumber()
  
  # Asks for cell line of cells. Something else that needs controlled vocabulary.
  cellLine <- function(){
    usrResponse <- readline(prompt="What cell line are these cells?: ")
    return(toupper(usrResponse))
  }
  cellLine <- cellLine()
  
  # Checks to see if it got the date right - if not, it asks the user for the date
  checkFileDate <- function(){
    usrResponse <- readline(prompt=paste("Is", fileDate, "the correct date for this file? (Y/N): "))
    return(toupper(toString(usrResponse)))
  }
  
  if (checkFileDate() == "N") {
    newFileDate <- function(){
      usrResponse <- readline(prompt=paste("Please write the correct file date (in the form MMDDYY): "))
      return(toupper(toString(usrResponse)))
    }
    fileDate <- newFileDate()
  }
  
  # Query user about dilution factor
  userConcentration <- function(){
    usrResponse <- readline(prompt="To account for dilution, divide by...(if no dilution, enter 1): ")
    return(as.numeric(as.character(usrResponse)))
  }
  userConcentrationValue <- userConcentration()
  
  } # Gathering "Metadata" through user querying

  groupCond <- as.tibble(elisaData) %>%
    group_by(V10) %>%
    mutate("numericMeanConc" = as.numeric(as.character(Mean_Conc))) %>%
    mutate("conditionMean" = mean(numericMeanConc, na.rm = TRUE)) %>%
    mutate("adjCondMean" = mean(conditionMean/userConcentrationValue, na.rm = TRUE))

  groupCond$everyOtherIndex <- rep_len(c(1,0), nrow(groupCond))
  obsPerCondition <- count(groupCond, V10)
  obsPerCondition <- obsPerCondition[,-1]
  sdRows <- list()
  for (i in 1:nrow(obsPerCondition)) {
    sdRows[[i]] <- seq(1, obsPerCondition[[i,1]], by = 2)
  }
  groupCond2 <- group_by(groupCond, V10) %>%
    filter(everyOtherIndex == 1) %>%
    mutate("condSD" = sd(numericMeanConc, na.rm = TRUE)) %>%
    mutate("condSEM" = condSD/sqrt(n())) %>%
    ungroup() %>%
    mutate("percentAboveGSIS" = (((conditionMean - first(conditionMean[V10 == "high glucose"]))/first(conditionMean[V10 == "high glucose"]) * 100)))

  groupCond$condSD <- rep(groupCond2$condSD, each = 2)
  groupCond$condSEM <- rep(groupCond2$condSEM, each = 2)
  groupCond$percentAboveGSIS <- rep(groupCond2$percentAboveGSIS, each = 2)
  
  # Drops unneeded columns
  groupCond <- groupCond[,-c(8,11,14)]
  View(groupCond)
  # Exports the files
  write.table(groupCond, file = paste(fileDate,"Tidied"), sep = "\t", col.names = NA)
  write.table(standardData, file = paste(fileDate, "Standard Data Tidied"), sep = "\t", col.names = NA)
  
  # Brings working directory back to before
  setwd("..")
}
