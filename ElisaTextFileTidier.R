# This is an R script that will take in an ELISA test file and turn it into a tidy, tabular file.
# This requires you to know which well has which condition, and you will have to label them based on a strict set of labels
# If a strict set of labels is not used, further analyses become very difficult

# TODO:
# Fallback if date is not used as file name
# Standard Data workup?
# Deal with "Range?"
# Database loading?
# Change static file system to dynamic
# Spit standard data out as is (don't forget naming conventions, ie DATE_Standards)
# Ask for date, passage number, and cell line
# Do calculations
# Name "[Date] ELISA file (tidy)"
# Spit out file
# Names to allcaps/allmin
# Rename these god awful variable names

# Assumptions:
# Assumes the numbers at the beginning of the file name are the date in the form of MMDDYY
# Assumes that 20 mL was used, and that 5 mL was supposed to be used (simple query can fix this)


insulinReader <- function(fileName){
  {
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
  elisaData <- allYouNeed[(endLocations3[1]+2):length(allYouNeed[,1]),] %>%
    fill(c(X1,X5,X6,X7,X8,X9)) %>%
    as.tibble()
  colnames(elisaData) <- as.character(unlist(elisaData[1,]))
  elisaData <- elisaData[-1,]
  } # File reading and wrangling code
  
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
    conditionNames <- paste("Condition ", rep(1:length(everyotherUnique), newLengthBetween))
  }
  # If not done in 'regular' duplicate, you have to enter in which unknowns are associated with which others.
  else{
      # Takes a specially formatted input from the user (whitespace independent). Don't add "Un"! Just the numbers.
      irregularConditionsList <- function(){
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
      splitByCondition <- str_replace_all(irregularConditionsList(), " ","")
      # Splits (first) by semicolon
      splitByCondition <- strsplit(splitByCondition, ";")
      splitByCondition <- matrix(unlist(splitByCondition)) #lapply instead?
      # Then, splits by comma (Now a 2D array)
      splitByUnknown <- apply(splitByCondition, MARGIN = 1, function(x) unlist(strsplit(x, ",")))
      # Checks to see if Un number is below 10, adds a leading 0 if it needs it.
      splitByUnknown <- lapply(splitByUnknown, function(x){
        sapply(x, function(y){
          if((as.integer(x) < 10) && (substring(x, 1, 1)!="0")){x <- paste0("0", x)}
        })
        z <- paste0("Un",x)
      })

      # Finds the rows in the original file with that Un associated with it, adds "Condition N" to each one, in order.
      # Needs to be able to accept names later!
      # List the condition names, sequentially from Condition 1.
      # Look to your records to make sure you are naming the same conditions the same thing. 
      # Because of this, make sure your naming is DESCRIPTIVE and CONSISTENT.
      # Make sure this thing actually is caps insensitive
      
      # Query user for condition names
      userConditionNames <- function(){
        usrResponse <- readline(prompt="Please enter in condition names, from condition 1 to condition 'n'. Caps insensitive. Separate by semicolons: ")
        return(as.character(usrResponse))
      }
      
      givenConditionNames <- userConditionNames()
      givenConditionNames <- unlist(strsplit(givenConditionNames, ";"))
      givenConditionNames <- trimws(givenConditionNames)
      
      for (i in 1:length(splitByUnknown)){
        associatedRows <- grep(paste(splitByUnknown[[i]], collapse ="|"), elisaData$Sample)
        elisaData[associatedRows,10] <- paste(givenConditionNames[i])
      }
      
    }
  
  userConcentration <- function(){
    usrResponse <- readline(prompt="To account for dilution, divide by...(if no dilution, enter 1): ")
    return(as.numeric(as.character(usrResponse)))
  }
  
  userConcentrationValue <- userConcentration()
  
  groupCond <- as.tibble(elisaData) %>%
    group_by(V10) %>%
    mutate("numericMeanConc" = as.numeric(as.character(Mean_Conc))) %>%
    mutate("conditionMean" = mean(numericMeanConc)) %>%
    mutate("adjCondMean" = mean(conditionMean/userConcentrationValue))

  groupCond$everyOtherIndex <- rep_len(c(1,0), nrow(groupCond))
  View(groupCond)
  obsPerCondition <- count(groupCond, V10)
  obsPerCondition <- obsPerCondition[,-1]
  sdRows <- list()
  for (i in 1:nrow(obsPerCondition)) {
    sdRows[[i]] <- seq(1, obsPerCondition[[i,1]], by = 2)
  }
  print(sdRows)
  groupCond2 <- group_by(groupCond, V10) %>%
    filter(everyOtherIndex == 1) %>%
    mutate("condSD" = sd(numericMeanConc))

  groupCond$condSD <- rep(groupCond2$condSD, each = 2)
  View(groupCond)
  # I'll still need the SEM and % above glucose
  # I'll likely want to do those operations on a table that only contains every other datapoint - or rather, maybe just one with the condition names and means.
  # Call it a 'reduced' tibble or something. Save the other data though.

}
