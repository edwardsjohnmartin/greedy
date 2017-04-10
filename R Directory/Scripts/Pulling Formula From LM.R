# gonna use this file with the purpose of figuring out how to extract a formula data type from a lm data type

printSpacing <- function(){
  print("")
  print("*************************************************************************************************")
  print("")
}

getDataReady <- function(){
  readFile <- "Datasets/fbData.csv"
  varData <- read.csv(file = readFile, header = FALSE, skip = 1)
  varData$V2 = NULL
  varData$V3 = NULL
  varData$V4 = NULL
  varData$V5 = NULL
  varData$V6 = NULL
  varData$V7 = NULL
  names(varData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
  
  printSpacing()
  
  #Create a formula model based on the data
  sampleFormula <- Y ~ 1 + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11
  
  print(paste("print sampleFormula which is type -", typeof(sampleFormula), "-", " with length of ", length(sampleFormula), sep = ""))
  print(sampleFormula)
  
  printSpacing()
  
  # Create a lm model to analyze
  sampleLM <- lm(sampleFormula, data = varData)
  
  print(paste("print sampleLM which is type -", typeof(sampleLM), "-", " with length of ", length(sampleLM), sep = ""))
  print(sampleLM)
  
  printSpacing()
  
  #printEverythingData(sampleFormula)
  pullFormulaFromLM(sampleLM)
}

printEverythingData <- function(f1){
  print(paste("print f1 which is type -", typeof(f1), "-", " with length of ", length(f1), sep = ""))
  print(f1)
  
  printSpacing()
  
  for(i in 1:length(f1)){
    printSpacing()
    
    print(paste("Here is f1[", i, "]", sep = ""))
    print(paste("print f1[i] which is type -", typeof(f1[i]), "-", " with length of ", length(f1[i]), sep = ""))
    print(f1[i])
  }
  
  
  printSpacing()
}

pullFormulaFromLM <- function(model){
  #printEverythingData(model)
  
  printSpacing()
  
  # mod.terms is a list of all the variables in the model passed into drops
  mod.terms <- attr(terms(model), "term.labels")
  print(paste("print mod.terms which is type -", typeof(mod.terms), "-", " with length of ", length(mod.terms), sep = ""))
  print(mod.terms)
  
  printSpacing()
  
  extractedFormula <- as.formula(paste(c(paste((attr(model$terms, "variables")[2]), "~ 1 "), mod.terms), collapse = " + "))
  print(paste("print extractedFormula which is type -", typeof(extractedFormula), "-", " with length of ", length(extractedFormula), sep = ""))
  print(extractedFormula)
  
  printSpacing()
  
  #printEverythingData(extractedFormula)
}

getDataReady()
