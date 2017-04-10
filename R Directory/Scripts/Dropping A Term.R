# gonna use this file with the purpose of figuring out how to extract a formula data type from a lm data type

#Create varData in the global sense so any function can use it
readFile <- "Datasets/fbData.csv"
varData <- read.csv(file = readFile, header = FALSE, skip = 1)
varData$V2 = NULL
varData$V3 = NULL
varData$V4 = NULL
varData$V5 = NULL
varData$V6 = NULL
varData$V7 = NULL
names(varData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")

inform = "AIC"

getTermToDrop <- function(model){
  coefMatrix <- coef(summary(model))[2:nrow(coef(summary(model))),]
  
  drop.this <- as.integer(which(abs(coefMatrix[, 3]) == min(abs(coefMatrix[, 3]))))
}

getNameOfDropped <- function(model){
  coefMatrix <- coef(summary(model))[2:nrow(coef(summary(model))),]
  
  dropped <- rownames(coefMatrix)[which(abs(coefMatrix[, 3]) == min(abs(coefMatrix[, 3])))]
}

dropTermFromLM <- function(model){
  coefMatrix <- coef(summary(model))[2:nrow(coef(summary(model))),]
  
  drop.this <- as.integer(which(abs(coefMatrix[, 3]) == min(abs(coefMatrix[, 3]))))
  
  termsData <- terms(as.formula(model$terms))
  
  temp <- drop.terms(termsData, drop.this, keep.response = TRUE)
  
  mod.terms <- attr(temp, "term.labels")
  
  extractedFormula <- as.formula(paste(c(paste((attr(model$terms, "variables")[2]), "~ 1 "), mod.terms), collapse = " + "))
  
  newLM <- lm(extractedFormula, data = varData)
}

getAICValue <- function(model, inform) {
  switch(inform,
         AIC = AIC(model),
         BIC = BIC(model),
         PRESS = PRESS(model, as.R2 = TRUE))
}

getFullFormula <- function(model){
  terms <- terms(as.formula(model$terms))
  
  X <- attr(terms, "term.labels")
  
  d <- attr(terms, "dataClasses")[2:(length(X) + 1)]
  
  Xn <- X[d == "numeric"]
  
  Xsq <- paste("I(", Xn, "^2)", sep = "")
  
  Xint <- outer(X, X, function(x, y) paste(x, ":", y, sep = ""))
  
  Xint <- Xint[upper.tri(Xint)]
  
  Xall <- c(X, Xint, Xsq)
  
  f <- as.formula(paste(c(paste((attr(model$terms, "variables")[2]), "~ 1 "), Xall), collapse = " + "))
}

createTable <- function(model, steps){
  tab <- matrix(nrow = steps, ncol = 3)
  colnames(tab) <- c("Model", "Drop", inform)
  
  add1 <- formula(model$terms)
  
  tab[1,] <- c(paste(add1[2], add1[1], add1[3], sep = ""), " ", round(getAICValue(model, inform), digits = 5))
  
  tab
}

getGreedy <- function(startFormula){
  startFormula <- Y ~ 1 + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11
 
  startLM <- lm(startFormula, data = varData)
  
  k <- length(attr(terms(startLM), "term.labels"))
  
  steps <- 1 + (k ^ 2 + 3 * k) / 2
  
  sampleFormula <- getFullFormula(startLM)
  
  newLM <- lm(sampleFormula, data = varData)
  
  gTab <- createTable(newLM, steps)
  
  onRow <- 2
  
  maxSize <- length(attr(terms(newLM), "term.labels")) + 1
  
  for(i in 1:(length(attr(terms(newLM), "term.labels")) - 1)){
    dropped <- getNameOfDropped(newLM)
    
    newLM <- dropTermFromLM(newLM)
  
    add1 <- formula(newLM$terms)
    
    gTab[onRow,] <- c(paste(add1[2], add1[1], add1[3], sep = ""), dropped, round(getAICValue(newLM, inform), digits = 5))
    
    onRow <- onRow + 1
  }
  new.mod <- update(newLM, ~ 1)
  
  gTab[maxSize,] <- c(paste(attr(new.mod$terms, "variables")[2], "~ 1", sep = ""), rownames(coef(summary(newLM)))[2], round(getAICValue(new.mod, inform), digits = 5))
  
  if (inform == "AIC" | inform == "BIC"){ 
    opt <- which(as.numeric(gTab[, 3]) == min(as.numeric(gTab[, 3])))
  }
  
  best <- gTab[opt,][1]
  best <- lm(noquote(best), data = varData)
  res <- list()
  res$out <- data.frame(gTab)
  res$method <- inform
  res$best <- best
  res$data <- data
  class(res) <- "greedy"
  res
}

formulaObject <- Y ~ 1 + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11

greedyObject <- getGreedy(formulaObject)

print(greedyObject$out)