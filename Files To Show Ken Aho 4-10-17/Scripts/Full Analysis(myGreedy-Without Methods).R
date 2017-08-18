#Will produce timing results of the greedy(Without Methods) and stepAIC algorithms
#Currently only fbData and parkinsonsData will work 

#----------Change this line only--------------
readFile <- "Datasets/parkinsonsData.csv"
#---------------------------------------------

varData <- read.csv(file = readFile, header = FALSE, skip = 1) 

if (readFile == "Datasets/case1202.csv") {
  varData$V2 = NULL
  varData$V3 = NULL 
  names(varData) <- c("Y", "V1", "V2", "V3", "V4")
} else if (readFile == "Datasets/case0902.csv") {
  varData$V1 = NULL
  names(varData) <- c("Y", "V1", "V2", "V3")
} else if (readFile == "Datasets/fbData.csv") {
  varData$V2 = NULL
  varData$V3 = NULL
  varData$V4 = NULL
  varData$V5 = NULL
  varData$V6 = NULL
  varData$V7 = NULL
  names(varData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
} else if (readFile == "Datasets/concreteData.csv") {
  names(varData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8")
} else if (readFile == "Datasets/whiteWineData.csv") {
  names(varData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11")
} else if (readFile == "Datasets/powerPlantData.csv") {
  names(varData) <- c("V1", "V2", "V3", "V4", "Y")
} else if (readFile == "Datasets/parkinsonsData.csv") {
  varData$V1 = NULL
  varData$V2 = NULL
  varData$V3 = NULL
  varData$V4 = NULL
  # This dataset has 2 predictor variables, uncomment the next line and comment the line after to change the predictor variable
  #varData$V5 = NULL
  varData$V6 = NULL
  names(varData) <- c( "Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16")
} else if (readFile == "Datasets/parkinsonsData(100).csv") {
  varData$V1 = NULL
  varData$V2 = NULL
  varData$V3 = NULL
  varData$V4 = NULL
  # This dataset has 2 predictor variables, uncomment the next line and comment the line after to change the predictor variable
  #varData$V5 = NULL
  varData$V6 = NULL
  names(varData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16")
}

greedy <- function(formula, data = NULL, center = NULL, inform = "AIC"){
  getAICValue <- function(model, inform) {
    switch(inform,
           AIC = AIC(model),
           BIC = BIC(model),
           PRESS = PRESS(model, as.R2 = TRUE))
  }
  
  startLM <- lm(formula, data = varData)
  
  k <- length(attr(terms(startLM), "term.labels"))
  
  steps <- 1 + (k ^ 2 + 3 * k) / 2
  
  terms <- terms(as.formula(startLM$terms))
  
  X <- attr(terms, "term.labels")
  
  d <- attr(terms, "dataClasses")[2:(length(X) + 1)]
  
  Xn <- X[d == "numeric"]
  
  Xsq <- paste("I(", Xn, "^2)", sep = "")
  
  Xint <- outer(X, X, function(x, y) paste(x, ":", y, sep = ""))
  
  Xint <- Xint[upper.tri(Xint)]
  
  Xall <- c(X, Xint, Xsq)
  
  sampleFormula <- as.formula(paste(c(paste((attr(startLM$terms, "variables")[2]), "~ 1 "), Xall), collapse = " + "))
  
  newLM <- lm(sampleFormula, data = varData)
  
  gTab <- matrix(nrow = steps, ncol = 3)
  colnames(gTab) <- c("Model", "Drop", inform)
  
  add1 <- formula(newLM$terms)
  
  gTab[1,] <- c(paste(add1[2], add1[1], add1[3], sep = ""), " ", round(getAICValue(newLM, inform), digits = 5))
  
  onRow <- 2
  
  maxSize <- length(attr(terms(newLM), "term.labels")) + 1
  
  for(i in 1:(length(attr(terms(newLM), "term.labels")) - 1)){
    coefMatrix <- coef(summary(newLM))[2:nrow(coef(summary(newLM))),]
    
    dropped <- rownames(coefMatrix)[which(abs(coefMatrix[, 3]) == min(abs(coefMatrix[, 3])))]
    
    drop.this <- as.integer(which(abs(coefMatrix[, 3]) == min(abs(coefMatrix[, 3]))))
    
    termsData <- terms(as.formula(newLM$terms))
    
    temp <- drop.terms(termsData, drop.this, keep.response = TRUE)
    
    mod.terms <- attr(temp, "term.labels")
    
    extractedFormula <- as.formula(paste(c(paste((attr(newLM$terms, "variables")[2]), "~ 1 "), mod.terms), collapse = " + "))
    
    newLM <- lm(extractedFormula, data = varData)
    
    add1 <- formula(newLM$terms)
    
    gTab[onRow,] <- c(paste(add1[2], add1[1], add1[3], sep = ""), dropped, round(getAICValue(newLM, inform), digits = 5))
    
    onRow <- onRow + 1
  }
  
  new.mod <- update(newLM, ~ 1)
  
  gTab[maxSize,] <- c(paste(attr(new.mod$terms, "variables")[2], "~ 1", sep = ""), rownames(coef(summary(newLM)))[2], round(getAICValue(new.mod, inform), digits = 5))
  
  if (inform == "AIC" | inform == "BIC"){ 
    opt <- which(as.numeric(gTab[, 3]) == min(as.numeric(gTab[, 3])))
  }
  
  print(formula)
  
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

times13 <- function() {
  require(MASS)
  #sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
  sg0 <-0.00
  ss0 <- system.time(lm(Y ~ 1, data = varData))[3]
  
  sg1 <- system.time(greedy(Y ~ V1, data = varData))[3]
  ss1 <- system.time(stepAIC(lm(Y ~ 1 + V1 + I(V1 ^ 2), data = varData), trace = F))[3]
  
  sg2 <- system.time(greedy(Y ~ V1 + V2, data = varData))[3]
  ss2 <- system.time(stepAIC(lm(Y ~ 1 + V1 + V2 + V1:V2 + I(V1 ^ 2) + I(V2 ^ 2), data = varData), trace = F))[3]
  
  sg3 <- system.time(greedy(Y ~ V1 + V2 + V3, data = varData))[3]
  ss3 <- system.time(stepAIC(lm(Y ~ 1 + I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3) ^ 2, data = varData), trace = F))[3]
  
  sg4 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4, data = varData))[3]
  ss4 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3 + V4) ^ 2, data = varData), trace = F))[3]
  
  sg5 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5, data = varData))[3]
  ss5 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + (V1 + V2 + V3 + V4 + V5) ^ 2, data = varData), trace = F))[3]
  
  sg6 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6, data = varData))[3]
  ss6 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6) ^ 2, data = varData), trace = F))[3]
  
  sg7 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7, data = varData))[3]
  ss7 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7) ^ 2, data = varData), trace = F))[3]
  
  sg8 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8, data = varData))[3]
  ss8 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8) ^ 2, data = varData), trace = F))[3]
  
  sg9 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9, data = varData))[3]
  ss9 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9) ^ 2, data = varData), trace = F))[3]
  
  sg10 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = varData))[3]
  ss10 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10) ^ 2, data = varData), trace = F))[3]
  
  sg11 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = varData))[3]
  ss11 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11) ^ 2, data = varData), trace = F))[3]
  
  sg12 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, data = varData))[3]
  ss12 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12) ^ 2, data = varData), trace = F))[3]
  
  sg <- c(sg0, sg1, sg2, sg3, sg4, sg5, sg6, sg7, sg8, sg9, sg10, sg11, sg12)
  ss <- c(ss0, ss1, sg2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11, ss12)
  list(sg = sg, ss = ss)
}

times17 <- function() {
  require(MASS)
  #sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
  sg0 <- 0.00
  ss0 <- system.time(lm(Y ~ 1, data = varData))[3]
  
  sg1 <- system.time(greedy(Y ~ V1, data = varData))[3]
  ss1 <- system.time(stepAIC(lm(Y ~ 1 + V1 + I(V1 ^ 2), data = varData), trace = F))[3]
  
  sg2 <- system.time(greedy(Y ~ V1 + V2, data = varData))[3]
  ss2 <- system.time(stepAIC(lm(Y ~ 1 + V1 + V2 + V1:V2 + I(V1 ^ 2) + I(V2 ^ 2), data = varData), trace = F))[3]
  
  sg3 <- system.time(greedy(Y ~ V1 + V2 + V3, data = varData))[3]
  ss3 <- system.time(stepAIC(lm(Y ~ 1 + I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3) ^ 2, data = varData), trace = F))[3]
  
  sg4 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4, data = varData))[3]
  ss4 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3 + V4) ^ 2, data = varData), trace = F))[3]
  
  sg5 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5, data = varData))[3]
  ss5 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + (V1 + V2 + V3 + V4 + V5) ^ 2, data = varData), trace = F))[3]
  
  sg6 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6, data = varData))[3]
  ss6 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6) ^ 2, data = varData), trace = F))[3]
  
  sg7 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7, data = varData))[3]
  ss7 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7) ^ 2, data = varData), trace = F))[3]
  
  sg8 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8, data = varData))[3]
  ss8 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8) ^ 2, data = varData), trace = F))[3]
  
  sg9 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9, data = varData))[3]
  ss9 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9) ^ 2, data = varData), trace = F))[3]
  
  sg10 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = varData))[3]
  ss10 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10) ^ 2, data = varData), trace = F))[3]
  
  sg11 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = varData))[3]
  ss11 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11) ^ 2, data = varData), trace = F))[3]
  
  sg12 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, data = varData))[3]
  ss12 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12) ^ 2, data = varData), trace = F))[3]
  
  sg13 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13, data = varData))[3]
  ss13 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + I(V13 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13) ^ 2, data = varData), trace = F))[3]
  
  sg14 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14, data = varData))[3]
  ss14 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + I(V13 ^ 2) + I(V14 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14) ^ 2, data = varData), trace = F))[3]
  
  sg15 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15, data = varData))[3]
  ss15 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + I(V13 ^ 2) + I(V14 ^ 2) + I(V15 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15) ^ 2, data = varData), trace = F))[3]
  
  sg16 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16, data = varData))[3]
  ss16 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + I(V13 ^ 2) + I(V14 ^ 2) + I(V15 ^ 2) + I(V16 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16) ^ 2, data = varData), trace = F))[3]
  
  sg <- c(sg0, sg1, sg2, sg3, sg4, sg5, sg6, sg7, sg8, sg9, sg10, sg11, sg12, sg13, sg14, sg15, sg16)
  ss <- c(ss0, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11, ss12, ss13, ss14, ss15, ss16)
  list(sg = sg, ss = ss)
}

nmodels <- function(k) {
  Mk <- matrix(ncol = 1, nrow = length(seq(0:k)))
  for (j in 0:k) {
    Mk[j + 1] <- choose(k, j) * 2 ^ ((j ^ 2 + j) / 2)
  }
  if (k == 0) Mk = 1
  sum(Mk)
}

NeedsToDoAllThis <- function(){
  len <- length(varData)
  sapply(0:len, nmodels)
  sim <- 1
  time <- array(dim = c(len, 2, sim), dimnames = list(0:(len-1), c("ss", "sg"), 1:sim))
  for (i in 1:sim) {
    if (len == 4) {
      timing <- times4()
      ss <- timing$ss
      sg <- timing$sg
    } else if (len == 5) {
      timing <- times5()
      ss <- timing$ss
      sg <- timing$sg
    } else if (len == 9) {
      timing <- times9()
      ss <- timing$ss
      sg <- timing$sg
    } else if (len == 12) {
      timing <- times12()
      ss <- timing$ss
      sg <- timing$sg
    } else if (len == 13) {
      timing <- times13()
      ss <- timing$ss
      sg <- timing$sg
    } else if (len == 17) {
      timing <- times17()
      ss <- timing$ss
      sg <- timing$sg
    }
    time[,, i][, 1] <- ss
    time[,, i][, 2] <- sg
  }
  print(paste(readFile, " - Without Methods in myGreedy", sep = ""))
  
  print('The time it takes to evaluate each regression line using the Greedy algorithm is')
  print(sg)
  print('The time it takes to evaluate each regression line using the StepAIC algorithm is')
  print(ss)
}

NeedsToDoAllThis()