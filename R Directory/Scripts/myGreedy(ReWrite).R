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