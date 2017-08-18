greedy <- function(formula, data = NULL, center = NULL, digits = 5, inform = "AIC"){
  require(asbio)
  data <- get_all_vars(formula, data = data)
  if(!is.null(center)){
    w <- which(names(data) == center)
    temp <- apply(data[,w], 2, function(x) x - mean(x))
    data[,w] <- temp
  } 
  m <- model.frame(formula, data = data)
  Y <- model.extract(m, "response")
  terms <- terms(m)
  X <- attr(terms,"term.labels")
  k <- length(X)
  d <- attr(terms,"dataClasses")[2:(length(X)+1)]
  Xn <- X[d == "numeric"]; Xnn <- length(Xn[!is.na(Xn)])
  Xc <- X[d == "factor"]; Xcn <- length(Xc[!is.na(Xc)])
  
  steps <- (1 + (k) + (k^2 - k)/2 + Xnn) 
  
  inf <- function(model, inform){ # choose information theoretic criterion
    switch(inform, 
           AIC = AIC(model),
           BIC = BIC(model),
           PRESS = press(model, as.R2 = TRUE))
  }
  
  tab <- matrix(nrow = steps, ncol = 3)
  colnames(tab) <- c("Model", "Drop", inform)
  test <- lm(formula, data = data)
  
  if(steps==1) tab[1,] = c(deparse(formula(test$terms)), " ", round(inf(test, inform), digits = digits)) # intercept only model
  
  else if(steps>=2){
    Xsq <- paste("I(", Xn, "^2)", sep="")
    Xint <- outer(X, X, function(x,y) paste(x,":",y,sep=""))
    Xint <- Xint[upper.tri(Xint)]
    Xall <- c(X, Xint, Xsq)
    
    # -- Categorical variables in lm --#
    if(Xcn >= 1){
      Xcm <- data.frame(data[,names(data)==Xc])
      names(Xcm) <- Xc
      nlXcm <- matrix(ncol = 1, nrow = Xcn)
      for(i in ncol(Xcm)) 
        nlXcm[i] <- nlevels(Xcm[,i])
      
      # This was already commented out. Haven't done anything with it yet
      # Xcmu <- paste(Xc, lXcm[1:nlXcm], sep="")
      
      Xcm1 <- list()
      for(i in 1:Xcn){
        Xcm1[[i]] <- paste(Xc[i],levels(data[,names(data)==Xc[i]]), sep="")[2:nlXcm[i]]
      }
      names(Xcm1) <- Xc
      
      Xclm <- 1:length(Xcm1) + 1
      
      for(i in 1:Xcn) {
        temp <- cbind(rep(names(Xcm1)[i], length(Xcm1[[i]])), Xcm1[[i]])
        if(Xcn > 1)
          Xclm <- rbind(Xclm, temp)
        if(Xcn == 1) 
          Xclm <- temp
        Xclm
      }
      
      if(nrow(Xclm)>1) 
        Xclm <- Xclm[-1,]
      Xsqlm <- cbind(as.matrix(Xsq), as.matrix(Xsq))
      Xdrop <- X[-(Xclm[,1] == X)]
      XnewL <- c(Xclm[,2], Xdrop)
      XnewF <- c(Xclm[,1], Xdrop)
      
      XintL <- outer(XnewL, XnewL, function(x,y) paste(x,":",y,sep=""))
      XintL <- XintL[upper.tri(XintL)]
      XintF <- outer(XnewF, XnewF, function(x,y) paste(x,":",y,sep=""))
      XintF <- XintF[upper.tri(XintF)]
      
      Xlm <- rbind(cbind(as.matrix(Xn), as.matrix(Xn)), Xclm, Xsqlm, cbind(as.matrix(XintF), as.matrix(XintL)))   
      colnames(Xlm) <- c("Predictor", "in_lm")
    }
    
    #-------------------------------#
    
    Yname <- names(m)[1]
    
    f <- as.formula(paste(c(paste(Yname,"~ 1 "), Xall), collapse = " + "))
    
    sat <- lm(f, data = data)
    
    # Finds the term with the lowest magnitude t-value that isn't a main effect
    # that still has instances of itself in the model
    getTermToDrop = function(sumsat, mod.terms){
      
      if(class(sumsat)=="numeric")
        sumsat = t(as.matrix(sumsat))
      
      if(nrow(sumsat)==1)
        rownames(sumsat) = mod.terms
      
      # Finds the term with the lowest magnitude t-value
      drop1 <- rownames(sumsat)[which(abs(sumsat[,3]) == min(abs(sumsat[,3])))]
      
      # If it is a main effect and has at least one other instance of itself in the model
      if(any(X==drop1) && length(grep(drop1, mod.terms[mod.terms!=drop1], fixed = TRUE))>0){
        drop2 <- getTermToDrop(sumsat[rownames(sumsat)!=drop1,], mod.terms[mod.terms!=drop1])
        drop1 = drop2
      }
      drop1
    }
    
    drops <-function(mod1, X, data){
      np <- nrow(coef(summary(mod1)))
      if(np == 2){
        new.mod <- update(mod1, ~ 1)
        res <- list(formula = paste(Yname,"~ 1"), model = new.mod, drop = attr(terms(mod1), "term.labels"), 
                    inf.crit = round(inf(new.mod, inform), digits = digits))
      }
      if(np > 2){
        # Create the coefficient matrix of model to observe t-values
        sumsat <- coef(summary(mod1))[2:np,]
        
        # Create a list of predictor variables in model
        mod.terms <- attr(terms(mod1), "term.labels")
        
        drop1 = getTermToDrop(sumsat, mod.terms)
        
        if(Xcn >= 1){
          drop1 <- Xlm[,1][Xlm[,2]==drop1]
          m <- match(rownames(sumsat),Xlm[,2]) 
          rownames(sumsat) <- Xlm[,1][m]
        }
        
        #--- find a legitimate model for dropping (there is a better way of doing this)" ----#
        
        f1 <- as.formula(paste(c(paste(Yname,"~ 1 "), mod.terms[mod.terms!=drop1]), collapse=" + "))
        
        new.mod <- lm(f1, data= data)
        res <- list(formula = paste(c(paste(Yname,"~ 1 "), mod.terms[mod.terms!=drop1]), collapse=" + "), 
                    model = new.mod, drop = drop1, inf.crit = round(inf(new.mod, inform), digits = digits))
      }
      res
    }
    
    j = 2; temp <- sat
    tab[1,] <-c(paste(c(paste(Yname,"~ 1"), Xall), collapse=" + "), " ", round(inf(temp, inform), digits = digits))
    
    while(j <= steps){
      temp <- drops(temp, X, data = data)
      tab[j,] <- c(temp$formula, temp$drop, temp$inf.crit) 
      print(tab)
      temp <- temp$model
      j = j + 1
    }
  }
  if(inform == "AIC" | inform == "BIC") opt <- which(as.numeric(tab[,3])== min(as.numeric(tab[,3])))
  if(inform == "PRESS") opt <- which(as.numeric(tab[,3]) == max(as.numeric(tab[,3])))
  best <- tab[opt,][1]
  best <- lm(noquote(best), data = data)
  res <- list()
  res$out <- data.frame(tab)
  res$method <- inform
  res$best <- best
  res$data <- data
  class(res) <- "greedy"
  res
}

readFile = "Datasets/Case0902.csv"
case0902 <- read.csv(file = readFile, header = FALSE, skip = 1)
names(case0902) <- c("Xspec", "Y", "Xb", "Xg", "Xl")
#grdy <- greedy(log(Y) ~ log(Xb) + Xg + Xl, data = case0902)
grdy <- greedy(log(Y) ~ log(Xb) + Xg + Xl + Xspec, data = case0902)

# readFile = "Datasets/FacebookData.csv"
# fbData <- read.csv(file = readFile, header = FALSE, skip = 1)
# # This is the column with categorical variable "Type" with levels - Photo, Status, Link, Video
# #fbData$V2 = NULL
# fbData$V3 = NULL
# fbData$V4 = NULL
# fbData$V5 = NULL
# fbData$V6 = NULL
# fbData$V7 = NULL
# names(fbData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13")
# grdy <- greedy(log(Y) ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13, data = fbData)

print(grdy$out)
