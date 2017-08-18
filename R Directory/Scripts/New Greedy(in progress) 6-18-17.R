# Agenda:
# cannot drop a main effect if a squared term of it exists (complete)
# cannot drop a squared term if it exists in another interaction (complete)
# Incorporate glim models
# Squared terms with interactions(mostly complete, more testing required)
# cubed terms with interactions
# library(nlme), other model to worry about, not a big deal yet
# put centering of data back in(complete)
# check to see if greedy will work with data where the response variable isn't the first column
# find out if centering is an all or nothing, or if some columns will be centered and some won't be

#6/21/17
#worked on getting the drop order correct
#now able to drop the squared term before the main effect
#ran into issue where if a squared term is dropped but it still exists in an interaction
#then an error would occur

#6/22/17
#finished getting the drop order correct
#worked on getting squared terms to work
#got it pretty much working all the way with quantitative variables
#working on doing interactions between categorical variables and squared quantitative variables
#currently, ints is not filling in properly
#it is showing 
#$`XC1:XQ1:I(XQ1^2)`
#[1] "XC1:XQ1:I(XQ1^2)"
#instead of doing it for each level of XC1 (levels = "a", "b", "c")
#this is an issue in get.higher.order

#6/24/2017
#got ints to work properly with squared terms 
#did this by adding a column to this.data in get.higher.order for each squared term
#this allows the levels(interaction(this.data)) function to work and properly fill ints
#works with models having multiple categorical variables, quantitative variables, and squared terms

#6/25/17
#worked on getting data centered
#ran into issue where having cat variables in data would break the centering(fixed)
#data is now centered correctly
#bug fixing that occured when testing various models
#changed some variable names
#removed some extraneous variables

#6/26/17
#started doing tests to see if greedy will work with real data
#had to change some things around for it to work
#centering data works as well with real data
#every dataset tried so far has worked(case1202, facebookData, concreteData)
#most notably, concreteData worked using all 7 columns interacted with each other plus 2 columns squared (585 steps)
#it takes a little while to run but it works
#also worth noting that the output is way too long to be displayed 
#timing of concrete data with 7 predvars and 2 squared terms all interacted
#elapsed      center    gcFirst
#62.61        TRUE      FALSE
#63.00        TRUE      TRUE
#62.79        FALSE     FALSE
#62.35        FALSE     TRUE

#6/27/17
#StepAIC timing with a model using concreteData and all 7 columns interacted with each other plus 2 columns squared
#elapsed 
#22980 
greedy = function(model, data = NULL, center = FALSE, inform = "AIC", digits = 5){
  inf <- function(model, inform){ # choose information theoretic criterion
    switch(inform, 
           AIC = AIC(model),
           BIC = BIC(model),
           PRESS = press(model, as.R2 = TRUE))
  }
  
  get.higher.order <- function(model, data) {
    res = list()
    
    terms <- terms(model)
    X <- attr(terms,"term.labels")
    d <- attr(terms,"dataClasses")[2:(length(X)+1)]
    d <- d[complete.cases(d)]
    Xn <- X[match(names(d)[d == "numeric"],X)]
    Xc <- X[match(names(d)[d == "factor"],X)]
    order.lvl <- attr(terms,"order")
    higher <- X[order.lvl > 1]

    this.data = data[-1]
    
    cat.col = which(is.na(match(names(this.data), Xn)))
    quan.col = which(!is.na(match(names(this.data), Xn)))

    if(length(cat.col) > 0){
      this.data[,cat.col] = mapply(paste, names(this.data)[cat.col], this.data[,cat.col], sep = "")
    }
    if(length(quan.col) > 0){
      this.data[,quan.col] = mapply(rep, names(this.data)[quan.col], nrow(this.data))
    }

    sq.vars = Xn[which(grepl("I", Xn) == TRUE)]
    start = (ncol(this.data) + 1)
    sq.col = rep(start:(start + length(sq.vars)), length.out = length(sq.vars))
    
    if(length(sq.vars) != 0){
      this.data[,sq.col] = mapply(rep, sq.vars, nrow(this.data))
      colnames(this.data)[sq.col] = sq.vars
    }
    
    names.data = names(this.data)
    
    for(i in 1:length(higher)){
      Xs <- unlist(strsplit(higher[i], ":"))
      loc = match(Xs, names.data)
      res[[i]] = levels(interaction(this.data[loc], sep = ":"))
    }
    names(res) = higher
    res
  }
  
  get.names.matrix = function(model, data){
    X <- attr(terms(model),"term.labels")
    d <- attr(terms(model),"dataClasses")[2:(length(X)+1)]
    d <- d[complete.cases(d)]
    Xn <- X[match(names(d)[d == "numeric"],X)]
    Xc <- X[match(names(d)[d == "factor"],X)]
    XnLen = length(Xn)
    XcLen = length(Xc)
    totalLen = XnLen + XcLen
    
    if(XcLen > 1)
      no.lvls = apply(data[,match(Xc, names(data))],2,function(x) nlevels(as.factor(x)))
    if(XcLen == 1)
      no.lvls <- nlevels(as.factor(data[,match(Xc, names(data))]))
    
    if(exists("no.lvls")){
      maxNumber = as.numeric(no.lvls[which(no.lvls == max(no.lvls))])
    } else {
      maxNumber = 1
    }
    
    tempMat = matrix(nrow = totalLen, ncol = maxNumber)
    rownames(tempMat) = c(Xc, Xn)
    
    if(XcLen > 0){
      for(j in 1:XcLen){
        lv = paste(Xc[j], levels(data[,colnames(data) == Xc[j]]), sep = "")
        tempMat[j,1:length(lv)] = lv
      }
    }
    
    if(XnLen > 0){
      tempMat[(XcLen + 1):totalLen, 1] = Xn
    }
    
    tempMat
  }
  
  get.clean.model = function(model){
    if(length(attr(terms(model), "term.labels")) > 1){
      sumsat <- summary(model)$coefficients
      sumsat = sumsat[2:nrow(sumsat),]
      base.terms = attr(terms(model), "term.labels")
      levels.terms = rownames(sumsat)
      
      for(i in 1:length(levels.terms)){
        if(grepl(":", levels.terms[i])){
          new.term = names(list.search(ints, any(. == levels.terms[i])))
        } else {
          name.row = match(levels.terms[i], namesMat)
          new.term = rownames(namesMat)[row(namesMat)[name.row]]
        }
        levels.terms[i] = new.term
      }
      corrected.terms = unique(levels.terms)
      na.terms = which(is.na(match(base.terms, corrected.terms)))
      
      if(length(na.terms) != 0){
        model = lm(drop.terms(terms(model), na.terms, keep.response = TRUE), data = data)
      } 
    }
    model
  }
  
  do.this.stuff = function(){
    sumsat <- summary(model)$coefficients
    sumsat = sumsat[2:nrow(sumsat),]
    
    if(class(sumsat)=="numeric"){
      sumsat = t(as.matrix(sumsat))
      rownames(sumsat) = attr(terms(model), "term.labels")
    }
    
    if(length(which(is.na(sumsat[,3]) == FALSE)) == 0){
      stop("Flawed model, all t-values are NA")
    }
    
    tempsat = sumsat

    all.pos = unique(c(grep(":", rownames(sumsat)), grep("I", rownames(sumsat))))
    all.names = rownames(sumsat)
    
    if(length(all.pos) != 0){
      main.terms = all.names[-all.pos]
      keep.List = rownames(sumsat)[all.pos]
      int.terms = keep.List[-grep(":", keep.List)]

      remove.terms = list()
      i = 1
      
      for(x in 1:length(main.terms)){
        if(any(grepl(main.terms[x], keep.List)) == TRUE){
          remove.terms[i] = which(all.names == main.terms[x])
          i = i + 1
        } 
      }
      if(length(int.terms) != 0){
        for(x in 1:length(int.terms)){
          if(any(grepl(int.terms[x], keep.List, fixed = TRUE)) == TRUE){
            remove.terms[i] = which(all.names == int.terms[x])
            i = i + 1
          } 
        }
      }
      
      remove.terms = unlist(remove.terms)
      
      if(!is.null(remove.terms))
        tempsat = tempsat[-remove.terms,]

      if(class(tempsat) == "numeric"){
        tempsat = t(as.matrix(tempsat))
        rownames(tempsat) = rownames(sumsat)[-remove.terms]
      }
    } 
    
    drop1 <- rownames(tempsat)[which(abs(tempsat[, 3]) == min(abs(tempsat[, 3])))]
    if(grepl(":", drop1)){
      if(exists("ints")){
        drop.lm <- names(list.search(ints, any(. == drop1)))
      } else {
        drop.lm = drop1
      }
    } else {
      if(nrow(tempsat) == 1){
        drop.lm = drop1
      } else {
        name.row = match(drop1, namesMat)
        name.row = row(namesMat)[name.row]
        drop.lm = rownames(namesMat)[name.row]
      }
    }
    
    if(length(attr(terms(model), "term.labels")) > 1){
      terms.new <- drop.terms(terms(model), which(attr(terms(model), "term.labels") == drop.lm), keep.response = TRUE)
      model = lm(terms.new, data = data)
    } else {
      model = lm(Y~1)
    }
    
    f = formula(model)
    res = list(formula = paste(f[2], f[1], f[3]), model = model, drop = drop.lm)
    res
  }
  
  # --------- Start of Function ---------- #
  library(rlist)
  
  if(is.null(data)){
    data <- get_all_vars(model)
  }

  if(center == TRUE){
    num.cols = which(sapply(data, is.numeric) == TRUE)
    num.cols = unname(num.cols[which(num.cols != attr(terms(model), "response"))])
    scaled.data = scale(data[num.cols], scale = FALSE)
    data[,num.cols] = scaled.data
  }

  terms <- terms(model.frame(model, data = data))
  X <- attr(terms,"term.labels")
  Y <- model.extract(model.frame(model, data = data), "response")
  d <- attr(terms,"dataClasses")[2:(length(X)+1)]
  d <- d[complete.cases(d)]
  Xn <- X[match(names(d)[d == "numeric"],X)]
  Xc <- X[match(names(d)[d == "factor"],X)]
  order.lvl <- attr(terms,"order")
  higher <- X[order.lvl > 1]
  
  namesMat = get.names.matrix(model, data)
  if(length(higher) > 0){
    ints = get.higher.order(model, data)
  }
  model = get.clean.model(model)
  
  steps = length(attr(model$terms, "term.labels")) + 1
  tab <- matrix(nrow = steps, ncol = 3)
  colnames(tab) <- c("Model", "Drop", inform)
  f = formula(model)
  tab[1,] = c(paste(f[2], f[1], f[3]), " ", round(inf(model, inform), digits = digits))

  if(steps > 1){
    for(j in 2:steps){
      temp = do.this.stuff()
      model = temp$model
      tab[j,] = c(temp$formula, temp$drop, round(inf(temp$model, inform), digits = digits))
    }
  }
  if(inform == "AIC" | inform == "BIC") 
    opt <- which(as.numeric(tab[,3])== min(as.numeric(tab[,3])))
  if(inform == "PRESS") 
    opt <- which(as.numeric(tab[,3]) == max(as.numeric(tab[,3])))
  best <- tab[opt,][1]
  best <- lm(noquote(best), data = data)
  res <- list()
  res$out <- data.frame(tab)
  res$method <- inform
  res$best <- best
  res$data <- data
  class(res) <- "greedy"
  # print(res$out)
  res
}

library("MASS")

Y <- c(rnorm(36),rnorm(36,2))+rep(c(1,1.5),length=72)
XC1 <- factor(rep(c("a","b","c"), each = 24))
XC2 <- factor(rep(c("m","n"), length = 72))
XC3 <- factor(rep(rep(c("x","y"), each = 18),2))
XC4 <- factor(rep(c("k","l","p","r"), length = 72))
XQ1 <- round(runif(72, 0, 10),0)
XQ2 <- round(runif(72, 14, 999),0)
XQ3 <- round(runif(72, 600, 1200),0)

#1 cat 1 quan 1 sq
# model <- lm(Y ~ XC1 * XQ1 * I(XQ1^2))
#1 cat 2 quan 1 sq
# model <- lm(Y ~ XC1 * XQ1 * XQ2 * I(XQ1^2))
#1 cat 2 quan 2 sq
# model <- lm(Y ~ XC1 * XQ1 * XQ2 * I(XQ1^2) * I(XQ2^2))
#2 cat 2 quan 1 sq
# model <- lm(Y ~ XC1 * XC2 * XQ1 * XQ2 * I(XQ1^2))
#2 cat 2 quan 2 sq
# model <- lm(Y ~ XC1 * XC2 * XQ1 * XQ2 * I(XQ1^2) * I(XQ2^2))
#2 quan 1 sq term ints
# model <- lm(Y ~ XQ1 * XQ2 * I(XQ1^2))
#2 quan 2 sq term ints
# model <- lm(Y ~ XQ1 * XQ2 * I(XQ1^2) * I(XQ2^2))
#3 quan 1 sq term ints
# model <- lm(Y ~ XQ1 * XQ2 * XQ3 * I(XQ1^2))
#3 quan 2 sq term ints
# model <- lm(Y ~ XQ1 * XQ2 * XQ3 * I(XQ1^2) * I(XQ2^2))
#3 quan 3 sq term ints
# model <- lm(Y ~ XQ1 * XQ2 * XQ3 * I(XQ1^2) * I(XQ2^2) * I(XQ3^2))
#3 quan 0 cat
# model <- lm(Y ~ XQ1 * XQ2 * XQ3)
# 0 quan 3 cat
# model <- lm(Y ~ XC1 * XC2 * XC3)
#3 quan 3 cat
# model <- lm(Y ~ XC1 * XC2 * XC3 * XQ1 * XQ2 * XQ3)
#2 quan 2 cat
# model <- lm(Y ~ XC1 * XC2 * XQ1 * XQ2)
#1 quan 1 cat
# model <- lm(Y ~ XC1 * XQ1)
#1 quan 0 cat
# model = lm(Y ~ XQ1)
#0 quan 1 cat
# model <- lm(Y ~ XC1)

#case1202 dataset
# case1202 = read.csv(file = "Datasets/Case1202.csv")
# model = lm(BSAL ~ SAL77 * SEX * SENIOR * AGE * EDUC * EXPER, data = case1202)
# g = greedy(model, data = case1202)

#concrete dataset initialize
# conData = read.csv(file = "Datasets/concreteData.csv")
# names(conData) = c("Y", "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8")
# model = lm(Y ~ I(X1^2) + I(X2^2) + I(X3^2) + I(X4^2) + I(X5^2) + I(X6^2) + I(X7^2) + I(X8^2) + (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8)^3, data = conData)

#concrete dataset tests
# concrete.greedy = greedy(model, data = conData)
# concrete.step = stepAIC(model, data = conData, trace = F)
# 
# system.time(greedy(model, data = conData))
# system.time(stepAIC(model, data = conData, trace = F))
# 
# AIC(concrete.greedy$best)
# AIC(concrete.step)
# 
# concrete.greedy$best
# concrete.step

#facebook dataset initialize
fb.names = c("Y" , "C1","C2","C3","C4","C5","C6","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12")
fb.colclasses = c("C1"="factor","C2"="factor","C3"="factor","C4"="factor","C5"="factor","C6"="factor")
fbData = read.csv(file = "Datasets/FacebookData.csv", colClasses = fb.colclasses, col.names = fb.names)
# model = lm(Y ~ C1 + C2 + C3 + C4 + C5 + C6, data = fbData)
# model = lm(Y ~ (C1 + C2 + C3 + C4 + C5 + C6)^2, data = fbData)
# model = lm(Y ~ C1 + C2 + C3 + C4 + C5 + C6 + V01 + V02 + V03 + V04 + V05 + V06 + V07 + V08 + V09 + V10 + V11 + V12, data = fbData)
model = lm(Y ~ C1 + C2 + C3 + C4 + C5 + C6 + (V01 + V02 + V03 + V04 + V05 + V06 + V07 + V08 + V09 + V10 + V11 + V12)^2, data = fbData)

#facebook dataset tests
fb.greedy = greedy(model, data = fbData)
fb.step = stepAIC(model, data = fbData, trace = F)

system.time(greedy(model, data = fbData))
system.time(stepAIC(model, data = fbData, trace = F))

AIC(fb.greedy$best)
AIC(fb.step)

fb.greedy$best
fb.step