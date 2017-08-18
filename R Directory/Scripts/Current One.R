greedy = function(model, method = "AIC", digits = 5){
  
  inf <- function(model, inform){ # choose information theoretic criterion
    switch(inform, 
           AIC = AIC(model),
           BIC = BIC(model),
           PRESS = press(model, as.R2 = TRUE))
  }
  
  get.higher.order <- function(x, Xn, data) {
    res = list()
    data = data[-1]
    # print(data)
    data[names(data) == Xn] <- rep(Xn, nrow(data))
    # print(data)
    # print(Xn)
    if(length(Xn) == 0){
      cols = which(names(data) != "")
    } else {
      cols = which(names(data) != Xn)  
    }
    # print(cols)
    # print(paste("cols is", cols))
    for(i in 1:length(cols)){
      data[,cols[i]] = paste(names(data[cols[i]]), data[,cols[i]], sep = "")
    }
    for(i in 1:length(x)){
      Xs <- unlist(strsplit(x[i], ":"))
      loc = match(Xs, names(data))
      res[[i]] = levels(interaction(data[loc], sep = ":"))
    }
    names(res) = x
    res
  }
  
  get.names.matrix = function(){
    
    data <- get_all_vars(model)
    m <- model.frame(model, data = data)
    terms <- terms(m)
    X <- attr(terms,"term.labels")
    d <- attr(terms,"dataClasses")[2:(length(X)+1)]
    d <- d[complete.cases(d)]
    Xc <- X[match(names(d)[d == "factor"],X)]
    Xclen = length(Xc)
    if(Xclen > 1)
      no.lvls = apply(data[,match(Xc, names(data))],2,function(x) nlevels(as.factor(x)))
    if(Xclen == 1)
      no.lvls <- nlevels(as.factor(data[,match(Xc, names(data))]))
    
    maxNumber = as.numeric(no.lvls[which(no.lvls == max(no.lvls))])
    
    tempMat = matrix(nrow = Xclen, ncol = maxNumber - 1)
    rownames(tempMat) = Xc
    
    for(j in 1:Xclen){
      lv = levels(data[,colnames(data) == Xc[j]])
      lv = lv[2:length(lv)]
      for(k in 1:length(lv)){
        tempMat[j,k] = paste(Xc[j], lv[k], sep = "")
      }
    }
    tempMat
  }
  
  do.this.stuff = function(){
    print("start of do stuff")
    data <- get_all_vars(model)
    m <- model.frame(model, data = data)
    Y <- model.extract(m, "response")
    terms <- terms(m)
    X <- attr(terms,"term.labels")
    # print(X)
    d <- attr(terms,"dataClasses")[2:(length(X)+1)]
    # d <- attr(terms,"dataClasses")
    # print("d before is")
    # print(d)
    d <- d[complete.cases(d)]
    # print("d after is")
    # print(d)
    Xn <- X[match(names(d)[d == "numeric"],X)]
    Xc <- X[match(names(d)[d == "factor"],X)]
    order.lvl <- attr(terms,"order")
    higher <- X[order.lvl > 1]
    
    if(length(higher) > 0) 
      ints <- get.higher.order(higher, Xn, data)
    
    sumsat <- summary(model)$coefficients
    # print(paste("nrow in sumsat is", nrow(sumsat)))
    sumsat = sumsat[2:nrow(sumsat),]
    
    if(class(sumsat)=="numeric"){
      sumsat = t(as.matrix(sumsat))
      rownames(sumsat) = X
    }
    
    checkThis = grep(":", rownames(sumsat))
    
    tempsat = sumsat
    
    # print("------")
    # print("tempsat before is")
    # print(paste("the number of rows in tempsat is ", nrow(tempsat)))
    # print(tempsat)
    # 
    if(length(checkThis) != 0){
      deleteList = rownames(sumsat)[checkThis]
      nextList = strsplit(deleteList, ":")
      finalList = unique(unlist(nextList))
      ugg = match(finalList, rownames(sumsat))
      ugg = ugg[!is.na(ugg)]
      tempsat = sumsat[-ugg,]
    }
    
    library(rlist)
    # 
    # print("------")
    # print("tempsat after is")
    # print(paste("the number of rows in tempsat is ", nrow(tempsat)))
    # print(tempsat)
    
    drop1 <- rownames(tempsat)[which(abs(tempsat[, 3]) == min(abs(tempsat[, 3])))]
    
    if(grepl(":", drop1)){
      drop.lm <- names(list.search(ints, any(. == drop1)))
    } else if(nrow(sumsat) > 1){
      nh = match(drop1, namesMat)
      
      if(!is.na(nh)){
        jgjg = row(namesMat)
        
        jgjg = jgjg[nh]
        
        jgjg = rownames(namesMat)[jgjg]
        
        drop.lm = jgjg
      } else {
        drop.lm = drop1
      }
      
      
      if(length(drop.lm) == 0){
        drop.lm = drop1
      }
    } else{
      drop.lm = drop1
    }
    
    # if(is.na(drop.lm)){
    #   print("--------drop.lm is NA--------")
    # }
    # print("------")
    # print(paste("--------", "drop.lm is", drop.lm, "----------"))
    # print(drop.lm)
    
    print(paste("--------", "drop.lm is", drop.lm, "----------"))
    
    if(length(attr(terms, "term.labels")) > 1){
      
      lf = which(attr(terms, "term.labels") == drop.lm)
      # print("lf is")
      # print(lf)
      
      # print("---------terms(model) is---------")
      # print(terms(model))
      
      # print(summary(model))
      # print(summary(model)$coefficients)
      
      fg = attr(terms, "term.labels")
      print("fg is")
      print(fg)
      
      fr = rownames(sumsat)
      # print("fr before is")
      # print(fr)
      
      for(i in 1:length(fr)){
        uy = fr[i]
        
        if(grepl(":", uy)){
          new.uy <- names(list.search(ints, any(. == uy)))
        } else {
          nh = match(uy, namesMat)
          
          if(!is.na(nh)){
            jgjg = row(namesMat)
            
            jgjg = jgjg[nh]
            
            jgjg = rownames(namesMat)[jgjg]
            
            new.uy = jgjg
          } else {
            new.uy = uy
          }
        }
        
        fr[i] = new.uy
      }  
    
      fr = unique(fr)
      
      print("fr after is")
      print(fr)
      
      pie = which(is.na(match(fg, fr)))
      # print("pie is ")
      # print(pie)

      if(length(pie) != 0){
        term.pie = drop.terms(terms(model), pie, keep.response = TRUE)
        print("term.pie is")
        print(term.pie)
        
        term.choco = drop.terms(term.pie, which(attr(term.pie, "term.labels") == drop.lm), keep.response = TRUE)
        print("term.choco is")
        print(term.choco)
        terms.new = term.choco
      } else {
        terms.new <- drop.terms(terms(model), which(attr(terms, "term.labels") == drop.lm),
                                keep.response = TRUE)
      }
      
      
      
      # py = terms(model)
      # print("py is")
      # print(py)
      
      # py = which(grep(fg, fr) == TRUE)
      
      
      # py = grepl(fg, fr)
      # py = fg %in% fr
      # py = which(attr(terms, "term.labels") != rownames(sumsat))
      # print("py is")
      # print(py)
      
      
      
      # gt = attr(terms.new, "term.labels")
      # print("gt is")
      # print(gt)
      
      # print("-----------terms.new is----------")
      # print(terms.new)
      
      # print(data)
      
      # print(model)
      
      model = lm(terms.new)
    } else {
      model = lm(Y~1)
    }
    
    # print("model is")
    # print(model)
    
    # print("working up to here")
    
    # print("------")
    # print("model is")
    # print(model)
    
    f = formula(model)
    
    # print("------")
    # print("f is")
    # print(f)
    
    
    res = list(formula = paste(f[2], f[1], f[3]), model = model, drop = drop.lm)
    
    # print("-------res is --------")
    # print(res)
    res
  }
  
  namesMat = get.names.matrix()
  
  print("-----")
  print("namesMat is")
  print(namesMat)
  
  steps = length(attr(model$terms, "term.labels")) + 1
  
  tab <- matrix(nrow = steps, ncol = 3)
  colnames(tab) <- c("Model", "Drop", method)
  
  tab[1,] = c(as.character(model$call[2]), " ", round(inf(model, method), digits = digits))
  
  
  
  for(j in 2:steps){
    temp = do.this.stuff()
    model = temp$model
    tab[j,] = c(temp$formula, temp$drop, round(inf(temp$model, method), digits = digits))
  }
  
  tab
}

Y <- c(rnorm(36),rnorm(36,2))+rep(c(1,1.5),length=72)
XC1 <- factor(rep(c("a","b","c"), each = 24))
XC2 <- factor(rep(c("m","n"), length = 72))
XC3 <- factor(rep(rep(c("x","y"), each = 18),2))
XC4 <- factor(rep(c("k","l","p","r"), length = 72))
XQ1 <- round(runif(72, 0, 10),0)

model <- lm(Y ~ XC4 * XC1 * XC2 * XC3 * XQ1)
terms = attr(terms(model),"term.labels")

g = greedy(model)
g