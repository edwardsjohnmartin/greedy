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
    
    dep.name = as.character(formula(model))[2]
    dep.col = which(colnames(data)==dep.name)
    this.data = data[-dep.col]
    
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
      
      test = grep(":", base.terms, invert = T)
      test1 = grep(":", corrected.terms, invert = T)
      test2 = test[-test1]
      
      if(length(test2) > 0){
        for(i in 1:length(test2)){
          test3 = grep(test2[i], corrected.terms)
          if(!is.null(test3)){
            test4 = which(na.terms == test2[i])
            na.terms = na.terms[-test4]
          }
        }
      }
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
    
    all.test.pos = unique(c(grep(":", rownames(sumsat)), grep("I", rownames(sumsat))))
    testsat = rownames(sumsat[-all.test.pos,])
    testsat1 = rownames(sumsat)[all.test.pos]
    
    test.list = list()
    
    if(!is.null(testsat)){
      for(x in 1:length(testsat)){
        name.row = match(testsat[x], namesMat)
        new.term = rownames(namesMat)[row(namesMat)[name.row]]
        testsat[x] = new.term
        
        test = grep(testsat[x], testsat1)
        test = test + length(testsat)
        
        if(length(test) == 0){
          test4 = grep(testsat[x], testsat)
          test.list[[x]] = test4
        } else {
          test.list[[x]] = test
        }
      }
    } else {
      test.list = unique(c(grep(":", rownames(sumsat), invert = T), grep("I", rownames(sumsat), invert = T)))
    } 
    
    test.list = sort(unique(unlist(test.list)))
    tempsat = sumsat[test.list,]
    if(class(tempsat)=="numeric"){
      tempsat = t(as.matrix(tempsat))
      # rownames(tempsat) = attr(terms(model), "term.labels")
      rownames(tempsat) = rownames(sumsat)[test.list]
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

      if(class(model)[1] == "lm"){
        model = lm(terms.new, data = data)
      } 
      else if(class(model)[1] == "glm"){
        model = glm(terms.new, family = "binomial", data = data)
        f = formula(model)
      }
    } else {
      if(class(model)[1] == "lm"){
        dep.name = as.character(formula(model))[2]
        new.form = as.formula(paste(dep.name, "~", 1, sep = ""))
        model = lm(new.form, data = data)
      } 
      else if(class(model)[1] == "glm"){
        dep.name = as.character(formula(model))[2]
        dep.col = which(colnames(data)==dep.name)
        
        new.form = as.formula(paste(dep.name, "~", 1, sep = ""))
        model = glm(new.form, family = "binomial", data = data)
        
        f = formula(model)
      }
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
  terms = terms(model)
  X <- attr(terms,"term.labels")
  d <- attr(terms,"dataClasses")[2:(length(X)+1)]
  d <- d[complete.cases(d)]
  Xn <- X[match(names(d)[d == "numeric"],X)]
  Xc <- X[match(names(d)[d == "factor"],X)]
  order.lvl <- attr(terms,"order")
  higher <- X[order.lvl > 1]
  if(length(X) > 0){
    namesMat = get.names.matrix(model, data)
  }
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
  
  if(class(model)[1] == "lm"){
    best <- lm(noquote(best), data = data)
  } 
  else if(class(model)[1] == "glm"){
    best <- glm(noquote(best), family = "binomial", data = data)
  }
  
  res <- list()
  res$out <- data.frame(tab)
  res$method <- inform
  res$best <- best
  res$data <- data
  class(res) <- "greedy"
  res
}

test.models = function(){
  library("MASS")
  
  # fb.names = c("Y","1C","2C","3C","4C","5C","6C","V01","V02","V03","V04","V05","V06","V07","V08","V09","V10","V11","V12")
  # fb.colclasses = c("X1C"="factor","X2C"="factor","X3C"="factor","X4C"="factor","X5C"="factor","X6C"="factor")
  # data = read.csv(file = "Datasets/FacebookData.csv", colClasses = fb.colClasses, col.names = fb.names)

  # titan.names = c("pclass", "survived", "name", "sex", "age", "sibsp", "parch", "ticket", "fare", "cabin", "embarked", "boat", "body", "homedest")
  # titan.colclasses = c("survived"="factor","sex"="factor","pclass"="factor","sibsp"="factor","parch"="factor", "embarked" = "factor", "boat" = "factor")
  # data = read.csv(file = "Datasets/titanicData-Binom.csv", colClasses = titan.colclasses, col.names = titan.names)

  # wine.names = c("qual", "vol.acid", "cit.acid", "res.sugar", "chlor", "free.sulfur", "total.sulfur", "dens", "PH", "sulph", "alcohol", "fixed.acid")
  # data = read.csv(file = "Datasets/whiteWineData.csv")
  
  # parkin.names = c("subject", "age", "sex", "test.time", "motor.updrs", "total.updrs", "jitter.percent", "jitter.abs", "jitter.rap", "jitter.ppq5", "jitter.ddp", "shim", "shim.db", "shim.apq3", "shim.apq5", "shim.apq11", "shim.dda", "nhr", "hnr", "rdpe", "dfa", "ppe")
  # parkin.colclasses = c("sex"="factor")
  # data = read.csv(file = "Datasets/parkinsonsData.csv", col.names = parkin.names, colClasses = parkin.colclasses)
  # model = lm(total.updrs ~ (age + sex + test.time + jitter.percent + jitter.abs + jitter.rap + jitter.ppq5 + jitter.ddp + shim + shim.db + shim.apq3 + shim.dda + nhr + hnr + rdpe + dfa)^3, data = data)
  
  # concrete.names = c("cement", "blast", "fly.ash", "water", "superplas", "coarse.agg", "fine.agg", "age", "compress.str")
  # data = read.csv(file = "Datasets/concreteData.csv", col.names = concrete.names)
  # model = lm(cement ~ (blast + fly.ash + water + superplas + coarse.agg + fine.agg + age + compress.str)^4, data = data)
  
  # power.names = c("temp", "vacuum", "ambient.pres", "relative.humid", "energy.out")
  # data = read.csv(file = "Datasets/powerPlantData.csv", col.names = power.names)
  # model = lm(energy.out ~ (temp + vacuum + ambient.pres + relative.humid), data = data)
  
  data = read.csv(file = "Datasets/OnlineNewsPopularityData.csv")
  model = lm(shares ~ (n_tokens_title + n_tokens_content + n_unique_tokens + n_non_stop_words + n_non_stop_unique_tokens + num_hrefs +  num_self_hrefs +  num_imgs +  num_videos +  average_token_length +  num_keywords)^3, data = data)
  
  # grades.colClasses = c("sex"="factor","address"="factor", "Pstatus"="factor", "Medu"="factor", "Fedu"="factor", "schoolsup" = "factor")
  # data = read.csv(file = "Datasets/StudentGrades.csv", colClasses = grades.colClasses)
  # model = lm(G1 ~ (sex + age + address + Pstatus + Medu + Fedu + traveltime + studytime + failures + schoolsup + absences)^3,data = data)
  
  # fires.colClasses = c("X"="factor", "Y"="factor", "month"="factor", "day"="factor")
  # data = read.csv(file = "Datasets/ForestFires.csv", colClasses = fires.colClasses)
  # model = lm(area ~ (X + Y + month + day + FFMC + DMC + DC + ISI + temp + RH + wind + rain)^3, data = data)
  
  # greedy.test = greedy(model, data = data)
  # # print("-----------greedy output--------------")
  # # print(greedy.test$out)
  # print("------------greedy AIC----------------")
  # print(AIC(greedy.test$best))
  # # print("-----------test for class-------------")
  # # print(class(greedy.test$best))
  # print("------------greedy Total Time---------")
  # print(system.time(greedy(model, data = data)))
  
  stepaic.test = stepAIC(model, data = data, trace = F)
  # print("-----------stepAIC output-------------")
  # print(stepaic.test)
  print("------------stepAIC AIC---------------")
  print(AIC(stepaic.test))
  print("------------stepAIC Total Time--------")
  print(system.time(stepAIC(model, data = data, trace = F)))
}
test.models()