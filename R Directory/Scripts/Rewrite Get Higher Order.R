get.higher.order <- function(model, data) {
  res = list()
  
  terms <- terms(model.frame(model, data = data))
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
  
  for(i in 1:length(higher)){
    Xs <- unlist(strsplit(higher[i], ":"))
    loc = match(Xs, names(this.data))
    res[[i]] = levels(interaction(this.data[loc], sep = ":"))
  }
  
  names(res) = higher
  res
}

Y <- c(rnorm(36),rnorm(36,2))+rep(c(1,1.5),length=72)
XC1 <- factor(rep(c("a","b","c"), each = 24))
XC2 <- factor(rep(c("m","n"), length = 72))
XC3 <- factor(rep(rep(c("x","y"), each = 18),2))
XC4 <- factor(rep(c("k","l","p","r"), length = 72))
XQ1 <- round(runif(72, 0, 10),0)
XQ2 <- round(runif(72, 14, 999),0)
XQ3 <- round(runif(72, 600, 1200),0)

#3 quan 0 cat
# model <- lm(Y ~ XQ1 * XQ2 * XQ3)
#0 quan 3 cat
model <- lm(Y ~ XC1 * XC2 * XC3)
#3 quan 3 cat
# model <- lm(Y ~ XC1 * XC2 * XC3 * XQ1 * XQ2 * XQ3)
#2 quan 2 cat
#1 quan 1 cat
#1 quan 0 cat
# model = lm(Y ~ XQ1)
#0 quan 1 cat
# model <- lm(Y ~ XC1)

data <- get_all_vars(model)

ints = get.higher.order(model, data)
ints