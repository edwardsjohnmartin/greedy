Y <- c(rnorm(36),rnorm(36,2))+rep(c(1,1.5),length=72)
XC1 <- rep(c("a","b","c"), each = 24)
XC2 <- factor(rep(c("m","n"), length = 72))
XC3 <- rep(rep(c("x","y"), each = 18),2)
XQ1 <- round(runif(72, 0, 10),0)

model <- lm(Y ~ XC1 * XC2)


get.higher.order <- function(x) {
  data1 <- data
  data1[names(data1) == Xn] <- rep(Xn, nrow(data1))
  res <- list()
  nhigher <- length(x)
  for (i in 1:nhigher) {
    temp <- list()
    higher.comp <- strsplit(x[i], ":")
    Xs <- unlist(higher.comp)
    mat <- match(unlist(higher.comp), names(data1))
    nlvl <- apply(data1[, mat], 2, function(x) nlevels(factor(x)))
    for (j in 1:length(mat)) {
      if (nlvl[j] == 1)
        temp[[j]] <- paste(Xs[j], levels(data1[, names(data1) == Xs[j]]),
                           sep = "")
      if (nlvl[j] > 1)
        temp[[j]] <- paste(Xs[j], levels(data1[, names(data1) == Xs[j]]),
                           sep = "")[2:nlvl[j]]
    }
    temp1 <- as.matrix(interaction(temp, sep = ":"))
    res[[i]] <- as.vector(temp1)
  }
  names(res) <- x
  res
}

data <- get_all_vars(model)
m <- model.frame(model, data = data)
Y <- model.extract(m, "response")
terms <- terms(m)
X <- attr(terms,"term.labels")
d <- attr(terms,"dataClasses")[2:(length(X)+1)]
d <- d[complete.cases(d)]
Xn <- X[match(names(d)[d == "numeric"],X)]
Xc <- X[match(names(d)[d == "factor"],X)]
order.lvl <- attr(terms,"order")
higher <- X[order.lvl > 1]
no.lvls <- apply(data[,match(Xc, names(data))],2,function(x)nlevels(as.factor(x)))

if(length(higher) > 0) ints <- get.higher.order(higher)
ints

sumsat <- summary(model)$coefficients


library(rlist)
n.main <- sum(no.lvls) + length(Xn)
n.coef <- nrow(sumsat)


drop1 <- rownames(sumsat)[n.main:n.coef][which(abs(sumsat[, 3][n.main:n.coef]) ==
                                                 min(abs(sumsat[, 3][n.main:n.coef])))]
drop1


drop.lm <- names(list.search(ints, any(. == drop1)))
drop.lm

terms.new <- drop.terms(terms(model), which(attr(terms, "term.labels") == drop.lm),
                        keep.response = TRUE)
lm(terms.new)
