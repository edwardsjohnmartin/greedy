greedy <- function(formula, data = NULL, center = NULL, digits = 5, inform = "AIC") {
    require(asbio)
    data <- get_all_vars(formula, data = data)
    if (!is.null(center)) {
        w <- which(names(data) == center)
        temp <- apply(data[, w], 2, function(x) x - mean(x))
        data[, w] <- temp
    }
    m <- model.frame(formula, data = data)
    Y <- model.extract(m, "response")
    terms <- terms(m)
    X <- attr(terms, "term.labels")
    k <- length(X)
    steps <- 1 + (k ^ 2 + 3 * k) / 2

    inf <- function(model, inform) {
        switch(inform,
        AIC = AIC(model),
        BIC = BIC(model),
        PRESS = PRESS(model, as.R2 = TRUE))
    }

    tab <- matrix(nrow = steps, ncol = 3)
    colnames(tab) <- c("Model", "Drop", inform)
    test <- lm(formula, data = data)

    if (steps == 1) tab[1,] = c(deparse(formula(test$terms)), " ", round(inf(test, inform), digits = digits))

    else if (steps >= 2) {
        d <- attr(terms, "dataClasses")[2:(length(X) + 1)]
        Xn <- X[d == "numeric"]
        Xsq <- paste("I(", Xn, "^2)", sep = "")
        Xint <- outer(X, X, function(x, y) paste(x, ":", y, sep = ""))
        Xint <- Xint[upper.tri(Xint)]
        Xall <- c(X, Xint, Xsq)
        if (!any(match(Xn, X))) Xall <- c(X, Xint)
        Yname <- names(m)[1]

        f <- as.formula(paste(c(paste(Yname, "~ 1 "), Xall), collapse = " + "))
        sat <- lm(f, data = data)

        redo <- function(drop1, sumsat) {
            sumsat1 <- sumsat[rownames(sumsat) != drop1,]
            if (class(sumsat1) == "numeric") sumsat1 = t(as.matrix(sumsat1))
            rn <- rownames(sumsat)[rownames(sumsat) != drop1]
            if (nrow(sumsat1) == 1) rownames(sumsat1) = rn
            drop2 <- rn[which(abs(sumsat1[, 3]) == min(abs(sumsat1[, 3])))]
            sumsat2 = sumsat[rownames(sumsat) != drop2,]
            res <- list(sumsat = sumsat2, drop1 = drop2)
            res
        }

        drops <- function(mod1, X, data) {
          np <- nrow(coef(summary(mod1)))
          if (np == 2) {
            print("np == 2 in the drops was run")
          }
          else if (np > 2){
            print("np > 2 in the drops was run")
            
            sumsat <- coef(summary(mod1))[2:np,]
            
            drop.this <- as.integer(which(abs(sumsat[, 3]) == min(abs(sumsat[, 3]))))
            print(paste("print drop.this which is type -", typeof(drop.this), "-", " with length of ", length(drop.this), sep = ""))
            print(drop.this)
            
            form <- as.formula(mod1$terms)
            print(paste("print form which is type -", typeof(form), "-", " with length of ", length(form), sep = ""))
            print(form)
            
            tt <- terms(form)
            
            temp <- drop.terms(tt, drop.this, keep.response = TRUE)
            
            temp <- reformulate(attr(temp, "term.labels"))
            
            mod.terms <- attr(terms(temp), "term.labels")
            
            f1 <- as.formula(paste(c(paste(Yname, "~ 1 "), mod.terms[mod.terms != drop.this]), collapse = " + "))
            
            new.mod <- lm(f1, data = data)
            
            res <- list(formula = f1, model = new.mod, drop = rownames(sumsat)[which(abs(sumsat[, 3]) == min(abs(sumsat[, 3])))],
                        inf.crit = round(inf(new.mod, inform), digits = digits))
            print(paste("print res which is type -", typeof(res), "-", " with length of ", length(res), sep = ""))
            print(res)
          }
          print("drop finished")
          res
        }

        j = 2;
        temp <- sat
        tab[1,] <- c(paste(c(paste(Yname, "~ 1"), Xall), collapse = " + "), " ", round(inf(temp, inform), digits = digits))
        while (j <= 2) {
            temp <- drops(temp, X, data = data)
            tab[j,] <- c(temp$formula, temp$drop, temp$inf.crit)
            temp <- temp$model
            
            print(paste("attempt # ", j, " was successful", sep = ""))
            
            j = j + 1
        }
    }
    
    print("Addition of the model and AIC value to the first row of tab")
    print(head(tab))
    
    if (inform == "AIC" | inform == "BIC") opt <- which(as.numeric(tab[, 3]) == min(as.numeric(tab[, 3])))
    if (inform == "PRESS") opt <- which(as.numeric(tab[, 3]) == max(as.numeric(tab[, 3])))
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

readFile <- "Datasets/fbData.csv"
varData <- read.csv(file = readFile, header = FALSE, skip = 1)
varData$V2 = NULL
varData$V3 = NULL
varData$V4 = NULL
varData$V5 = NULL
varData$V6 = NULL
varData$V7 = NULL
names(varData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")

#print("print varData")
#typeof(varData)
#print(head(varData))

#greedyData <- greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = varData)[3]
#greedyData <- greedy(Y ~ 1, data = varData)[3]
greedyData <- greedy(Y ~ V1 + V2, data = varData)[3]

#print("greedyData")
#greedyData