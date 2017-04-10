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
                new.mod <- update(mod1, ~ 1)
                res <- list(formula = paste(Yname, "~ 1"), model = new.mod, drop = attr(terms(mod1), "term.labels"),
                      inf.crit = round(inf(new.mod, inform), digits = digits))
            }
            if (np > 2) {
                sumsat <- coef(summary(mod1))[2:np,]
                drop1 <- rownames(sumsat)[which(abs(sumsat[, 3]) == min(abs(sumsat[, 3])))]
                mod.terms <- attr(terms(mod1), "term.labels")

                if (any(X == drop1) & length(grep(drop1, mod.terms[mod.terms != drop1])) > 0) {
                    drop1 <- redo(drop1, sumsat)$drop1
                    if (any(X == drop1) & length(grep(drop1, mod.terms[mod.terms != drop1]) > 0)) {
                        drop1 <- redo(drop1, sumsat)$drop1
                        if (any(X == drop1) & length(grep(drop1, mod.terms[mod.terms != drop1]) > 0)) {
                            drop1 <- redo(drop1, sumsat)$drop1
                        }
                    }
                }
                f1 <- as.formula(paste(c(paste(Yname, "~ 1 "), mod.terms[mod.terms != drop1]), collapse = " + "))
                new.mod <- lm(f1, data = data)
                res <- list(formula = paste(c(paste(Yname, "~ 1 "), mod.terms[mod.terms != drop1]), collapse = " + "),
                      model = new.mod, drop = drop1, inf.crit = round(inf(new.mod, inform), digits = digits))
            }
            res
        }

        j = 2;
        temp <- sat
        tab[1,] <- c(paste(c(paste(Yname, "~ 1"), Xall), collapse = " + "), " ", round(inf(temp, inform), digits = digits))
        while (j <= steps) {
            temp <- drops(temp, X, data = data)
            tab[j,] <- c(temp$formula, temp$drop, temp$inf.crit)
            temp <- temp$model
            j = j + 1
        }
        print(tab)
    }
    if (inform == "AIC" | inform == "BIC") opt <- which(as.numeric(tab[, 3]) == min(as.numeric(tab[, 3])))
    if (inform == "PRESS") opt <- which(as.numeric(tab[, 3]) == max(as.numeric(tab[, 3])))
    
    #---------added for testing----------
    print(paste("print opt which is type -", typeof(opt), "-", " with length of ", length(opt), sep = ""))
    print(opt)
    #---------added for testing----------
    
    best <- tab[opt,][1]
    #---------added for testing----------
    print(paste("print best which is type -", typeof(best), "-", " with length of ", length(best), sep = ""))
    print(best)
    #---------added for testing----------
    
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
names(varData) <- c("Y", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
print(head(varData))

greedyData <- greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = varData)[3]

print("greedyData")
greedyData