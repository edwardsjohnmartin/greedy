# I have taken all the formatting code out of the RMD file and only included where the various functions and variables get defined
# I have also only included the case1202 file in the script and excluded the case0902 one for now
# Running everything in this script will pull the data from the csv and will output the plot at the end
# No matter what I do though, one error persists. After running all the code, the message
# "Error: object 'k' not found"
# always comes up. I have not been able to figure out why or how to fix it. 

# Notes for 2/27/2017
# Test running time against other packages such as SAS(JMP PRO), -> Use Proc Reg(fits model to the data)
# Maybe use SPSS if there is time



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
    }
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
    print(paste("print res which is type -", typeof(res), "-", " with length of ", length(res), sep = ""))
    print(res)
    res
}

print.greedy <- function(x, ...) {
    cat("\n")
    out <- structure(x$out)
    print(out)
    invisible(x)
}

library(MASS)
library(car)

# This is where you put the location of the csv file you will be examining
readFile <- "Datasets/fbData.csv"
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

len <<- length(varData)
print("Preview of the dataset being used")
head(varData)

print("How many columns are in the dataset")
len

#names(varData) <- c("Y", "Xsal", "Xsex", "Xsen", "Xage", "Xed", "Xexp", "x1", "x2", "x3", "x4", "x5")
#g <- greedy(Y ~ Xsen + Xage + Xed + Xexp, data = varData)
#g <- greedy(Y ~ V1 + V2 + V3 + V4, data = varData)
#print("line 171 'g'")
#g
#print("line 173 'g$best'")
#g$best
#print("line 175 'vif(g$best)'")
#vif(g$best)
#print("line 177 'stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3 + V4) ^ 2, data = varData), trace = FALSE)'")
#stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3 + V4) ^ 2, data = varData), trace = FALSE)

times4 <- function() {
    require(MASS)
    # Commented out these lines presumably so that it will use the data directly in the spreadsheet as opposed to randomly created numbers
    sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
    ss0 <- system.time(lm(Y ~ 1, data = varData))[3]

    sg1 <- system.time(greedy(Y ~ V1, data = varData))[3]
    ss1 <- system.time(stepAIC(lm(Y ~ 1 + V1 + I(V1 ^ 2), data = varData), trace = F))[3]

    sg2 <- system.time(greedy(Y ~ V1 + V2, data = varData))[3]
    ss2 <- system.time(stepAIC(lm(Y ~ 1 + V1 + V2 + V1:V2 + I(V1 ^ 2) + I(V2 ^ 2), data = varData), trace = F))[3]

    sg3 <- system.time(greedy(Y ~ V1 + V2 + V3, data = varData))[3]
    ss3 <- system.time(stepAIC(lm(Y ~ 1 + I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + (V1 + V2 + V3) ^ 2, data = varData), trace = F))[3]

    #This is where I'm going to store the data for checking to make sure the models returned by each algorithm is the same
    greedyCoef <<- greedy(Y ~ V1 + V2 + V3, data = varData)
    stepAICCoef <<- stepAIC(lm(Y ~ 1 + I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + (V1 + V2 + V3) ^ 2, data = varData), trace = F)

    sg <- c(sg0, sg1, sg2, sg3)
    ss <- c(ss0, ss1, sg2, ss3)
    list(sg = sg, ss = ss)
}

times5 <- function() {
    require(MASS)
    # Commented out these lines presumably so that it will use the data directly in the spreadsheet as opposed to randomly created numbers
    sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
    ss0 <- system.time(lm(Y ~ 1, data = varData))[3]

    sg1 <- system.time(greedy(Y ~ V1, data = varData))[3]
    ss1 <- system.time(stepAIC(lm(Y ~ 1 + V1 + I(V1 ^ 2), data = varData), trace = F))[3]

    sg2 <- system.time(greedy(Y ~ V1 + V2, data = varData))[3]
    ss2 <- system.time(stepAIC(lm(Y ~ 1 + V1 + V2 + V1:V2 + I(V1 ^ 2) + I(V2 ^ 2), data = varData), trace = F))[3]

    sg3 <- system.time(greedy(Y ~ V1 + V2 + V3, data = varData))[3]
    ss3 <- system.time(stepAIC(lm(Y ~ 1 + I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3) ^ 2, data = varData), trace = F))[3]

    sg4 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4, data = varData))[3]
    ss4 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3 + V4) ^ 2, data = varData), trace = F))[3]

    #This is where I'm going to store the data for checking to make sure the models returned by each algorithm is the same
    greedyCoef <<- greedy(Y ~ V1 + V2 + V3 + V4, data = varData)
    stepAICCoef <<- stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3 + V4) ^ 2, data = varData), trace = F)

    sg <- c(sg0, sg1, sg2, sg3, sg4)
    ss <- c(ss0, ss1, sg2, ss3, ss4)
    list(sg = sg, ss = ss)
}

times9 <- function() {
    require(MASS)
    # Commented out these lines presumably so that it will use the data directly in the spreadsheet as opposed to randomly created numbers
    sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
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

    #This is where I'm going to store the data for checking to make sure the models returned by each algorithm is the same
    greedyCoef <<- greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8, data = varData)
    stepAICCoef <<- stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8) ^ 2, data = varData), trace = F)

    sg <- c(sg0, sg1, sg2, sg3, sg4, sg5, sg6, sg7, sg8)
    ss <- c(ss0, ss1, sg2, ss3, ss4, ss5, ss6, ss7, ss8)
    list(sg = sg, ss = ss)
}

times12 <- function() {
    require(MASS)
    sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
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

    #This is where I'm going to store the data for checking to make sure the models returned by each algorithm is the same
    greedyCoef <<- greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = varData)
    stepAICCoef <<- stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11) ^ 2, data = varData), trace = F)

    sg <- c(sg0, sg1, sg2, sg3, sg4, sg5, sg6, sg7, sg8, sg9, sg10, sg11)
    ss <- c(ss0, ss1, sg2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11)
    list(sg = sg, ss = ss)
}

times13 <- function() {
    require(MASS)
    sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
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

    #This is where I'm going to store the data for checking to make sure the models returned by each algorithm is the same
    greedyCoef <<- greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, data = varData)
    stepAICCoef <<- stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12) ^ 2, data = varData), trace = F)

    sg <- c(sg0, sg1, sg2, sg3, sg4, sg5, sg6, sg7, sg8, sg9, sg10, sg11, sg12)
    ss <- c(ss0, ss1, sg2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11, ss12)
    list(sg = sg, ss = ss)
}

times17 <- function() {
    require(MASS)
    print("start sg/ss 0")
    sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
    ss0 <- system.time(lm(Y ~ 1, data = varData))[3]

    print("start sg/ss 1")
    sg1 <- system.time(greedy(Y ~ V1, data = varData))[3]
    ss1 <- system.time(stepAIC(lm(Y ~ 1 + V1 + I(V1 ^ 2), data = varData), trace = F))[3]

    print("start sg/ss 2")
    sg2 <- system.time(greedy(Y ~ V1 + V2, data = varData))[3]
    ss2 <- system.time(stepAIC(lm(Y ~ 1 + V1 + V2 + V1:V2 + I(V1 ^ 2) + I(V2 ^ 2), data = varData), trace = F))[3]

    print("start sg/ss 3")
    sg3 <- system.time(greedy(Y ~ V1 + V2 + V3, data = varData))[3]
    ss3 <- system.time(stepAIC(lm(Y ~ 1 + I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 4")
    sg4 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4, data = varData))[3]
    ss4 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + (V1 + V2 + V3 + V4) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 5")
    sg5 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5, data = varData))[3]
    ss5 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + (V1 + V2 + V3 + V4 + V5) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 6")
    sg6 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6, data = varData))[3]
    ss6 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 7")
    sg7 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7, data = varData))[3]
    ss7 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 8")
    sg8 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8, data = varData))[3]
    ss8 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 9")
    sg9 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9, data = varData))[3]
    ss9 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 10")
    sg10 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10, data = varData))[3]
    ss10 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 11")
    sg11 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = varData))[3]
    ss11 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 12")
    sg12 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12, data = varData))[3]
    ss12 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12) ^ 2, data = varData), trace = F))[3]

    print("start sg 13")
    #sg13 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13, data = varData))[3]
    print("start ss 13")
    ss13 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + I(V13 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 14")
    sg14 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14, data = varData))[3]
    ss14 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + I(V13 ^ 2) + I(V14 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 15")
    sg15 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15, data = varData))[3]
    ss15 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + I(V13 ^ 2) + I(V14 ^ 2) + I(V15 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15) ^ 2, data = varData), trace = F))[3]

    print("start sg/ss 16")
    sg16 <- system.time(greedy(Y ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16, data = varData))[3]
    ss16 <- system.time(stepAIC(lm(Y ~ I(V1 ^ 2) + I(V2 ^ 2) + I(V3 ^ 2) + I(V4 ^ 2) + I(V5 ^ 2) + I(V6 ^ 2) + I(V7 ^ 2) + I(V8 ^ 2) + I(V9 ^ 2) + I(V10 ^ 2) + I(V11 ^ 2) + I(V12 ^ 2) + I(V13 ^ 2) + I(V14 ^ 2) + I(V15 ^ 2) + I(V16 ^ 2) + (V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11 + V12 + V13 + V14 + V15 + V16) ^ 2, data = varData), trace = F))[3]

    sg <- c(sg0, sg1, sg2, sg3, sg4, sg5, sg6, sg7, sg8, sg9, sg10, sg11, sg12, sg13, sg14, sg15, sg16)
    ss <- c(ss0, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9, ss10, ss11, ss12, ss13, ss14, ss15, ss16)
    list(sg = sg, ss = ss)
}

times <- function() {
    require(MASS)
    sg0 <- system.time(greedy(Y ~ 1, data = varData))[3]
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

    sg <- c(sg0, sg1, sg2, sg3, sg4, sg5, sg6, sg7, sg8, sg9)
    ss <- c(ss0, ss1, ss2, ss3, ss4, ss5, ss6, ss7, ss8, ss9)
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

gsteps <- 1 + (k ^ 2 + 3 * k) / 2

sapply(0:len, nmodels)

sim <- 2
#time <- array(dim = c(10, 2, sim), dimnames = list(0:9, c("ss", "sg"), 1:sim))
# Had to change this to be for the general case and not hard-coded in for the original times length
time <- array(dim = c(len, 2, sim), dimnames = list(0:(len-1), c("ss", "sg"), 1:sim))

# This is the only place times gets called
for (i in 1:sim) {
    if (len == 4) {
        ss <- times4()$ss
        sg <- times4()$sg
    }
    else if (len == 5) {
        ss <- times5()$ss
        sg <- times5()$sg
    } else if (len == 9) {
        ss <- times9()$ss
        sg <- times9()$sg
    } else if (len == 12) {
        ss <- times12()$ss
        sg <- times12()$sg
    } else if (len == 13) {
        ss <- times13()$ss
        sg <- times13()$sg
    } else if (len == 17) {
        ss <- times17()$ss
        sg <- times17()$sg
    }
    time[,, i][, 1] <- ss
    time[,, i][, 2] <- sg
}

#This code will produce the output showing if the coeffiecients of each algorithm is the same

print("Show me the coefficients for the greedy algorithm")
print(greedyCoef$best)

print("Show me the coefficients for the stepAIC algorithm")
print(stepAICCoef)

#This code will produce the output showing the times each algorithm takes to evaluate the regression
print('The time it takes to evaluate each regression line using the Greedy algorithm is')
print(sg)
print('The time it takes to evaluate each regression line using the StepAIC algorithm is')
print(ss)

ss <- matrix(ncol = sim, nrow = len)

for (i in 1:sim) {
    ss[, i] <- time[,, i][, 1][1:len]
}

sg <- matrix(ncol = sim, nrow = len)
for (i in 1:sim) {
    sg[, i] <- time[,, i][, 2][1:len]
}

#Here is where the attempt to log transform the times will go
ssLog <- log(ss + 1)
sgLog <- log(sg + 1)

ss.means <- apply(ss, 1, mean)
ss.sds <- apply(ss, 1, sd)

print("print ss.means")
print(ss.means)
print("print ss.sds")
print(ss.sds)

sg.means <- apply(sg, 1, mean)
sg.sds <- apply(sg, 1, sd)

print("print sg.means")
print(sg.means)
print("print sg.sds")
print(sg.sds)

ssLog.means <- apply(ssLog, 1, mean)
ssLog.sds <- apply(ssLog, 1, sd)

print("print ssLog.means")
print(ssLog.means)
print("print ssLog.sds")
print(ssLog.sds)

sgLog.means <- apply(sgLog, 1, mean)
sgLog.sds <- apply(sgLog, 1, sd)

print("print sgLog.means")
print(sgLog.means)
print("print sgLog.sds")
print(sgLog.sds)

nmodels <- function(k) {
    Mk <- matrix(ncol = 1, nrow = length(seq(0:k)))
    for (j in 0:k) {
        Mk[j + 1] <- choose(k, j) * 2 ^ ((j ^ 2 + j) / 2)
    }
    res <- sum(Mk)
    if (k == 0) res = 1
    res
}

print("ss is")
print(ss)

print("ss log transformed is")
print(ssLog)

print("sg is")
print(sg)

print("sg log transformed is")
print(sgLog)

if (ss[len] >= sg[len]) {
    maxTime <- ss[len]
} else {
    maxTime <- sg[len]
}

if (ssLog[len] >= sgLog[len]) {
    logMaxTime <- ssLog[len]
} else {
    logMaxTime <- sgLog[len]
}

print("maxTime equals")
print(maxTime)

print("logMaxTime equals")
print(logMaxTime)

CreatePlotNormalTimes <- function() {
    plot.new()
    tiff(file = "tempNorm.tiff", width = 1600, height = 1600, units = "px", res = 200)
    # All of the lines below controls making the graphs
    par(mar = c(5, 4, 4, 1.2))
    # In the line below, ylim = c(0,3) controls the range of the y-axis. Originally the 3 was a .2 but the range was too small after I added data in the csv
    plot(0:(len - 1), ss.means, type = "l", col = 1, ylab = "Accumulated system time (sec)", xlab = "Number of exp. variables", xlim = c(0, (len - 1)), ylim = c(0, maxTime))
    axis(side = 1, at = 0:(len - 1), cex.axis = 1)
    points(0:(len - 1), sg.means, type = "l", col = 2)
    #segments(0:(len - 1), ss.means - ss.sds, 0:(len - 1), ss.means + ss.sds)
    #segments(0:(len - 1), sg.means - sg.sds, 0:(len - 1), sg.means + sg.sds, col = 2)
    legend("topleft", lty = 1, col = 1:2, legend = c("stepAIC", "Greedy algorithm"))
    #axis(3, at = 0:9, labels = c(1, 3, 13, 95, 1337, 38619, expression(2.3 %*% 10 ^ 6), expression(2.8 %*% 10 ^ 8), expression(7.0 %*% 10 ^ 10), expression(3.6 %*% 10 ^ 13)), cex.axis = .75)
    mtext(side = 3, outer = F, paste("Timing of ", readFile, " by algorithm with normal times"), line = 2)
    dev.off()
}

CreatePlotLogTimes <- function() {
    plot.new()
    tiff(file = "tempLog.tiff", width = 1600, height = 1600, units = "px", res = 200)
    # All of the lines below controls making the graphs
    par(mar = c(5, 4, 4, 1.2))
    # In the line below, ylim = c(0,3) controls the range of the y-axis. Originally the 3 was a .2 but the range was too small after I added data in the csv
    plot(0:(len - 1), ssLog.means, type = "l", col = 1, ylab = "Accumulated system time (sec)", xlab = "Number of exp. variables", xlim = c(0, (len - 1)), ylim = c(0, logMaxTime))
    axis(side = 1, at = 0:(len - 1), cex.axis = 1)
    points(0:(len - 1), sgLog.means, type = "l", col = 2)
    #segments(0:(len - 1), ssLog.means - ssLog.sds, 0:(len - 1), ssLog.means + ssLog.sds)
    #segments(0:(len - 1), sgLog.means - sgLog.sds, 0:(len - 1), sgLog.means + sgLog.sds, col = 2)
    legend("topleft", lty = 1, col = 1:2, legend = c("stepAIC", "Greedy algorithm"))
    #axis(3, at = 0:9, labels = c(1, 3, 13, 95, 1337, 38619, expression(2.3 %*% 10 ^ 6), expression(2.8 %*% 10 ^ 8), expression(7.0 %*% 10 ^ 10), expression(3.6 %*% 10 ^ 13)), cex.axis = .75)
    mtext(side = 3, outer = F, paste("Timing of ", readFile, " by algorithm with y-axis log transformed"), line = 2)
    dev.off()
}

CreatePlotNormalTimes()

#CreatePlotLogTimes()