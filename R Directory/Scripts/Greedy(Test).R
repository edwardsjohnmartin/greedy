printSpacing <- function(){
  print("")
  print("*************************************************************************************************")
  print("")
}

greedy <- function(formula, data = NULL, center = NULL, digits = 5, inform = "AIC") {
    print("start of greedy function")

    # This is the passed in data, I'm not sure what the get_all_vars method does because the data looks the exact same as before it was run
    data <- get_all_vars(formula, data = data)
    print(paste("print data which is type -", typeof(data), "-", " with length of ", length(data), sep = ""))
    print(head(data))

    printSpacing()
    
    # I haven't seen this get run so I'm not sure what it does
    if (!is.null(center)) {
        w <- which(names(data) == center)
        temp <- apply(data[, w], 2, function(x) x - mean(x))
        data[, w] <- temp
        print("if statement was run")
    }

    # This just essentially injects the formula into the data
    m <- model.frame(formula, data = data)
    print(paste("print m which is type -", typeof(m), "-", " with length of ", length(m), sep = ""))
    print(head(m))

    printSpacing()
    
    # This is an array of integers that is a list of all the entries in the response column
    Y <- model.extract(m, "response")
    #print("print Y")
    #print(typeof(Y))
    #print(Y)

    # Not sure what this does yet
    terms <- terms(m)
    print(paste("print terms which is type -", typeof(terms), "-", " with length of ", length(terms), sep = ""))
    print(terms)
    
    printSpacing()
    
    # X is just an array of all the predictor variable names
    X <- attr(terms, "term.labels")
    print(paste("print X which is type -", typeof(X), "-", " with length of ", length(X), sep = ""))
    print(X)

    printSpacing()
    
    # k is a numeric variable to hold the length of X or the total number of predictor variables
    k <- length(X)

    # steps is a numeric variable which holds a calculation on the length of X, will be used to create the rows of the tab matrix
    steps <- 1 + (k ^ 2 + 3 * k) / 2
    print(paste("print steps which is type -", typeof(steps), "-", " with value of ", steps, sep = ""))

    printSpacing()
    
    # function that takes in a linear regression model and applies a formula to it and returns the result as a double
    # our greedy function defines inform to be AIC, so its the only formula used 
    inf <- function(model, inform) {
        switch(inform,
        AIC = AIC(model),
        BIC = BIC(model),
        PRESS = PRESS(model, as.R2 = TRUE))
    }

    # my own created function that take in a model and attempt to do the AIC formula on it
    myAIC <- function(model) {
        -2 * as.numeric(logLik(model)) + 2 * (length(model$coefficients) + 1)
    }

    # creates a matrix table with <steps> rows and 3 columns
    tab <- matrix(nrow = steps, ncol = 3)
    print(paste("print tab which is type -", typeof(tab), "-", " with size of ", steps, "x", 3, sep = ""))
    #print(tab)

    printSpacing()
    
    # Sets the column labels of tab
    colnames(tab) <- c("Model", "Drop", inform)
    print("tab with updated column headers")
    print(head(tab))
    
    printSpacing()
    
    # lm performs least square regression using the formula and data passed into the greedy function
    test <- lm(formula, data = data)
    print(paste("print test which is type -", typeof(test), "-", " with length of ", length(test), sep = ""))
    print(test)

    printSpacing()
    
    # print all the hidden values within test, attributes(test) provides a list with the information you can retrieve from it
    #print(attributes(test))
    #for (j in 1:length(test)) {
    #    print(paste("test[", j, "]", sep = ""))
    #    print(test[j])
    #}

    # If steps == 1 (only way for this to happen is if there are no predictor variables)
    if (steps == 1) {
        # Fill the table matrix with the formula model being used, nothing for drop, and the AIC value rounded to 5 digits
        tab[1,] = c(deparse(formula(test$terms)), " ", round(inf(test, inform), digits = digits))
        #tab[1,] = c(deparse(formula(test$terms)), " ", round(myAIC(test), digits = digits))
        print(tab)
    }
    # If steps >= 2 (this will happen if there are 1 or more predictor variables)
    else if (steps >= 2) {
      # puts the predictor variable names from terms into an numeric array of the same length as pred. variables, also puts "numeric" in each element
      d <- attr(terms, "dataClasses")[2:(length(X) + 1)]
      print(paste("print d which is type -", typeof(d), "-", " with length of ", length(d), sep = ""))
      print(d)
      
      printSpacing()
      
      # Seems to create the same array of all the predictor variable names as what the variable X is
      # the d == "numeric" part in the arguement may make it function different, but I haven't discovered the extra functionality yet
      Xn <- X[d == "numeric"]
      print(paste("print Xn which is type -", typeof(Xn), "-", " with length of ", length(Xn), sep = ""))
      print(Xn)
      
      printSpacing()
      
      # Takes each element in the Xn character array and adds "I(" to the front and "^2)" to the end
      Xsq <- paste("I(", Xn, "^2)", sep = "")
      print(paste("print Xsq which is type -", typeof(Xsq), "-", " with length of ", length(Xsq), sep = ""))
      print(Xsq)
      
      printSpacing()
      
      # Xint is a character array of length [len(X), len(X)] that puts "X[x]:X[x]" into each element of the array
      Xint <- outer(X, X, function(x, y) paste(x, ":", y, sep = ""))
      print(paste("print xint which is type -", typeof(Xint), "-", " with length of ", length(Xint), sep = ""))
      print(Xint)
      
      printSpacing()
      
      # Creates a new array of the upper diagonal of the Xint array
      Xint <- Xint[upper.tri(Xint)]
      print(paste("print xint which is type -", typeof(Xint), "-", " with length of ", length(Xint), sep = ""))
      print(Xint)
      
      printSpacing()
      
      # Takes evertything that was in the X variable, Xint variable, or Xsq variable and and puts them all into the same array
      Xall <- c(X, Xint, Xsq)
      print(paste("print xall which is type -", typeof(Xall), "-", " with length of ", length(Xall), sep = ""))
      print(Xall)
      
      printSpacing()
      
      #
      if (!any(match(Xn, X))){ 
        printSpacing()
        Xall <- c(X, Xint)
        print("Xall was changed")
      }
      
      # Yname stores the name of the predictor variable column
      Yname <- names(m)[1]
      print(paste("print Yname which is type -", typeof(Yname), "-", " with length of ", length(Yname), sep = ""))
      print(Yname)
      
      printSpacing()
      
      # f is a formula variable that creates a formula with the form -> Yname ~ 1 + (all the items in Xall)
      f <- as.formula(paste(c(paste(Yname, "~ 1 "), Xall), collapse = " + "))
      print(paste("print f which is type -", typeof(f), "-", " with length of ", length(f), sep = ""))
      print(f)
      
      printSpacing()
      
      # sat is a linear fit model obtained by using the formula f on the data passed into Greedy
      sat <- lm(f, data = data)
      print(paste("print sat which is type -", typeof(sat), "-", " with length of ", length(sat), sep = ""))
      print(sat)
      
      printSpacing()
      
      #
      redo <- function(drop1, sumsat) {
        sumsat1 <- sumsat[rownames(sumsat) != drop1,] 
      }
      
      #
      drops <- function(mod1, X, data) {
        # np is the number of rows in the coefficient matrix of the lm model passed into the function
        np <- nrow(coef(summary(mod1)))
        print(paste("print np which is type -", typeof(np), "-", " with value of ", np, sep = ""))
        
        printSpacing()
        
        # I haven't figured out how to get this part to run yet
        if (np == 2) {
          print("np == 2 in the drops was run")
          
          #
          new.mod <- update(mod1, ~ 1)
          print(paste("print new.mod which is type -", typeof(new.mod), "-", " with length of ", length(new.mod), sep = ""))
          print(new.mod)
          
          printSpacing()
          
          #
          res <- list(formula = paste(Yname, "~ 1"), model = new.mod, drop = attr(terms(mod1), "term.labels"),
                      inf.crit = round(inf(new.mod, inform), digits = digits))
        }
        #
        else if (np > 2){
          print("np > 2 in the drops was run")
          
          # This is a summary of the coeffiecients based on the model passed into drops
          # It includes the Estimate, Std. Error, t value, and Pr(>|t|)
          sumsat <- coef(summary(mod1))[2:np,]
          print(paste("print sumsat which is type -", typeof(sumsat), "-", " with length of ", length(sumsat), sep = ""))
          print(sumsat)
          
          printSpacing()
          
          # drop1 is the name of the variable with the lowest magnitude t value
          drop1 <- rownames(sumsat)[which(abs(sumsat[, 3]) == min(abs(sumsat[, 3])))]
          print(paste("print drop1 which is type -", typeof(drop1), "-", " with length of ", length(drop1), sep = ""))
          print(drop1)
          
          printSpacing()
          
          # mod.terms is a list of all the variables in the model passed into drops
          mod.terms <- attr(terms(mod1), "term.labels")
          print(paste("print mod.terms which is type -", typeof(mod.terms), "-", " with length of ", length(mod.terms), sep = ""))
          print(mod.terms)
        
          printSpacing()
      
          # if drop1 cannot be found in the mod.terms list then this will evaluate to true, i believe its used as error checking
          if (any(X == drop1) & length(grep(drop1, mod.terms[mod.terms != drop1])) > 0) {
            print("if (any(X == drop1) & length(grep(drop1, mod.terms[mod.terms != drop1])) > 0) evaluated to true")
            
            printSpacing()
            
            # haven't checked this part yet because i haven't got the if statement to evaluate to true yet
            drop1 <- redo(drop1, sumsat)$drop1
            print(paste("print drop1 which is type -", typeof(drop1), "-", " with length of ", length(drop1), sep = ""))
            print(drop1)
          }
          
          # f1 is a new formula made from mod.terms disincluding drop1
          f1 <- as.formula(paste(c(paste(Yname, "~ 1 "), mod.terms[mod.terms != drop1]), collapse = " + "))
          print(paste("print f1 which is type -", typeof(f1), "-", " with length of ", length(f1), sep = ""))
          print(f1)
          
          printSpacing()
          
          # new.mod is a linear model with the new f1 formula
          new.mod <- lm(f1, data = data)
          print(paste("print new.mod which is type -", typeof(new.mod), "-", " with length of ", length(new.mod), sep = ""))
          print(new.mod)
          
          printSpacing()
          
          # res is a list with length 4 with the new formula, new linear model, the term that was dropped, and the AIC value of the new formula
          res <- list(formula = paste(c(paste(Yname, "~ 1 "), mod.terms[mod.terms != drop1]), collapse = " + "),
                      model = new.mod, drop = drop1, inf.crit = round(inf(new.mod, inform), digits = digits))
          print(paste("print res which is type -", typeof(res), "-", " with length of ", length(res), sep = ""))
          print(res)
          
          printSpacing()
        }
        res
      }
      
      j = 2;
      
      # lmTemp is equal to our variable that holds the lm model
      lmTemp <- sat
      
      # Replaces the 1st row of the tab matrix with the entire formula and its AIC value
      tab[1,] <- c(paste(c(paste(Yname, "~ 1"), Xall), collapse = " + "), " ", round(inf(lmTemp, inform), digits = digits))
      print("Addition of the model and AIC value to the first row of tab")
      print(head(tab))
      
      printSpacing()
      
      # This loop will from 2 to steps times where steps = 1 + (k ^ 2 + 3 * k) / 2 where k is the number of predictor variables
      # ex if k = 8 => steps = 45, if k = 11 => steps = 78, etc...
      while (j <= 2) {
        
        # temp is a list with length 4 with the new formula, new linear model, the term that was dropped, and the AIC value of the new formula
        temp <- drops(lmTemp, X, data = data)
        print(paste("print temp which is type -", typeof(temp), "-", " with length of ", length(temp), sep = ""))
        print(temp)
        
        printSpacing()
        
        # Adds the new formaula, the term that was dropped, and the new formula's AIC value
        tab[j,] <- c(temp$formula, temp$drop, temp$inf.crit)
        print("Addition of the model and AIC value to the first row of tab")
        print(head(tab))
        
        printSpacing()
        
        # temp becomes just the linear model to be passed into drops again
        temp <- temp$model
        print(paste("print temp which is type -", typeof(temp), "-", " with length of ", length(temp), sep = ""))
        print(temp)
        
        printSpacing()
        
        j = j + 1
      }
    }

    print("end of greedy function")
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
greedyData <- greedy(Y ~ 1 + V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V11, data = varData)[3]

#print("greedyData")
#greedyData