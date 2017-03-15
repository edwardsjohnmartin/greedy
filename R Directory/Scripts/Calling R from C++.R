set.seed(42)
x <- rnorm(1e5)

typeof(x)

fivenum(x)

library(Rcpp)
sourceCpp("C++ files/RfromC.cpp")

callFunction(x, fivenum)

DoubleMe(5)

DoubleMe(20)

y <- 100

DoubleMe(y)