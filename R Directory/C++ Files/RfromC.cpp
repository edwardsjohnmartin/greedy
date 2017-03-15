#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector callFunction(NumericVector x, Function f) {
    NumericVector res = f(x);
    return res;
}

// [[Rcpp::export]]
int DoubleMe(int y){
	return y * 2;
}	

// [[Rcpp::export]]
int AddMe(int y){
	return y + 2;
}