\name{WhiteWineData}
\alias{WhiteWineData}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
White Wine Quality
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The data is related to white variants of the Portuguese "Vinho Verde" wine.
The inputs include objective tests (e.g. PH values) and the output is based on sensory data
  (median of at least 3 evaluations made by wine experts). Each expert graded the wine quality 
  between 0 (very bad) and 10 (very excellent).
The dataset can be viewed as classification or regression tasks.
}
\usage{data("WhiteWineData")}

\usage{data("WhiteWineData")}
\format{
  A data frame with 4898 observations on the following 12 variables.
  \describe{
    Output Variable
    \item{\code{Quality}}{a numeric vector}
    Input Variables
    \item{\samp{Volatile Acidity}}{a numeric vector}
    \item{\samp{Citric Acid}}{a numeric vector}
    \item{\samp{Residual Sugar}}{a numeric vector}
    \item{\code{Chlorides}}{a numeric vector}
    \item{\samp{Free Sulfur Dioxide}}{a numeric vector}
    \item{\samp{Total Sulfur Dioxide}}{a numeric vector}
    \item{\code{Density}}{a numeric vector}
    \item{\code{PH}}{a numeric vector}
    \item{\code{Sulphates}}{a numeric vector}
    \item{\code{Alcohol}}{a numeric vector}
    \item{\samp{Fixed Acidity}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
  Created by: Paulo Cortez (Univ. Minho), Antonio Cerdeira, Fernando Almeida, Telmo Matos and Jose Reis (CVRVV) @ 2009
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Past Usage:

  P. Cortez, A. Cerdeira, F. Almeida, T. Matos and J. Reis. 
  Modeling wine preferences by data mining from physicochemical properties.
  In Decision Support Systems, Elsevier, 47(4):547-553. ISSN: 0167-9236.
}
\examples{
data(WhiteWineData)
## maybe str(WhiteWineData) ; plot(WhiteWineData) ...
}
\keyword{datasets}

