\name{PowerPlantData}
\alias{PowerPlantData}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The data was collected from a Combined Cycle Power Plant over 6 years (2006-2011), when the power plant was set to work with full load. 
A combined cycle power plant (CCPP) is composed of gas turbines (GT), steam turbines (ST) and heat recovery steam generators. 
}
\usage{data("PowerPlantData")}

\usage{data("PowerPlantData")}
\format{
  The dataset contains 9568 data points, on the following 5 variables. 
  \describe{
    Input Variables
    \item{\code{Ambient Temperature (AT)}}{a numeric vector}
    \item{\code{Exhaust Vacuum (V)}}{a numeric vector}
    \item{\code{Ambient Pressure (AP)}}{a numeric vector}
    \item{\code{Relative Humidity (RH)}}{a numeric vector}
    Output Variable
    \item{\code{Energy Output (EP)}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Pinar Tüfekci, Prediction of full load electrical power output of a base load operated combined cycle power plant using machine learning methods, International Journal 
of Electrical Power & Energy Systems, Volume 60, September 2014, Pages 126-140, ISSN 0142-0615, http://dx.doi.org/10.1016/j.ijepes.2014.02.027.
(http://www.sciencedirect.com/science/article/pii/S0142061514000908)

Heysem Kaya, 
Pinar Tüfekci , Sadik Fikret Gürgen: Local and Global Learning Methods for Predicting Power of a Combined Gas & Steam Turbine, Proceedings of the International Conference on Emerging Trends in Computer and Electronics Engineering 
ICETCEE 2012, pp. 13-18 (Mar. 2012, Dubai)

}
\examples{
data(PowerPlantData)
## maybe str(PowerPlantData) ; plot(PowerPlantData) ...
}
\keyword{datasets}

