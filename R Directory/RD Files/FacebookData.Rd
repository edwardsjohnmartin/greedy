\name{FacebookData}
\alias{FacebookData}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Facebook Performance Metrics
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The data is related to posts' published during the year of 2014 on the Facebook's page of a renowned cosmetics brand.
   This dataset contains 500 of the 790 rows and part of the features analyzed by Moro et al. (2016). The remaining were 
   omitted due to confidentiality issues.
}
\usage{data("FacebookData")}

\usage{data("FacebookData")}
\format{
  A data frame with 500 observations on the following 19 variables.
  \describe{
    \item{\code{Page Total Likes}}{a numeric vector}
    \item{\code{Type}}{a factor with levels \code{Link}, \code{Photo}, \code{Status}, \code{Video}}
    \item{\code{Category}}{a numeric vector}
    \item{\code{Post Month}}{a numeric vector}
    \item{\code{Post Weekday}}{a numeric vector}
    \item{\code{Post Hour}}{a numeric vector}
    \item{\code{Paid}}{a numeric vector}
    \item{\code{Lifetime Post Total Reach}}{a numeric vector}
    \item{\code{Lifetime Post Total Impressions}}{a numeric vector}
    \item{\code{Lifetime Engaged Users}}{a numeric vector}
    \item{\code{Lifetime Post Consumers}}{a numeric vector}
    \item{\code{Lifetime Post Consumptions}}{a numeric vector}
    \item{\code{Lifetime Post Impressions by People Who Have Liked Your Page}}{a numeric vector}
    \item{\code{Lifetime Post Reach by People Who Like Your Page}}{a numeric vector}
    \item{\code{Lifetime People Who Have Liked Your Page and Engaged With Your Post}}{a numeric vector}
    \item{\code{Comment}}{a numeric vector}
    \item{\code{Like}}{a numeric vector}
    \item{\code{Share}}{a numeric vector}
    \item{\code{Total Interactions}}{a numeric vector}
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
   This dataset is public available for research. The details are described in (Moro et al., 2016).
   Please include this citation if you plan to use this database: 

   (Moro et al., 2016) S. Moro, P. Rita and B. Vala. Predicting social media performance metrics and evaluation 
   of the impact on brand building: A data mining approach. Journal of Business Research, Elsevier, In press.
 
   Available at: http://dx.doi.org/10.1016/j.jbusres.2016.02.010
   The full dataset was described and analyzed in:

   S. Moro, P. Rita and B. Vala. Predicting social media performance metrics and evaluation of the impact on 
   brand building: A data mining approach. Journal of Business Research, Elsevier, In press, Available online 
   since 28 February 2016.
}
\examples{
data(FacebookData)
## maybe str(FacebookData) ; plot(FacebookData) ...
}
\keyword{datasets}

