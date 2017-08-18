\name{ConcreteData}
\alias{ConcreteData}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
Concrete Compressive Strength -
The actual concrete compressive strength (MPa) for a given mixture under a 
specific age (days) was determined from laboratory. Data is in raw form (not scaled). 
}
\usage{data("ConcreteData")}

\usage{data("ConcreteData")}
\format{
  A data frame with 1030 observations on the following 9 variables.
  \describe{
    Input Variables
    \item{\code{Cement}}{kg in a m3 mixture, a numeric vector}
    \item{\code{Blast Furnace Slag}}{kg in a m3 mixture, a numeric vector}
    \item{\code{Fly Ash}}{kg in a m3 mixture, a numeric vector}
    \item{\code{Water}}{kg in a m3 mixture, a numeric vector}
    \item{\code{Superplasticizer}}{kg in a m3 mixture, a numeric vector}
    \item{\code{Coarse Aggregate}}{kg in a m3 mixture, a numeric vector}
    \item{\code{Fine Aggregate}}{kg in a m3 mixture, a numeric vector}
    \item{\code{Age}}{day(1~365), a numeric vector}
    Output Variable
    \item{\code{Concrete compressive strength}}{MPa, a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
Given is the variable name, variable type, the measurement unit and a brief description. 
The concrete compressive strength is the regression problem. The order of this listing 
corresponds to the order of numerals along the rows of the database. 
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Original Owner and Donor
  Prof. I-Cheng Yeh
  Department of Information Management 
  Chung-Hua University, 
  Hsin Chu, Taiwan 30067, R.O.C.
  e-mail:icyeh@chu.edu.tw
  TEL:886-3-5186511
}
\references{
%%  ~~ possibly secondary sources and usages ~~

Past Usage: 

Main
1. I-Cheng Yeh, "Modeling of strength of high performance concrete using artificial 
neural networks," Cement and Concrete Research, Vol. 28, No. 12, pp. 1797-1808 (1998).

Others
2. I-Cheng Yeh, "Modeling Concrete Strength with Augment-Neuron Networks," J. of 
Materials in Civil Engineering, ASCE, Vol. 10, No. 4, pp. 263-268 (1998).

3. I-Cheng Yeh, "Design of High Performance Concrete Mixture Using Neural Networks,"  
J. of Computing in Civil Engineering, ASCE, Vol. 13, No. 1, pp. 36-42 (1999).

4. I-Cheng Yeh, "Prediction of Strength of Fly Ash and Slag Concrete By The Use of 
Artificial Neural Networks," Journal of the Chinese Institute of Civil and Hydraulic 
Engineering, Vol. 15, No. 4, pp. 659-663 (2003).

5. I-Cheng Yeh, "A mix Proportioning Methodology for Fly Ash and Slag Concrete Using 
Artificial Neural Networks," Chung Hua Journal of Science and Engineering, Vol. 1, No. 
1, pp. 77-84 (2003).

6. Yeh, I-Cheng, "Analysis of strength of concrete using design of experiments and 
neural networks,": Journal of Materials in Civil Engineering, ASCE, Vol.18, No.4, 
pp.597-604 ?2006?.

Acknowledgements, Copyright Information, and Availability:

NOTE: Reuse of this database is unlimited with retention of copyright notice for 
Prof. I-Cheng Yeh and the following published paper:

I-Cheng Yeh, "Modeling of strength of high performance concrete using artificial 
neural networks," Cement and Concrete Research, Vol. 28, No. 12, pp. 1797-1808 (1998)
}
\examples{
data(ConcreteData)
## maybe str(ConcreteData) ; plot(ConcreteData) ...
}
\keyword{datasets}

