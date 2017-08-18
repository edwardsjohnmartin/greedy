\name{ParkinsonsData}
\alias{}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
  This dataset is composed of a range of biomedical voice measurements from 42
  people with early-stage Parkinson's disease recruited to a six-month trial of
  a telemonitoring device for remote symptom progression monitoring. The
  recordings were automatically captured in the patient's homes.
}
\usage{
  data("parkinsons")
}

\format{
  A data frame with 5875 observations on the following 22 variables.
  \describe{
    \item{\code{Subject}}{Integer that uniquely identifies each subject, a numeric vector}
    \item{\code{Age}}{Subject age, a numeric vector}
    \item{\code{Sex}}{Subject gender '0' - male, '1' - female, a numeric vector}
    \item{\code{Test Time}}{Time since recruitment into the trial. The integer part is the 
                            number of days since recruitment, a numeric vector}
    \item{\code{Motor UPDRS}}{Clinician's motor UPDRS score, linearly interpolated, a numeric vector}
    \item{\code{Total UPDRS}}{Clinician's total UPDRS score, linearly interpolated, a numeric vector}
    \item{\code{Jitter, Jitter.Abs, Jitter.RAP, Jitter.PPQ5, Jitter.DDP}}{Several measures of variation in fundamental frequency, a numeric vector}
    \item{\code{Shimmer, Shimmer.dB, Shimmer.APQ3, Shimmer.APQ5, Shimmer.APQ11, Shimmer.DDA}}{Several measures of variation in amplitude, a numeric vector}
    \item{\code{NHR, HNR}}{Two measures of ratio of noise to tonal components in the voice, a numeric vector}
    \item{\code{RPDE}}{A nonlinear dynamical complexity measure, a numeric vector}
    \item{\code{DFA}}{Signal fractal scaling exponent, a numeric vector}
    \item{\code{PPE}}{A nonlinear measure of fundamental frequency variation, a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
Further details are contained in the following reference -- if you use this 
dataset, please cite:
Athanasios Tsanas, Max A. Little, Patrick E. McSharry, Lorraine O. Ramig (2009),
'Accurate telemonitoring of Parkinson.s disease progression by non-invasive 
speech tests',
IEEE Transactions on Biomedical Engineering (to appear).

Further details about the biomedical voice measures can be found in:
Max A. Little, Patrick E. McSharry, Eric J. Hunter, Lorraine O. Ramig (2009),
'Suitability of dysphonia measurements for telemonitoring of Parkinson's 
disease',
IEEE Transactions on Biomedical Engineering, 56(4):1015-1022 
}
\references{
%%  ~~ possibly secondary sources and usages ~~
Little MA, McSharry PE, Hunter EJ, Ramig LO (2009),
'Suitability of dysphonia measurements for telemonitoring of Parkinson's 
disease',
IEEE Transactions on Biomedical Engineering, 56(4):1015-1022

Little MA, McSharry PE, Roberts SJ, Costello DAE, Moroz IM.
'Exploiting Nonlinear Recurrence and Fractal Scaling Properties for Voice 
Disorder Detection',
BioMedical Engineering OnLine 2007, 6:23 (26 June 2007) 

If you use this dataset, please cite the following paper:
A Tsanas, MA Little, PE McSharry, LO Ramig (2009)
'Accurate telemonitoring of Parkinson.s disease progression by non-invasive 
speech tests',
IEEE Transactions on Biomedical Engineering (to appear). 
}
\examples{
data(ParkinsonsData)
## maybe str(ParkinsonsData) ; plot(ParkinsonsData) ...
}
\keyword{datasets}

