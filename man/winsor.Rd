\name{winsor}
\alias{winsor}
\title{Winsor normalization with fitting to <0,1> range.}
\usage{
winsor(input_vector, winsor_level = 0.9, only_top = FALSE)
}
\arguments{
\item{input_vector}{Vector with values to be Winsorized}

\item{winsor_level}{Winsorization level. Bottom outliers will be set to
(1-winsor_level)/2 quantile and top outliers to (1+winsor_level)/2 quantile.}

\item{only_top}{If TRUE then bottom values are not Winsorized and the lowest
is set to 0.}
}
\value{
Vector of numerics within <0,1>.
}
\description{
Function performs Winsor normalization of a supplied vector. Steps:
1. Calcualate top winsor value [(1+winsor_level)/2 quantile], and bottom
winsor value ((1-winsor_level)/2 quantile)
2. Each value below bottom winsor value set to bottom winsor value; each
value above top winsor value set to top winsor value
3. Transform linearly all the values to [0,1] range
}
\examples{

data_set <- runif(1:100)*100
plot(winsor(data_set, winsor_level=0.8) ~ data_set)

}
\references{
Hastings, Cecil; Mosteller, Frederick; Tukey, John W.; Winsor,
Charles P. Low Moments for Small Samples: A Comparative Study of Order
Statistics. The Annals of Mathematical Statistics 18 (1947), no. 3,
413--426.
}
\author{
Lukasz Jan Kielpinski
}
\keyword{~winsorising}
