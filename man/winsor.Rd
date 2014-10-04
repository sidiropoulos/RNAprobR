\name{winsor}
\alias{winsor}
\title{
Winsor normalization with fitting to <0,1> range.
}
\description{
Function performs Winsor normalization of a supplied vector, sets bottom (winsorized) value to 0, top 1 and linearly transforms all the values to the <0,1> range.
}
\usage{
winsor(input_vector, winsor_level = 0.9, only_top = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_vector}{
Vector with values to be Winsorized
}
  \item{winsor_level}{
Winsorization level. Bottom outliers will be set to (1-winsor_level)/2 quantile and top outliers to (1+winsor_level)/2 quantile.
}
  \item{only_top}{
If TRUE then bottom values are not Winsorized and the lowest is set to 0.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Vector of numerics within <0,1>.
}
\references{
Hastings, Cecil; Mosteller, Frederick; Tukey, John W.; Winsor, Charles P. Low Moments for Small Samples: A Comparative Study of Order Statistics. The Annals of Mathematical Statistics 18 (1947), no. 3, 413--426. 
}
\author{
Lukasz Jan Kielpinski
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data_set <- runif(1:100)*100
plot(winsor(data_set, winsor_level=0.8) ~ data_set)


## The function is currently defined as
function (input_vector, winsor_level = 0.9, only_top = F) 
{
    bounds <- quantile(input_vector, c((1 - winsor_level)/2, 
        1 - (1 - winsor_level)/2), names = F, na.rm = T)
    if (only_top) {
        bounds[1] <- 0
    }
    input_vector[input_vector < bounds[1]] <- bounds[1]
    input_vector[input_vector > bounds[2]] <- bounds[2]
    if ((bounds[2] - bounds[1]) > 0) {
        input_vector <- (input_vector - bounds[1])/(bounds[2] - 
            bounds[1])
    }
    else {
        input_vector[1:length(input_vector)] <- NA
    }
    return(input_vector)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 } winsorising
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
