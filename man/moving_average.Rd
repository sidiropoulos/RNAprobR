\name{moving_average}
\alias{moving_average}
\title{
Calculate moving average
}
\description{
Moving average calculation. NA values ignored.
}
\usage{
moving_average(input_vector, window_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_vector}{
input_vector
}
  \item{window_size}{
window_size
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
numeric vector. Each position represents the mean of values in window_size positions centred at this position.
}
\references{
%% ~put references to the literature/web site here ~
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

moving_average(1:10, 5)

## The function is currently defined as
function (input_vector, window_size) 
{
    window_side <- window_size/2 - 0.5
    return(colMeans(construct_smoothing_matrix(input_vector, 
        window_size), na.rm = T)[(window_side + 1):(length(input_vector) + 
        window_side)])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
