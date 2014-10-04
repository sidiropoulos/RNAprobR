\name{construct_smoothing_matrix}
\alias{construct_smoothing_matrix}
\title{
Make smoothing matrix
}
\description{
Auxiliary function speeding up calculations in functions operating on sliding windows.
}
\usage{
construct_smoothing_matrix(input_vector, window_size)
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
Matrix with (length(input_vector)+window_size-1) columns and window_size rows. Each row contains input vector values starting at column equal to row number.
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
construct_smoothing_matrix(1:10, 3)

## The function is currently defined as
function (input_vector, window_size) 
{
    vector_length <- length(input_vector)
    my_mat <- matrix(nrow = window_size, ncol = (vector_length + 
        window_size - 1))
    for (i in 1:(window_size)) {
        my_mat[i, i:(vector_length + i - 1)] <- input_vector
    }
    return(my_mat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
