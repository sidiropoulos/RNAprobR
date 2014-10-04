\name{swinsor_vector}
\alias{swinsor_vector}
\title{
Smooth Winsor Normalization
}
\description{
Function performs Winsor normalization (see winsor() function) of each window of specified window_size, sliding in a given vector by 1 position, and reports a list of (1) mean Winsorized values for each vector position (mean of Winsorized value for a given position as calculated within each overlapping window) and (2) standard deviation of those Winsorized values.
}
\usage{
swinsor_vector(input_vector, window_size, winsor_level = 0.9, only_top = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input_vector}{
Vector with values to be smooth-Winsorized
}
  \item{window_size}{
Size of a sliding window.
}
  \item{winsor_level}{
Winsorization level. Bottom outliers will be set to (1-winsor_level)/2 quantile and top outliers to (1+winsor_level)/2 quantile.
}
  \item{only_top}{
If TRUE then bottom values are not Winsorized and are set to 0.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{comp1 }{Vector with mean Winsorized values for each input_vector position}
\item{comp2 }{Vector with standard deviation of Winsorized values for each input_vector position}

%% ...
}
\references{
SHAPES publication. Poulsen, Kielpinski, Salama, Krogh, Vinther. In review as for 2nd Oct 2014.
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
plot(swinsor_vector(data_set, window_size=71, winsor_level=0.8)[[1]] ~ data_set)


## The function is currently defined as
function (input_vector, window_size, winsor_level = 0.9, only_top = F) 
{
    align_smoothing_matrix <- function(input_matrix) {
        vector_length <- ncol(input_matrix) - nrow(input_matrix) + 
            1
        my_mat <- matrix(nrow = nrow(input_matrix), ncol = vector_length)
        for (i in 1:nrow(input_matrix)) {
            my_mat[i, ] <- input_matrix[i, i:(vector_length + 
                i - 1)]
        }
        return(my_mat)
    }
    my_matrix <- construct_smoothing_matrix(input_vector, window_size)
    winsorized_matrix <- apply(my_matrix, 2, winsor, winsor_level = winsor_level, 
        only_top = only_top)
    winsorized_matrix[, c(1:(window_size - 1), (length(input_vector) + 
        1):ncol(winsorized_matrix))] <- NA
    aligned_winsorized_matrix <- align_smoothing_matrix(winsorized_matrix)
    means_vector <- colMeans(aligned_winsorized_matrix, na.rm = T)
    sds_vector <- apply(aligned_winsorized_matrix, 2, FUN = sd, 
        na.rm = T)
    return(list(means_vector, sds_vector))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
