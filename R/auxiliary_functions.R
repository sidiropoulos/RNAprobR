###Auxiliary functions (need to be read-in with a package):
#Stacking offset vectors, used for moving average and moving sum calculations. It forms the matrix with number of rows equivalent to length of sliding window and in each row there is input vector offset by one position:






#' Auxiliary functions for package RNAstr
#' 
#' Auxiliary functions for package RNAstr
#' 
#' 




#' Auxiliary functions for package RNAstr
#' 
#' Auxiliary functions for package RNAstr
#' 
#' Make smoothing matrix
#' 
#' Auxiliary function speeding up calculations in functions operating on
#' sliding windows.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param input_vector input_vector
#' @param window_size window_size
#' @return Matrix with (length(input_vector)+window_size-1) columns and
#' window_size rows. Each row contains input vector values starting at column
#' equal to row number.
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' construct_smoothing_matrix(1:10, 3)
#' 
#' ## The function is currently defined as
#' function (input_vector, window_size)
#' {
#'     vector_length <- length(input_vector)
#'     my_mat <- matrix(nrow = window_size, ncol = (vector_length +
#'         window_size - 1))
#'     for (i in 1:(window_size)) {
#'         my_mat[i, i:(vector_length + i - 1)] <- input_vector
#'     }
#'     return(my_mat)
#'   }
#' 
#' @export construct_smoothing_matrix
construct_smoothing_matrix <- function(input_vector, window_size){
	vector_length <- length(input_vector)
	my_mat <- matrix(nrow=window_size, ncol=(vector_length+window_size-1)) #Create empty matrix.
	for(i in 1:(window_size)){my_mat[i,i:(vector_length+i-1)] <- input_vector}   #Fill the matrix, offsetting by 1 in each cycle.
	return(my_mat)
}

#calculating moving average of a vector:










#' Calculate moving average
#' 
#' Moving average calculation. NA values ignored.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param input_vector input_vector
#' @param window_size window_size
#' @return numeric vector. Each position represents the mean of values in
#' window_size positions centred at this position.
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' moving_average(1:10, 5)
#' 
#' ## The function is currently defined as
#' function (input_vector, window_size)
#' {
#'     window_side <- window_size/2 - 0.5
#'     return(colMeans(construct_smoothing_matrix(input_vector,
#'         window_size), na.rm = T)[(window_side + 1):(length(input_vector) +
#'         window_side)])
#'   }
#' 
#' @export moving_average
moving_average <- function(input_vector, window_size){
	window_side <- window_size/2-0.5
	return(colMeans(construct_smoothing_matrix(input_vector, window_size), na.rm=T)[(window_side+1):(length(input_vector)+window_side)])
}

#Winsorize a vector:










#' Winsor normalization with fitting to <0,1> range.
#' 
#' Function performs Winsor normalization of a supplied vector, sets bottom
#' (winsorized) value to 0, top 1 and linearly transforms all the values to the
#' <0,1> range.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param input_vector Vector with values to be Winsorized
#' @param winsor_level Winsorization level. Bottom outliers will be set to
#' (1-winsor_level)/2 quantile and top outliers to (1+winsor_level)/2 quantile.
#' @param only_top If TRUE then bottom values are not Winsorized and the lowest
#' is set to 0.
#' @return Vector of numerics within <0,1>.
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references Hastings, Cecil; Mosteller, Frederick; Tukey, John W.; Winsor,
#' Charles P. Low Moments for Small Samples: A Comparative Study of Order
#' Statistics. The Annals of Mathematical Statistics 18 (1947), no. 3,
#' 413--426.
#' @keywords winsorising ~kwd1 ~kwd2
#' @examples
#' 
#' data_set <- runif(1:100)*100
#' plot(winsor(data_set, winsor_level=0.8) ~ data_set)
#' 
#' 
#' ## The function is currently defined as
#' function (input_vector, winsor_level = 0.9, only_top = F)
#' {
#'     bounds <- quantile(input_vector, c((1 - winsor_level)/2,
#'         1 - (1 - winsor_level)/2), names = F, na.rm = T)
#'     if (only_top) {
#'         bounds[1] <- 0
#'     }
#'     input_vector[input_vector < bounds[1]] <- bounds[1]
#'     input_vector[input_vector > bounds[2]] <- bounds[2]
#'     if ((bounds[2] - bounds[1]) > 0) {
#'         input_vector <- (input_vector - bounds[1])/(bounds[2] -
#'             bounds[1])
#'     }
#'     else {
#'         input_vector[1:length(input_vector)] <- NA
#'     }
#'     return(input_vector)
#'   }
#' 
#' @export winsor
winsor <- function(input_vector, winsor_level=0.9, only_top=F){ #!!! Function changes integers to numerics!
  bounds <- quantile(input_vector, c((1-winsor_level)/2, 1-(1-winsor_level)/2), names=F, na.rm=T)
  if(only_top){bounds[1] <- 0}
  input_vector[input_vector < bounds[1]] <- bounds[1]
  input_vector[input_vector > bounds[2]] <- bounds[2]
  if((bounds[2]-bounds[1])>0){input_vector <- (input_vector-bounds[1])/(bounds[2]-bounds[1])}else{input_vector[1:length(input_vector)] <- NA}
  return(input_vector)
}

#Sliding-winsorization of a vector:










#' Smooth Winsor Normalization
#' 
#' Function performs Winsor normalization (see winsor() function) of each
#' window of specified window_size, sliding in a given vector by 1 position,
#' and reports a list of (1) mean Winsorized values for each vector position
#' (mean of Winsorized value for a given position as calculated within each
#' overlapping window) and (2) standard deviation of those Winsorized values.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param input_vector Vector with values to be smooth-Winsorized
#' @param window_size Size of a sliding window.
#' @param winsor_level Winsorization level. Bottom outliers will be set to
#' (1-winsor_level)/2 quantile and top outliers to (1+winsor_level)/2 quantile.
#' @param only_top If TRUE then bottom values are not Winsorized and are set to
#' 0.
#' @return \item{comp1}{Vector with mean Winsorized values for each
#' input_vector position} \item{comp2 }{Vector with standard deviation of
#' Winsorized values for each input_vector position}
#' 
#' %% ...
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references SHAPES publication. Poulsen, Kielpinski, Salama, Krogh, Vinther.
#' In review as for 2nd Oct 2014.
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' data_set <- runif(1:100)*100
#' plot(swinsor_vector(data_set, window_size=71, winsor_level=0.8)[[1]] ~ data_set)
#' 
#' 
#' ## The function is currently defined as
#' function (input_vector, window_size, winsor_level = 0.9, only_top = F)
#' {
#'     align_smoothing_matrix <- function(input_matrix) {
#'         vector_length <- ncol(input_matrix) - nrow(input_matrix) +
#'             1
#'         my_mat <- matrix(nrow = nrow(input_matrix), ncol = vector_length)
#'         for (i in 1:nrow(input_matrix)) {
#'             my_mat[i, ] <- input_matrix[i, i:(vector_length +
#'                 i - 1)]
#'         }
#'         return(my_mat)
#'     }
#'     my_matrix <- construct_smoothing_matrix(input_vector, window_size)
#'     winsorized_matrix <- apply(my_matrix, 2, winsor, winsor_level = winsor_level,
#'         only_top = only_top)
#'     winsorized_matrix[, c(1:(window_size - 1), (length(input_vector) +
#'         1):ncol(winsorized_matrix))] <- NA
#'     aligned_winsorized_matrix <- align_smoothing_matrix(winsorized_matrix)
#'     means_vector <- colMeans(aligned_winsorized_matrix, na.rm = T)
#'     sds_vector <- apply(aligned_winsorized_matrix, 2, FUN = sd,
#'         na.rm = T)
#'     return(list(means_vector, sds_vector))
#'   }
#' 
#' @export swinsor_vector
swinsor_vector <- function(input_vector, window_size, winsor_level=0.9, only_top=F){
	#Define function for aligning matrix generated by construct_smoothing_matrix
	align_smoothing_matrix <- function(input_matrix){
		vector_length <- ncol(input_matrix)-nrow(input_matrix)+1 #Length of the vector used to construct the matrix in construct_smoothing_matrix function.
		my_mat <- matrix(nrow=nrow(input_matrix), ncol=vector_length) #Create empty matrix.
		for(i in 1:nrow(input_matrix)){my_mat[i,] <- input_matrix[i,i:(vector_length+i-1)]}   #Fill the matrix, offseting by -1 in each cycle.
		return(my_mat)
	}
	#End of defining function.
	my_matrix <- construct_smoothing_matrix(input_vector, window_size)
	winsorized_matrix <- apply(my_matrix, 2, winsor, winsor_level=winsor_level, only_top=only_top)
	winsorized_matrix[,c(1:(window_size-1), (length(input_vector)+1):ncol(winsorized_matrix))] <- NA #To control the border behaviour
	aligned_winsorized_matrix <- align_smoothing_matrix(winsorized_matrix)
	means_vector <- colMeans(aligned_winsorized_matrix, na.rm=T)
	sds_vector <- apply(aligned_winsorized_matrix,2, FUN=sd, na.rm=T)
	return(list(means_vector, sds_vector))
}

#Convert GRanges to norm data frame:










#' Function to make data frame out of GRanges output of normalizing functions
#' (dtcr(), slograt(), swinsor(), compdata()) for all or a set of chosen
#' transcripts in the file.
#' 
#' Simple convenience function.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param norm_GR GRanges object made by other normalization function (dtcr(),
#' slograt(), swinsor(), compdata()) from which data is to be extracted
#' @param RNAid Transcript identifiers of transcripts that are to be extracted
#' @param norm_methods Names of the columns to be extracted.
#' @return Data frame object with columns: RNAid, Pos and desired metadata
#' columns (e.g. nt, dtcr)
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski, Nikos Sidiropoulos
#' @seealso \code{\link{norm_df2GR}}, \code{\link{dtcr}},
#' \code{\link{swinsor}}, \code{\link{slograt}}, \code{\link{compdata}}
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' dummy_euc_GR_treated <- GRanges(seqnames="DummyRNA", IRanges(start=round(runif(100)*100), width=round(runif(100)*100+1)), strand="+", EUC=round(runif(100)*100))
#' dummy_comp_GR_treated <- comp(dummy_euc_GR_treated)
#' dummy_swinsor <- swinsor(dummy_comp_GR_treated)
#' GR2norm_df(dummy_swinsor)
#' 
#' 
#' ## The function is currently defined as
#' function (norm_GR, RNAid="all", norm_methods="all")
#' {
#'     if(norm_methods=="all" || missing(norm_methods))
#'     {
#'         norm_methods <- names(mcols(norm_GR))
#'     }
#' 
#'     if(RNAid=="all" || missing(RNAid))
#'     {
#'         RNAid <- levels(seqnames(norm_GR))
#'     }
#' 
#'     norm_GR <- norm_GR[seqnames(norm_GR) %in% RNAid, norm_methods]
#'     data.frame(RNAid = as.character(seqnames(norm_GR)), Pos = as.integer(start(norm_GR)),
#'         mcols(norm_GR))
#'   }
#' 
#' @export GR2norm_df
GR2norm_df <- function(norm_GR, RNAid="all", norm_methods="all"){
  
  if(norm_methods=="all" || missing(norm_methods)){norm_methods <- names(mcols(norm_GR))}
  if(RNAid=="all" || missing(RNAid)){RNAid <- levels(seqnames(norm_GR))}
  
  # Parsing normalized file based on the selected RNAid(s) and normalization method(s).
  norm_GR <- norm_GR[seqnames(norm_GR) %in% RNAid, norm_methods]
	data.frame(RNAid=as.character(seqnames(norm_GR)), Pos=as.integer(start(norm_GR)), mcols(norm_GR))
}


#Convert norm data frame to GRanges:










#' Function to convert norm_df data frame (product of GR2norm_df()) to GRanges
#' 
#' Function to convert norm_df data frame (product of GR2norm_df()) to GRanges.
#' 
#' %% ~~ If necessary, more details than the description above ~~
#' 
#' @param norm_df norm_df data frame needs to have columns: RNAid (equivalent
#' to seqnames in GRanges) and Pos (equivalent to start in GRanges) and
#' metadata
#' @return GRanges compatible with objects created by normalizing functions
#' (dtcr(), slograt(), swinsor(), compdata())
#' @note %% ~~further notes~~
#' @author Lukasz Jan Kielpinski
#' @seealso %% ~~objects to See Also as \code{\link{help}}, ~~~
#' @references %% ~put references to the literature/web site here ~
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' dummy_norm_df <- data.frame(RNAid="dummyRNA", Pos=1:100, my_data1=runif(1:100))
#' norm_df2GR(dummy_norm_df)
#' 
#' ## The function is currently defined as
#' function (norm_df)
#' {
#'     GRanges(seqnames = norm_df$RNAid, IRanges(start = norm_df$Pos,
#'         width = 1), strand = "+", norm_df[names(norm_df) !=
#'         "RNAid" & names(norm_df) != "Pos"])
#'   }
#' 
#' @export norm_df2GR
norm_df2GR <- function(norm_df){
	GRanges(seqnames=norm_df$RNAid, IRanges(start=norm_df$Pos, width=1), strand="+", norm_df[names(norm_df)!="RNAid" & names(norm_df)!="Pos"])
}
