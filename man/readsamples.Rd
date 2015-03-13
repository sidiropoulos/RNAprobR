\name{readsamples}
\alias{readsamples}
\title{Import of tables prepared by Galaxy workflow to R environment}
\usage{
readsamples(samples, euc = "counts", m = "", k2n_files = "")
}
\arguments{
\item{samples}{vector with paths to unique_barcodes files to be combined}

\item{euc}{method of calculating estimated unique counts (default: "counts")}

\item{m}{random barcode complexity (required if and only if euc="Fu")}

\item{k2n_files}{vector with paths to k2n files corresponding to files given
in samples (required if and only if euc="HRF-Seq"; order important!).
Recycled if necessary}
}
\value{
euc_GR GRanges containing information: 1) seqnames (sequence name;
RNAid) 2) Start, 3) End, 4) EUC value of a given fragment
}
\description{
Function readsamples() reads the output of read processing and mapping
workflow which has to consist of 4 columns 1) RNAid, 2)Insert start,
3)Insert end, 4)Unique barcode count. It combines separate files coming from
the same treatment (e.g. controls) and calculates estimated
unique counts (EUCs) by either (a) keeping unique counts (euc="counts"), (b)
using formula from Fu GK et al. PNAS 2011 (binomial distribution
calculation) (euc="Fu") or (c) using method described in Kielpinski and
Vinther, NAR 2014 (euc="HRF-Seq") If euc="Fu" then the count of all possible
barcodes is required (m), e.g. if you use 7 nucleotide, fully degenerate
random barcodes (NNNNNNN) then m=16384 (m=4**7) If euc="HRF-Seq" then the
path to a precomputed k2n file is required (generate using k2n_calc()
function)(default: "counts")
}
\examples{
write("DummyRNA\\t1\\t2\\t3",file="dummy_unique_barcode")
readsamples(samples = "dummy_unique_barcode", euc = "counts")
}
\author{
Lukasz Jan Kielpinski, Nikos Sidiropoulos
}
\references{
Fu, G.K., Hu, J., Wang, P.H., and Fodor, S.P. (2011). Counting
individual DNA molecules by the stochastic attachment of diverse labels.
Proc Natl Acad Sci U S A 108, 9026-9031.
Kielpinski, L.J., and Vinther, J. (2014).
Massive parallel-sequencing-based hydroxyl radical probing of RNA
accessibility. Nucleic Acids Res.
}
\seealso{
\code{\link{comp}}, \code{\link{plotReads}}, \code{\link{k2n_calc}}
}

