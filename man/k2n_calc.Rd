\name{k2n_calc}
\alias{k2n_calc}
\title{Calculate number of Estimated Unique Counts (EUC's)
corresponding to given number of observed unique barcodes.}
\usage{
k2n_calc(merged_file, unique_barcode_file, output_file)
}
\arguments{
\item{merged_file}{path to merged_temp file containing 4 column: 1) RNAid,
2) Start, 3) End, 4) Barcode sequence (required)}

\item{unique_barcode_file}{character with path to unique_barcode file
(required)}

\item{output_file}{name of a file to be generated (if specified
[recommended] function will write a file, if not - function will return a
vector)}
}
\value{
If output_file specified function writes a file, if not - returns a
vector.
}
\description{
Function calculates EUC's for each number of observed barcodes accounting
for differential ligation probability of different barcodes. Function
k2n_calc() writes file with a vector in which an i-th element is an
estimated unique count given observing i unique barcodes.
}
\examples{

write(c("DummyRNA\\t1\\t1\\tA", "DummyRNA\\t1\\t1\\tC", "DummyRNA\\t2\\t2\\tG",
        "DummyRNA\\t2\\t2\\tT"),file="dummy_merged_file")
write(c("DummyRNA\\t1\\t1\\t2", "DummyRNA\\t2\\t2\\t2"),
        file="dummy_unique_barcode")
k2n_calc(merged_file = "dummy_merged_file",
        unique_barcode_file = "dummy_unique_barcode")

}
\references{
Kielpinski, L.J., and Vinther, J. (2014). Massive
parallel-sequencing-based hydroxyl radical probing of RNA accessibility.
Nucleic Acids Res.
}
\seealso{
\code{\link{readsamples}}
}
\author{
Lukasz Jan Kielpinski, Nikos Sidiropoulos
}
