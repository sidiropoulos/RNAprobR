\name{correct_oversaturation}
\alias{correct_oversaturation}
\title{Correcting EUC of oversaturated fragments.}
\usage{
correct_oversaturation(euc_GR, read_counts_file)
}
\arguments{
\item{euc_GR}{GRanges produced by readsamples() function}

\item{read_counts_file}{path to a file with observed read counts.}
}
\value{
euc_GR GRanges analogous to the readsamples() function output,
but with finite EUCs where infinity was present.
}
\description{
If for a given fragment the number of observed unique barcodes is equal to
the total barcode complexity (all combinations of barcodes are associated
with a given fragment), then the readsamples function assignes infinite EUC.
This can be corrected by the function correct_oversaturation(). By comparing
observed read counts with EUCs for other fragments it calculates the
correction factor.
Then, for the oversaturated fragments it multiplies the observed read counts
by the correction factor to estimate EUC. The assumption behind this
correction is that fragments have similar rate of PCR duplicates production.
}
\examples{
write(c("DummyRNA\\t1\\t2\\t1000", "DummyRNA\\t3\\t4\\t1024"),
     file="dummy_unique_barcode")
write(c("DummyRNA\\t1\\t2\\t5000", "DummyRNA\\t3\\t4\\t10000"),
      file="dummy_read_counts")
my_EUCs <- readsamples(samples = "dummy_unique_barcode", euc = "Fu", m=1024)
correct_oversaturation(euc_GR = my_EUCs,
                       read_counts_file = "dummy_read_counts")

}
