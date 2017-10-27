\name{BED2txDb}
\alias{BED2txDb}
\title{Bedgraph to TranscriptDb object}
\usage{
BED2txDb(input_bed_path)
}
\arguments{
\item{input_bed_path}{Path to BED file. If 12 column BED provided, function
is splice aware. If 6 column BED provided, function assumes no splicing.}
}
\value{
TranscriptDb object
}
\description{
Function to transform BED format file to Bioconductor TranscriptDb object
}
\examples{

txDb <- BED2txDb(system.file("extdata", "dummy.bed", package = "RNAprobR",
                             mustWork = TRUE))

}
\author{
Lukasz Jan Kielpinski, Nikos Sidiropoulos
}
