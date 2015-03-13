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
write(strwrap("chr1\\t134212702\\t134229870\\tENSMUST00000072177\\t0\\t+
             \\t134212806\\t134228958\\t0\\t8\\t347,121,24,152,66,120,133,1973,
             \\t0,8827,10080,11571,12005,13832,14433,15195,", width = 300),
      file="dummy.bed")
BED2txDb("dummy.bed")
}
\author{
Lukasz Jan Kielpinski, Nikos Sidiropoulos
}

