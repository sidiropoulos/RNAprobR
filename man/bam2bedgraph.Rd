\name{bam2bedgraph}
\alias{bam2bedgraph}
\title{Function converts bam file to bedgraph by counting number of reads starting
at each position (termination counts).
It creates two-track bedgraph file (one track for each strand).}
\usage{
bam2bedgraph(bam_path, allowed_flags = 0:4095, maxMemory = 8000,
  genome_build, bedgraph_out_file = "out_file", track_name = "Track_name",
  track_description = "Track_description")
}
\arguments{
\item{bam_path}{path to a bam file to be converted}

\item{allowed_flags}{integer vector with SAM flags should be kept, see
https://broadinstitute.github.io/picard/explain-flags.html for explanation}

\item{maxMemory}{maxMemory of scanBam function used internally}

\item{genome_build}{character specifying which UCSC genome build should data
be displayed in, e.g. "mm9"}

\item{bedgraph_out_file}{character specifying prefix of output file.
Generated file name is: prefix.bedgraph; if file with such a name already
exists new tracks will be appended.}

\item{track_name}{character specifying track name}

\item{track_description}{character specifying track description}
}
\value{
NULL. Creates a two-track bedgraph file (one track for each strand).
}
\description{
Function converts bam file to bedgraph by counting number of reads starting
at each position (termination counts).
It creates two-track bedgraph file (one track for each strand).
}
\author{
Lukasz Jan Kielpinski
}
