#Testing GRanges based workflow:
#Example on published HRF-Seq data:

library("RNAstr")

#Read in functions:
#sapply(list.files(pattern="[.]R$", path="./R", full.names=TRUE), source)

#Calculate k2n vectors for each samples:
# k2n_calc(merged_file="merged_temp14.gz", unique_barcode_file="unique_barcodes14.gz", output_file="k2n_14")
# k2n_calc(merged_file="merged_temp16.gz", unique_barcode_file="unique_barcodes16.gz", output_file="k2n_16")
# k2n_calc(merged_file="merged_temp22.gz", unique_barcode_file="unique_barcodes22.gz", output_file="k2n_22")
# k2n_calc(merged_file="merged_temp24.gz", unique_barcode_file="unique_barcodes24.gz", output_file="k2n_24")

#Read in the samples (as GRanges), converting unique barcodes count to EUC with HRF-Seq method:
treated <- c(paste(path.package("RNAstr"),"/inst/extdata/unique_barcodes14.gz", sep=""), 
             paste(path.package("RNAstr"),"/inst/extdata/unique_barcodes22.gz", sep=""))
control <- c(paste(path.package("RNAstr"),"/inst/extdata/unique_barcodes16.gz", sep=""), 
             paste(path.package("RNAstr"),"/inst/extdata/unique_barcodes24.gz", sep=""))

k2n_treated <- c(paste(path.package("RNAstr"),"/inst/extdata/k2n_14",sep=""),paste(path.package("RNAstr"),"/inst/extdata/k2n_22",sep=""))
k2n_control <- c(paste(path.package("RNAstr"),"/inst/extdata/k2n_16",sep=""),paste(path.package("RNAstr"),"/inst/extdata/k2n_24",sep=""))

control_euc <- readsamples(control, euc="HRF-Seq", k2n_files=k2n_control)
treated_euc <- readsamples(treated, euc="HRF-Seq", k2n_files=k2n_treated)

#Summarize terminations, priming and coverage at each position:
comp_control <- comp(control_euc, cutoff=101, fasta_file=paste(path.package("RNAstr"),"/inst/extdata/hrfseq.fa",sep="")) 
comp_treated <- comp(treated_euc, cutoff=101, fasta_file=paste(path.package("RNAstr"),"/inst/extdata/hrfseq.fa",sep=""))

#Normalize using deltaTCR method, smoothing with window_size=3
normalized <- dtcr(comp_control, comp_treated, window_size=3, nt_offset=1)
normalized <- swinsor(comp_control, add_to=normalized)
normalized <- slograt(comp_control, comp_treated, add_to=normalized)
normalized <- compdata(comp_treated, add_to=normalized)
normalized <- compdata(comp_control, add_to=normalized)

###############################
#Example on LNA-Stop-Seq (much bigger dataset)

R

#Read in functions:
sapply(list.files(pattern="[.]R$", path="/data/lukasz/projects/RNA_R_package/", full.names=TRUE), source)

control=c("/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/ind11/unique_barcodes.gz","/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/ind12/unique_barcodes.gz")
treated=c("/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/ind9/unique_barcodes.gz","/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/ind10/unique_barcodes.gz")

k2n_control=c("/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/ind11/Uf_to_Mf", "/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/ind12/Uf_to_Mf")
k2n_treated=c("/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/ind9/Uf_to_Mf", "/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/ind10/Uf_to_Mf")

control_euc <- readsamples(control, euc="HRF-Seq", k2n_files=k2n_control)
treated_euc <- readsamples(treated, euc="HRF-Seq", k2n_files=k2n_treated)

#Summarize terminations, priming and coverage at each position:
comp_control <- comp(control_euc, fasta_file="/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/choice_of_txs/chosen_for_analysis.fa") 
comp_treated <- comp(treated_euc, fasta_file="/seqdata/krogh/lukasz/130815/trimmed_adapters/Project_oligo_selection/transcript_based/choice_of_txs/chosen_for_analysis.fa")

normalized <- slograt(comp_control, comp_treated)

plotRNA(normalized, RNAid="ENSMUST00000037811", "slograt", "slograt.p")


