#' ce11_all_REs
#'
#' Regulatory elements annotated in C. elegans (ce11) according to Serizay et 
#' al. 2020, "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce11_all_REs)
#'
#' @format GRanges
#'
#' @keywords datasets
#'
#' @references Serizay et al. 2020, 
#' "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#' (\href{https://doi.org/10.1101/2020.02.20.958579}{DOI})
#'
#' @source \href{https://doi.org/10.1101/2020.02.20.958579}{BiorXiv}
#'
#' @examples
#' data(ce11_all_REs)
#' table(ce11_all_REs$regulatory_class)
#' table(ce11_all_REs$which.tissues)
"ce11_all_REs"

#' ce11_proms
#'
#' Promoters annotated in C. elegans (ce11) according to Serizay et 
#' al. 2020, "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#'
#' @docType data
#'
#' @usage data(ce11_proms)
#'
#' @format An object of class \code{"GRanges"}.
#'
#' @keywords datasets
#'
#' @references Serizay et al. 2020, 
#' "Tissue-specific profiling reveals distinctive regulatory 
#' architectures for ubiquitous, germline and somatic genes", BiorXiv.
#' (\href{https://doi.org/10.1101/2020.02.20.958579}{DOI})
#'
#' @source \href{https://doi.org/10.1101/2020.02.20.958579}{BiorXiv}
#'
#' @examples
#' data(ce11_proms)
#' table(ce11_proms$which.tissues)
"ce11_proms"

#' REB1_sacCer3
#'
#' High-confidence REB1 binding sites, from 
#' Rossi, Lai & Pugh 2018 Genome Research
#'
#' @docType data
#'
#' @usage data(REB1_sacCer3)
#'
#' @format An object of class \code{"GRanges"}.
#'
#' @keywords datasets
#'
#' @references Rossi, Lai & Pugh 2018 Genome Research
#' 
#' @examples
#' data(REB1_sacCer3)
#' REB1_sacCer3
"REB1_sacCer3"

#' CTCF_hg38
#'
#' high-score CTCF binding motifs, obtained from JASPAR 
#'
#' @docType data
#'
#' @usage data(CTCF_hg38)
#'
#' @format An object of class \code{"GRanges"}.
#'
#' @keywords datasets
#' 
#' @examples
#' data(CTCF_hg38)
#' CTCF_hg38
"CTCF_hg38"

#' bam_test
#'
#' A .bam file sample
#'
#' @docType data
#'
#' @usage data(bam_test)
#'
#' @format An object of class \code{"GRanges"}.
#'
#' @keywords datasets
#' 
#' @examples
#' data(bam_test)
#' bam_test
"bam_test"

