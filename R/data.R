# $$$$$$$\   $$$$$$\ $$$$$$$$\  $$$$$$\
# $$  __$$\ $$  __$$\\__$$  __|$$  __$$\
# $$ |  $$ |$$ /  $$ |  $$ |   $$ /  $$ |
# $$ |  $$ |$$$$$$$$ |  $$ |   $$$$$$$$ |
# $$ |  $$ |$$  __$$ |  $$ |   $$  __$$ |
# $$ |  $$ |$$ |  $$ |  $$ |   $$ |  $$ |
# $$$$$$$  |$$ |  $$ |  $$ |   $$ |  $$ |
# \_______/ \__|  \__|  \__|   \__|  \__|


#' Ordered mutation types based on Supek,Lehner Cell 2017
#'
#' A vector containing the ordered mutation types at k=1
#' coming from Supek,Lehner Cell 2017
#'
#' @format A vector with 96 mutation types
#' @source \url{https://www.sciencedirect.com/science/article/pii/S0092867417307742}
"order_ms96_supekCell2017"


#' Ordered mutation types based on Cosmic signatures website
#'
#' A vector containing the ordered mutation types at k=1
#' from the sanger cosmic website
#'
#' @format A vector with 96 mutation types
#' @source \url{https://cancer.sanger.ac.uk/signatures_v2/Signature-2.png}
"order_ms96_cosmicSignatures"


#' Colors used to plot mutational signatures
#'
#' A vector containing the HEX codes for each mutation types at k=0.
#'
#' It can directly be used in scale_manual from ggplot2
#'
#' @format A names vector with 6 hex codes.
#' @source \url{https://cancer.sanger.ac.uk/signatures_v2/Signature-2.png}
"tr_colors"



#' IUPAC DNA codes
#'
#' A list of each DNA letter and its representative set in ACTG letters.
#'
#' @format A list with all IUPAC DNA codes and its meanings
#' @source \url{https://www.bioinformatics.org/sms/iupac.html}
"dna_codes"
