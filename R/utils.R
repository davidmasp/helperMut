# $$\   $$\ $$$$$$$$\ $$$$$$\ $$\       $$$$$$\
# $$ |  $$ |\__$$  __|\_$$  _|$$ |     $$  __$$\
# $$ |  $$ |   $$ |     $$ |  $$ |     $$ /  \__|
# $$ |  $$ |   $$ |     $$ |  $$ |     \$$$$$$\
# $$ |  $$ |   $$ |     $$ |  $$ |      \____$$\
# $$ |  $$ |   $$ |     $$ |  $$ |     $$\   $$ |
# \$$$$$$  |   $$ |   $$$$$$\ $$$$$$$$\\$$$$$$  |
#  \______/    \__|   \______|\________|\______/




.lenunique <- function(x){
  length(unique(x))
}

#' Jaccard Index
#'
#' @description
#'
#' Computes the jaccard index between two string vectors.
#'
#' @param x a character ector
#' @param y a second character vector
#'
#' @return the Jaccard index
#' @export
#'
#' @examples
#'
#' set1 = c("a","b","c")
#' set2 = c("b","d","t")
#' jaccard(set1,set2)
#'
jaccard <- function(x,y){
  length(intersect(x,y))/length(union(x,y))
}



# STATS ===============

#' pseudo-vectorized binomial test
#'
#' This function is a wrapper of the \link[stats]{binom.test}
#' function that can be used with vectors.
#'
#' @param x x for binom test
#' @param n n for binom test
#' @param p success probability
#' @param ... further params for binom.test
#'
#' @return a dataframe with statistics from the binom.test.
#'
#' @export
#'
#' @examples
#'
#' binom_test(x = c(1,2,3),n = c(2,4,6))
#'
binom_test <- function(x,n,p=NULL,...){

  if (is.null(p)){
    res = purrr::map_df(1:length(x),function(i){
      broom::tidy(stats::binom.test(x = x[i],n = n[i],...))
    })
  } else if (!is.null(n) ){
    res = purrr::map_df(1:length(x),function(i){
      tryCatch({
        test = stats::binom.test(x = x[i],n = n[i],p = p[i],...)
        if (is.logical(test$p.value)){
          test$p.value = as.numeric(test$p.value)
        }
        broom::tidy(test)
      }, error = function(e) {browser()})

    })
  }
  return(res)
}

#' Get k-mer frequency in mutation-neighbouring regions
#'
#' From mutations in a VRanges object, obtain the k-mer oligonucleotide
#' frequencies at a given window length.
#'
#' @param vr a VRanges object with SNVs
#' @param k the extension of mutation context in one direction
#' @param wl the window length from which each mutation would be expanded
#' @param gr a GRanges object delimiting the region of interest
#' @param genome a BSgenome object
#' @param ... arguments for oligonucleotideFrequency
#'
#' @return
#' A integer named vector with the oligonucleotides selected at a given k.
#' @export
#'
#' @examples
#' get_k_freq(vr_example,
#'           k = 1,
#'           wl = 1000,
#'           genome = genome_selector("Hsapiens.1000genomes.hs37d5"))
get_k_freq <- function(vr,k,wl,genome = genome_selector(),...){

  up = round(wl/2,0)
  down = round(wl/2,0)
  extended_region = extend(x = vr,upstream = up, downstream = down)

  # we reduce to avoid overlaps
  extended_region = GenomicRanges::reduce(extended_region)
  GenomeInfoDb::seqlevels(extended_region) = GenomeInfoDb::seqlevels(genome)
  gL = GenomeInfoDb::seqlengths(genome)
  GenomeInfoDb::seqlengths(extended_region) <- gL
  extended_region = GenomicRanges::trim(extended_region)

  # we now compute the frequency
  seqN = BSgenome::getSeq(genome, extended_region)
  K = (k*2)+1
  tri_freqs = Biostrings::oligonucleotideFrequency(seqN, width=K,...)
  tri_counts = apply(tri_freqs,2, sum)
  return(tri_counts)
}

#' @describeIn get_k_freq Get k-mer frequency in regions
get_k_freq_fromRegion <- function(k,gr,genome = genome_selector(),...){

  # this is not really needed but just in case
  extended_region = extend(x = gr,upstream = k, downstream = k)
  extended_region = GenomicRanges::reduce(extended_region)
  GenomeInfoDb::seqlevels(extended_region) = GenomeInfoDb::seqlevels(genome)
  gL = GenomeInfoDb::seqlengths(genome)
  GenomeInfoDb::seqlengths(extended_region) <- gL
  extended_region = GenomicRanges::trim(extended_region)

  # we now compute the frequency
  seqN = BSgenome::getSeq(genome, extended_region)
  K = (k*2)+1
  tri_freqs = Biostrings::oligonucleotideFrequency(seqN, width=K,...)
  tri_counts = apply(tri_freqs,2, sum)
  return(tri_counts)
}


#' Compute mutation motif enrichemnt
#'
#' Given a VRanges object which can contain either a sample, a set of samples
#' or a fraction of a sample, calculate the enrichment of a motif compared
#' to a control motif.
#'
#' A given window size (k_offset*2+1) is used to compute ntps at risk.
#'
#' Implemets the the enrichment measure used in
#' Roberts Nature Genetics (2013) for apobec enrichment in mutation
#' clusters.
#'
#' In brief, the enrichment is calculated as:
#'
#' \deqn{E = ( muts_set * ctx_control ) / (muts_control * ctx_set)}
#'
#' Note that this value is similar to the odds ratio, however, in the
#' fisher test and its odds ratio, the control set excludes the test set
#' while it's included in the enrichment measure.
#'
#' @param vr a VRanges object with mutation data
#' @param genome a BSgenome object
#' @param k_offset k for window offset, used to control for ntp at risk
#' @param set a mutation set to test
#' @param control_set a mutation set to control for
#'
#' @return
#'
#' A dataframe with tidied results of the fisher test. The dataframe contains
#' the enrichment value as a column.
#'
#' @export
#'
#' @examples
compute_motif_enrichment <- function(vr,
                                     genome,
                                     k_offset = 20,
                                     set= "TCW>K",
                                     control_set = "NCN>D") {

  if (!is(object = genome,"BSgenome")){
    stop(message = "Genome is not a BSgenome object")
  }

  if (!is(object = vr,"VRanges")){
    stop(message = "VR is not a VRanges object")
  }

  set_test = make_set(set)
  set_control = make_set(control_set)
  set_control_excluded = set_control[!set_control %in% set_test]
  set_params = identify_mut_aestetics(set_test)

  ctx_test = substr(set_test,1,set_params$K)
  ctx_control = substr(set_control,1,set_params$K)
  ctx_control_excluded = ctx_control[!ctx_control %in% ctx_test]

  # I extend the muts
  gr = as(object = vr,Class = "GRanges")
  # this is equivalent to a total wl of k_offset*2+1 (default is 41)
  gr_extended = extend(x = gr,upstream = k_offset,downstream = k_offset)
  seqs = Biostrings::getSeq(genome,gr_extended)
  oligos = Biostrings::oligonucleotideFrequency(seqs,width = set_params$K)

  ms = get_MS_VR(x = vr,k=set_params$k,
                 genome = genome,
                 simplify_set = set_params$mset_ref)
  ms = ms[!is.na(ms)]
  ms = ms[ms %in% c(set_test,set_control)]

  muts_test = sum(ms %in% set_test)
  muts_control = sum(ms %in% set_control)
  muts_control_excluded = sum(ms %in% set_control_excluded)

  oligos_sum = apply(oligos, 2, sum)
  nnames = simplify_ctx(names(oligos_sum),
                        simplify_set = set_params$mset_ref)
  oligos_sum = oligos_sum[!is.na(nnames)]
  names(oligos_sum) = nnames[!is.na(nnames)]
  oligos_sum %>% split(names(oligos_sum)) %>% purrr::map_dbl(sum) -> oligos_sum

  oligos_test = sum(oligos_sum[ctx_test])
  oligos_control = sum(oligos_sum[ctx_control])
  oligos_control_ex = sum(oligos_sum[ctx_control_excluded])

  enrichment = (muts_test * oligos_control) / (muts_control * oligos_test)

  matrix(c(muts_test,
           oligos_test,
           muts_control_excluded,
           oligos_control_ex),nrow = 2,ncol = 2) -> cnt_table

  fisher.test(x = cnt_table) %>% broom::tidy() -> res_df
  res_df$enrichment = enrichment

  return(res_df)

}

#' Compute apobec enrichment locally for a set of mutations
#'
#' Implemented from the methods in Seplyarskiy 2016 Genome Research.
#'
#' It computes the number of APOBEC mutations defined a TCW>K and it compares
#' these numbers to nonAPOBEC mutations in the VCW>K context surrounding the
#' given set of mutations. The mutations are extended with a wl described in
#' the arguments. The ratios are corrected by context occurence.
#'
#' @param vr a VRanges object with a set of mutations
#' @param genome a BSgenome object
#' @param wl the window lenght to explore
#'
#' @return
#' A single integer indicating the ratio of rates. Higher values indicate
#' more apobec3a/b enrichment.
#' @export
#'
#' @examples
compute_rAPOBEC <- function(vr,genome = genome_selector(),wl = 1000000){
  # defined as in Seplyarskiy_2016_GR and Roberts 2013 Nature Genetics
  # but locally (1mb bins)

  ms = get_MS_VR(x = vr,k=1)

  set_up = make_set("TCW>K")
  set_down = make_set("VCW>K")

  apobec_muts = sum(ms %in% set_up)
  non_apobec_muts = sum(ms %in% set_down)

  tri_counts = get_k_freq(vr = vr,k = 1,wl = wl,genome = genome)

  set_up_ctx = set_up %>% stringr::str_sub(1,3) %>% unique()
  set_down_ctx = set_down %>% stringr::str_sub(1,3) %>% unique()

  ctx_up_freq = sum(tri_counts[set_up_ctx])
  ctx_down_freq = sum(tri_counts[set_down_ctx])

  # we compute the total number of positions
  #total_lenght = sum(width(extended_region))
  #this is not needed.

  ratio_up = apobec_muts / non_apobec_muts

  ratio_down = ctx_up_freq / ctx_down_freq

  rr = ratio_up / ratio_down

  return(rr)

}
