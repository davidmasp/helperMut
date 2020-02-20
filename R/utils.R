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
  seqlevels(extended_region) = seqlevels(genome)
  gL = seqlengths(genome)
  seqlengths(extended_region) <- gL
  extended_region = trim(extended_region)

  # we now compute the frequency
  seqN = BSgenome::getSeq(genome, extended_region)
  K = (k*2)+1
  tri_freqs = oligonucleotideFrequency(seqN, width=K,...)
  tri_counts = apply(tri_freqs,2, sum)
  return(tri_counts)
}

#' @describeIn get_k_freq Get k-mer frequency in regions
get_k_freq_fromRegion <- function(k,gr,genome = genome_selector(),...){

  # this is not really needed but just in case
  extended_region = extend(x = gr,upstream = k, downstream = k)
  extended_region = GenomicRanges::reduce(extended_region)
  seqlevels(extended_region) = seqlevels(genome)
  gL = seqlengths(genome)
  seqlengths(extended_region) <- gL
  extended_region = trim(extended_region)

  # we now compute the frequency
  seqN = BSgenome::getSeq(genome, extended_region)
  K = (k*2)+1
  tri_freqs = oligonucleotideFrequency(seqN, width=K,...)
  tri_counts = apply(tri_freqs,2, sum)
  return(tri_counts)
}


compute_rapo <- function(vr,genome = genome_selector(),wl = 1000000){
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
