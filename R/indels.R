# INDELS
#
# This code should wrap up indel related funs.
#
# Code comes originally from http://fsupeksvr.irbbarcelona.pcb.ub.es/gitlab/dmas/HealthyTissues/blob/master/utils.R


# indel_microhomology_VR(dat[1:100],
#                        genome = genome_selector(alias = "Hsapiens.1000genomes.hs37d5"))


#' Get microhomology of deletions
#'
#' @param vr a VRanges object with deletions
#' @param genome a BSgenome object with tha assembly
#'
#' @return a VRanges object with added meta columns
#' @export
#'
#' @examples
del_microhomology_VR <- function(vr,genome = genome_selector()){
  # TODO, add a insertion check. This is more problematic than u may think.

  # this extends the sequence to the right and the left by the same
  # amount of length. This generates 3 segments with equal length,
  # which is the same as the original one.
  vr_ext = extend(vr,
                                     upstream = width(vr),
                                     downstream = width(vr))

  # we obtain the sequence for the segments.
  vr_seq = Biostrings::getSeq(genome,
                              vr_ext)

  # see below.
  dframe = detect_microhomology(as.character(vr_seq))

  # we put the results back to the vr
  mcols(vr) = cbind(mcols(vr),dframe)

  return(vr)
}


detect_microhomology <- function(seq_vector){
  # check that all regions have same width
  stopifnot(all(nchar(seq_vector)%%3 == 0))
  l = nchar(seq_vector)

  del_l = l /3 # this works because of the previous expansion.

  # this generates a vector per seq stored in a list
  seq_vector_split = stringr::str_split(string = seq_vector,pattern = "")

  purrr::map2_df(.x = seq_vector_split,.y = del_l,.f = function(str,len){
    # We have to check both directions because the homology can be in
    # both strands. There is no complementarity problems because homology
    # is in both strands.
    # first the positive strand
    mh_pos = compute_mh_length(str[1:len],
                               str[(len+1):(2*len)],
                               direction = "+")
    # neagtive strand
    mh_neg = compute_mh_length(str[(2*len + 1):(3*len)],
                               str[(len+1):(2*len)],
                               direction = "-")

    dframe = data.frame(
      mh_pos = mh_pos,
      mh_neg = mh_neg
    )
    return(dframe)
  }) -> res

  return(res)
}

compute_mh_length <- function(vec1,vec2,direction ="+"){

  # it checks the similiraty of the vectors per element and then
  # we transform to RLE to detect how long the matching sequences are.

  if (direction == "+"){
    vec1 = rev(vec1)
    vec2 = rev(vec2)
  }

  vec_rle = rle(vec1 == vec2)
  if (vec_rle$values[1]){
    len_mh = vec_rle$lengths[1]
  } else {
    len_mh = 0
  }

  return(len_mh)
}


#
# # tests
# seqvec = c("CTACTCGTTCGA",
#            "TTTCCTACTATTTCCGGGTAGGGGGGGGTA")
# detect_microhomology(seqvec)

# 0 3
# 2 0

library(VariantAnnotation)
library(helperMut)
library(magrittr)
fn = Sys.getenv("HARTWIG_example")

dat_vcf = VariantAnnotation::readVcf(fn)

dat_vr = as(object = dat_vcf,Class = "VRanges")
dat_indels = dat_vr[dat_vr$set == "indels"]



## this is assuming that the first letter in REF is shared with ALT and
## corresponds to the refseq

vr = dat_indels
mcols(vr) = NULL

params_check_format = T
if(params_check_format){
  params_genome = helperMut::genome_selector(alias = "Hsapiens.1000genomes.hs37d5")
  refSeq = getSeq(params_genome,vr)
  stopifnot(refSeq == ref(vr))
}

ref_l = nchar(ref(vr))
alt_l = nchar(alt(vr))

stopifnot(!any(ref_l == alt_l))
ins_mask = ref_l < alt_l

vr_ins = vr[ins_mask]
vr_del = vr[!ins_mask]

# this is a bit stupid but how my brain understands this
vr_indel = c(vr_ins,vr_del)
indel_type = c(rep("insertion",length(vr_ins)),
               rep("deletion", length(vr_del)))
params_maxRep = 5

### INSERTIONS ###
ins_affected_seq = stringr::str_sub(alt(vr_ins),start = 2)
ins_length = nchar(ins_affected_seq)

### DELETIONS ###
del_affected_seq = stringr::str_sub(ref(vr_del),start = 2)
del_lenght = nchar(del_affected_seq)

## JOIN, same order than the VR!!! ##
indel_affected_seq = c(ins_affected_seq,del_affected_seq)
indel_length = c(ins_length,del_lenght)

## first I compute the first feature the INSERTION SIZE
base_pair = dplyr::case_when(
  indel_length == 1 ~ indel_affected_seq,
  indel_length >= params_maxRep ~ as.character(glue::glue("{params_maxRep}+bp")),
  TRUE ~ as.character(glue::glue("{indel_length}bp"))
)

## this don't work directly in the case_when because it evaluates all
## the elements in the vector and then extracts the valid ones.
## For now I think it's best to leave it like this
base_pair[indel_length == 1] = simplify_ctx(base_pair[indel_length == 1])

## now I compute the second feature, REPEAT SIZE
indel_ampli_set = nchar(indel_affected_seq) * params_maxRep

# this is assuming left-aligned because it only looks to the right.
# (I think this is okay, but maybe we should check)
# the specific number depends if indel is ins or deletion
upstream_expansion
vr_indel_ex = helperMut::extend(x = vr_indel,
                              upstream = -nchar(ref(vr_indel)),
                              downstream = indel_ampli_set)

flank_seq = getSeq(params_genome,vr_indel_ex)

split_pattern = glue::glue("(?<=.{<<nchar(indel_affected_seq)>>})",
                           .open = "<<",.close = ">>")
flank_seq_spl = strsplit(as.character(flank_seq), split_pattern, perl = TRUE)

purrr::map2_dbl(.x = flank_seq_spl,
            .y = indel_affected_seq,
            function(flank,indel){
              vec_rle = Rle(indel == flank)
              if (vec_rle@values[1]){
                rpt_size = vec_rle@lengths[1]
              } else {
                rpt_size = 0
              }
              return(rpt_size)
}) -> indel_rpt_size


table(indel_rpt_size,base_pair,indel_type) %>% as.data.frame()






# indels can be left or right aligned (I think)
indels_clssifier <- function(vr,genome = genome_selector()){

}

