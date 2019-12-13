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


