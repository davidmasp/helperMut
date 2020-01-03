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




### LOCAL RUN TEST!
# library(helperMut)
# fn = Sys.getenv("HARTWIG_example")
# dat_vcf = VariantAnnotation::readVcf(fn)
# dat_vr = as(object = dat_vcf,Class = "VRanges")
# dat_indels = dat_vr[dat_vr$set == "indels"]
# genome = helperMut::genome_selector(alias = "Hsapiens.1000genomes.hs37d5")
# indels_classifier(vr = dat_indels,genome = genome)

## this is assuming that the first letter in REF is shared with ALT and
## corresponds to the refseq
# indels NEED TO BE left-aligned
indels_classifier <- function(vr,
                             maxRep = 5,
                             genome = genome_selector(),
                             check_format = TRUE,
                             groupVar = "group"){

  if (length(groupVar) == 1){
    groupVar = rep(groupVar,length(vr))
  } else if (length(groupVar) != length(vr)){
    stop("wrong group var inputed")
  }

  mcols(vr) = NULL

  if(check_format){
    refSeq = Biostrings::getSeq(genome,vr)
    stopifnot(refSeq == VariantAnnotation::ref(vr))
  }

  ref_l = nchar(VariantAnnotation::ref(vr))
  alt_l = nchar(VariantAnnotation::alt(vr))

  stopifnot(!any(ref_l == alt_l))
  ins_mask = ref_l < alt_l

  vr_ins = vr[ins_mask]
  vr_del = vr[!ins_mask]


  # this is a bit stupid but how my brain understands this
  vr_indel = c(vr_ins,vr_del)
  indel_type = c(rep("insertion",length(vr_ins)),
                 rep("deletion", length(vr_del)))

  ### INSERTIONS ###
  ins_affected_seq = stringr::str_sub(VariantAnnotation::alt(vr_ins),start = 2)
  ins_length = nchar(ins_affected_seq)

  ### DELETIONS ###
  del_affected_seq = stringr::str_sub(VariantAnnotation::ref(vr_del),start = 2)
  del_lenght = nchar(del_affected_seq)

  ## MICROHOMOLOGY ##
  # mh can only be computed in +2bp
  candidate_mh = del_lenght >= 2
  vr_del_mh = extend(x = vr_del,
                     upstream = del_lenght-1,
                     downstream = del_lenght)
  seq_vec = as.character(Biostrings::getSeq(genome,vr_del_mh))
  mh_det = detect_microhomology(seq_vector = seq_vec)
  mh_det = apply(mh_det,1, max)
  actually_mh = candidate_mh & mh_det > 0
  ## this values will be uesd later to over-right the repeat values
  ## from when computing the repeats. It is important to check that their
  ## repeat length has to be 0!

  ## JOIN, same order than the VR!!! ##
  indel_affected_seq = c(ins_affected_seq,del_affected_seq)
  indel_length = c(ins_length,del_lenght)
  # this is also important as mh is only computed in deletions
  actually_mh = c(rep(FALSE,length(ins_affected_seq)),actually_mh)
  mh_det = c(rep(NA,length(ins_affected_seq)),mh_det)

  ## first I compute the first feature the INSERTION SIZE
  base_pair = dplyr::case_when(
    indel_length == 1 ~ indel_affected_seq,
    indel_length >= maxRep ~ as.character(glue::glue("{maxRep}+bp")),
    TRUE ~ as.character(glue::glue("{indel_length}bp"))
  )

  ## this don't work directly in the case_when because it evaluates all
  ## the elements in the vector and then extracts the valid ones.
  ## For now I think it's best to leave it like this
  base_pair[indel_length == 1] = simplify_ctx(base_pair[indel_length == 1])

  ## now I compute the second feature, REPEAT SIZE
  indel_ampli_set = nchar(indel_affected_seq) * maxRep

  # this is assuming left-aligned because it only looks to the right.
  # (I think this is okay, but maybe we should check)
  # the specific number depends if indel is ins or deletion
  # upstream_expansion
  vr_indel_ex = extend(x = vr_indel,
                                  upstream = -nchar(VariantAnnotation::ref(vr_indel)),
                                  downstream = indel_ampli_set)

  flank_seq = Biostrings::getSeq(genome,vr_indel_ex)

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


  ## adding mh info,
  ## A repeat will be catched as MH
  to_change_mask = actually_mh & indel_rpt_size == 0

  ## I also can change to strings now
  indel_rpt_size = ifelse(indel_rpt_size >= maxRep,
                          yes =  glue::glue("+{maxRep}"),
                          no = glue::glue("{indel_rpt_size}"))

  mh_det_str =  ifelse(mh_det >= maxRep,
                       yes =  glue::glue("+{maxRep}"),
                       no = glue::glue("{mh_det}"))

  indel_rpt_size[to_change_mask] = mh_det_str[to_change_mask]
  indel_type[to_change_mask] = "deletion-mh"

  table(indel_rpt_size = factor(indel_rpt_size),
        base_pair = factor(base_pair),
        indel_type = factor(indel_type),
        group = groupVar) -> res_table

  as.data.frame(res_table) -> res

  # impossible cases
  # I have to use factors to get 0 in the table step but then this generates
  # 0 in impossible cases that have to be removed

  mh_1 = res$indel_type == "deletion-mh" & (res$base_pair %in% c("A","C"))

  irs = as.numeric(stringr::str_extract(res$indel_rpt_size,"[:digit:]+"))
  baseP = as.numeric(stringr::str_extract(res$base_pair,"[:digit:]+"))
  baseP[is.na(baseP)] = 1
  mh_2 =  res$indel_type == "deletion-mh" & (irs >= baseP) & (baseP != maxRep)
  mh_3 = res$indel_type == "deletion-mh" &  irs == 0

  impossible_vec = res[mh_1 | mh_2 | mh_3,][["Freq"]]
  stopifnot(all(impossible_vec == 0))
  res = res[!(mh_1 | mh_2 | mh_3),]
  rownames(res) = NULL

  res
}




