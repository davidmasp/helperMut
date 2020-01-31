#                               $$\
#                               \__|
#  $$$$$$\   $$$$$$\   $$$$$$\  $$\  $$$$$$\  $$$$$$$\   $$$$$$$\
# $$  __$$\ $$  __$$\ $$  __$$\ $$ |$$  __$$\ $$  __$$\ $$  _____|
# $$ |  \__|$$$$$$$$ |$$ /  $$ |$$ |$$ /  $$ |$$ |  $$ |\$$$$$$\
# $$ |      $$   ____|$$ |  $$ |$$ |$$ |  $$ |$$ |  $$ | \____$$\
# $$ |      \$$$$$$$\ \$$$$$$$ |$$ |\$$$$$$  |$$ |  $$ |$$$$$$$  |
# \__|       \_______| \____$$ |\__| \______/ \__|  \__|\_______/
#                     $$\   $$ |
#                     \$$$$$$  |
#                      \______/



#' Shufle BED
#'
#' Similar to \href{https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html}{shufleBed} from bedtools.
#'
#'
#' @param regions regions object to be shufled (Granges)
#' @param mask mask regions where to shufle (Granges)
#'
#' @return
#' The function returns a GR object with the regions located in new random
#' positions. The index of the output should match index of the input.
#'
#' Chromosome / seqnames are not maintained.
#'
#' Resulting regions may overlap
#'
#' @export
#'
#' @examples
#'
shufle_regions <- function(regions,mask) {
  # hereI should check that regions have a unique match in the mask
  # although I am not sure if I care

  gr_reg = regions
  mask_gr = mask

  # here I compute the width of the region to shuffle
  width_reg_vec = BiocGenerics::width(gr_reg)
  width_mask_vec = BiocGenerics::width(mask_gr)

  # now I choose to which instance in the mask the new region will be sent
  # chromosomes / seqnames are not maintained! see split
  mask_instance = integer(length = length(gr_reg))


  for (i in seq(from = 1, to = length(gr_reg),by = 1)){

    width_r = width_reg_vec[i]
    mask_of_masks = width_r < width_mask_vec
    # here I compute a unified probability, False will be 0!
    probs = as.numeric(mask_of_masks) / sum(mask_of_masks)

    # first sampling, instance in the mask
    res = sample(x = length(mask_gr),size = 1,prob = probs)
    mask_instance[i] = res
    # so here the problem is that probs which is the width mask needs to
    # be defined per instance so then I cannot vectorize this fully.
  }

  # now I sample which window instance will every region be sent, the max value
  # is the one for the widest instance in the mask
  max_value = max(width_mask_vec)

  # I sample from 1 to this number
  sampled_values = sample(max_value,size = length(gr_reg),replace = TRUE)

  # each gr has been asigned to a new instance so I compute the width of the
  # new instance
  window_width_values = width_mask_vec[mask_instance]

  # the max start which is equal to number of possible bins that you will find
  # in a instance is equal to the width of the mask instanec minus the lenght
  # of the region. These are tiling windows so its the same as the width
  max_start =  window_width_values - width_reg_vec

  # because we sampled with the max value in the whole set of instances, this
  # value can be bigger that the number of possible instances to solve this
  # we compute the residual of each division. Then we will get the instance
  # value if we filled all the possibles ones before.
  wind_idx = sampled_values %% max_start


  # generate the output

  # now I select the seqname of the mask
  seqnames_rnd = GenomeInfoDb::seqnames(mask_gr)[mask_instance]

  # now I select the start
  start_rnd = BiocGenerics::start(mask_gr)[mask_instance] + wind_idx

  # I create the output
  random_gr = GenomicRanges::GRanges(
    seqnames = seqnames_rnd,
    ranges = IRanges::IRanges(start = start_rnd,
                     width = width_reg_vec)
  )

  return(random_gr)

}







## regions


#' Shufle BED
#'
#' Similar to \href{https://bedtools.readthedocs.io/en/latest/content/tools/shuffle.html}{shufleBed} from bedtools.
#'
#'
#' @param regions regions object to be shufled (Granges)
#' @param mask mask regions where to shufle (Granges)
#'
#' @return
#' The function returns a GR object with the regions located in new random
#' positions. The index of the output should match index of the input.
#'
#' Chromosome / seqnames are not maintained.
#'
#' Resulting regions may overlap
#'
#' @export
#'
#' @examples
#'
shufle_regions2 <- function(regions,mask) {
  # hereI should check that regions have a unique match in the mask
  # although I am not sure if I care

  gr_reg = regions
  mask_gr = mask

  # here I compute the width of the region to shuffle
  width_reg_vec = BiocGenerics::width(gr_reg)
  width_mask_vec = BiocGenerics::width(mask_gr)

  # now I choose to which instance in the mask the new region will be sent
  # chromosomes / seqnames are not maintained! see split
  mask_instance = sample(x = length(mask_gr),
                         size = length(gr_reg),
                         replace = TRUE)

  wider_mask = width_reg_vec > width_mask_vec[mask_instance]

  if (sum(wider_mask) > 0){
    for (i in which(wider_mask)){
      width_r = width_reg_vec[i]
      mask_of_masks = width_r < width_mask_vec
      # if there is no place for this region to shufle, move it!
      stopifnot(sum(mask_of_masks) != 0)
      # here I compute a unified probability, False will be 0!
      probs = as.numeric(mask_of_masks) / sum(mask_of_masks)
      # first sampling, instance in the mask
      res = sample(x = length(mask_gr),size = 1,prob = probs)
      mask_instance[i] = res
      # so here the problem is that probs which is the width mask needs to
      # be defined per instance so then I cannot vectorize this fully.
    }
  }

  # now I sample which window instance will every region be sent, the max value
  # is the one for the widest instance in the mask
  max_value = max(width_mask_vec)

  # I sample from 1 to this number
  sampled_values = sample(max_value,size = length(gr_reg),replace = TRUE)

  # each gr has been asigned to a new instance so I compute the width of the
  # new instance
  window_width_values = width_mask_vec[mask_instance]

  # the max start which is equal to number of possible bins that you will find
  # in a instance is equal to the width of the mask instanec minus the lenght
  # of the region. These are tiling windows so its the same as the width
  max_start =  window_width_values - width_reg_vec

  # because we sampled with the max value in the whole set of instances, this
  # value can be bigger that the number of possible instances to solve this
  # we compute the residual of each division. Then we will get the instance
  # value if we filled all the possibles ones before.
  wind_idx = sampled_values %% max_start


  # generate the output

  # now I select the seqname of the mask
  seqnames_rnd = GenomeInfoDb::seqnames(mask_gr)[mask_instance]

  # now I select the start
  start_rnd = BiocGenerics::start(mask_gr)[mask_instance] + wind_idx

  # I create the output
  random_gr = GenomicRanges::GRanges(
    seqnames = seqnames_rnd,
    ranges = IRanges::IRanges(start = start_rnd,
                              width = width_reg_vec)
  )

  return(random_gr)

}

