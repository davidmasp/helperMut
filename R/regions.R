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


#' Obtain genomic regions from gene lists
#'
#' Wrapper for a workflow that obtains a given genomic feature
#' from some gene identifier.
#'
#' Basically it comprises 2 steps, first, transform input
#' gene ids (such as SYMBOL) to the gene id type for the
#' txdb object. Second, it extracts the given feature filtering
#' for those genes inputed.
#'
#' @param feat the genomic feature to extract from (such as genes, transcripts...) see ?GenomicFeatures::transcripts
#' @param annotation_pkg An OrgDb object
#' @param transcript_db A TxDb object
#' @param identifiers the genes that we want to retrieve
#' @param indentifier_type the input identifier column
#' @param transcript_db_id_type the identifier type for the txdb object
#' @param transcript_db_columns extra columns to retrieve from the txdb object
#' @param ... extra options for feat function
#'
#' @export
#' @examples
#'
#' \dontrun{
#' # retrieving from ensembl
#' human = org.Hs.eg.db::org.Hs.eg.db
#' txdb_ensembl = makeTxDbFromEnsembl(organism = "homo_sapiens")
#' obtain_genomic_feature(feat = "transcripts",
#'                        annotation_pkg = human,
#'                        transcript_db = txdb_ensembl,
#'                        identifiers = c("BRCA1","BRCA2"),
#'                        indentifier_type = "SYMBOL",
#'                        transcript_db_id_type = "ENSEMBL",
#'                        transcript_db_columns = "TXNAME")
#' }
#'
#' # retrieving from UCSC (faster)
#'
#' human = org.Hs.eg.db::org.Hs.eg.db
#' txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
#' obtain_genomic_feature(feat = "transcripts",
#'                        annotation_pkg = human,
#'                        transcript_db = txdb,
#'                        identifiers = c("BRCA1","BRCA2"),
#'                        indentifier_type = "SYMBOL",
#'                        transcript_db_columns = "TXNAME")
#'
obtain_genomic_feature <- function(feat,
                                   annotation_pkg,
                                   transcript_db,
                                   identifiers,
                                   indentifier_type,
                                   transcript_db_id_type =  "ENTREZID",
                                   transcript_db_columns = NULL,
                                   ...){
  ## because idk about more annotations than txdb I think for human at least
  ## we are save moving to ENTTREZID which is then used to retrieve the
  ## cds / whatever
  if (indentifier_type != transcript_db_id_type ){
    AnnotationDbi::mapIds(x = annotation_pkg,
                          keys = identifiers,
                          keytype = indentifier_type,
                          column = transcript_db_id_type) -> txdb_ready_ids
  } else {
    txdb_ready_ids = identifiers
  }

  # will stop if there are missing genes
  stopifnot(all(!is.na(txdb_ready_ids)))

  ## Here I select the feature I want to extract
  fun = getFromNamespace(feat,ns = "GenomicFeatures")

  filtering_list = list(gene_id = txdb_ready_ids)
  fun(x = transcript_db,
      filter = filtering_list,
      columns = transcript_db_columns,...) -> res_gr

  return(res_gr)

}


