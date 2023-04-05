# $$\      $$\ $$\   $$\ $$$$$$$$\  $$$$$$\
# $$$\    $$$ |$$ |  $$ |\__$$  __|$$  __$$\
# $$$$\  $$$$ |$$ |  $$ |   $$ |   $$ /  \__|
# $$\$$\$$ $$ |$$ |  $$ |   $$ |   \$$$$$$\
# $$ \$$$  $$ |$$ |  $$ |   $$ |    \____$$\
# $$ |\$  /$$ |$$ |  $$ |   $$ |   $$\   $$ |
# $$ | \_/ $$ |\$$$$$$  |   $$ |   \$$$$$$  |
# \__|     \__| \______/    \__|    \______/



#' GRanges extension
#'
#' @description
#' The function extend all the ranges in a \link[GenomicRanges]{GRanges} object by the number
#' of nucleotides defined in the upstream or downstream of each range.
#'
#' It also works for \link[VariantAnnotation]{VRanges} objects.
#'
#' extracted from \url{https://support.bioconductor.org/p/78652/}
#
#'
#' @param x A GRanges or VRanges object
#' @param upstream Number of nucleotides to extend upstream-wise
#' @param downstream Number of nucleotides to extend downstream-wise
#'
#' @return
#'
#' A GRanges or VRanges object
#'
#' @export
#'
#' @examples
extend <- function(x, upstream, downstream) {

  # I have to do this because they don't let me use the strand in a VRanges
  # object, therefore we need to translate to GR before extending
  # see here -> https://support.bioconductor.org/p/127590/
  if (is(object = x,class2 = "VRanges")){
    x = as(object = x,Class = "GRanges")
  }

  if (any(BiocGenerics::strand(x) == "*")){
    warning("'*' ranges were treated as '+'")
  }

  on_plus <- BiocGenerics::strand(x) == "+" | BiocGenerics::strand(x) == "*"
  new_start <- BiocGenerics::start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- BiocGenerics::end(x) + ifelse(on_plus, downstream, upstream)
  IRanges::ranges(x) <- IRanges::IRanges(new_start, new_end)
  x = IRanges::trim(x) # this removes the negative Ranges
  return(x)
}


#' get mutation subtypes from VR
#'
#' @param x a VRanges object
#' @param sep a separator string to divide context and alternative allele
#' @param k number of nucleotides to expand the context in one direction
#' @param genome a BSgenome object
#' @param simplify_set a set of DNA bases which to simplify the mutation strings
#' @param keep_strand true if no simplification is needed
#'
#' @return
#'
#' A string vector with mutation strings corresponding to each position in the
#' VRanges object. Should have the same length as the original inputed
#' object.
#'
#' @export
#'
#' @examples
#'
#' library(VariantAnnotation)
#' library(helperMut)
#' genome = genome_selector("Hsapiens.1000genomes.hs37d5")
#' vr = VRanges(seqnames = c(1,2,3),
#'             ranges = IRanges(start = c(1e6,1e6,1e6),width = 1),
#'             ref = "N", alt = c("A","T","A"))
#' get_MS_VR(vr,genome = genome)
#' get_MS_VR(vr,genome = genome,sep = "_")
#' get_MS_VR(vr,genome = genome,sep = "_",keep_strand = TRUE)
#' get_MS_VR(vr,genome = genome,simplify_set = c("C","T"))
#'
get_MS_VR <- function(x,
                      sep=">",
                      k = 1,
                      genome = genome_selector(),
                      simplify_set = c("C","A"),
                      keep_strand = FALSE){

  # general objects needed
  requireNamespace("stringr",quietly = TRUE)
  requireNamespace("Biostrings",quietly = TRUE)
  requireNamespace("VariantAnnotation",quietly = TRUE)
  requireNamespace("GenomeInfoDb",quietly = TRUE)

  if (is(object = genome,class2 = "BSgenome")) {
    genome = genome
  } else {
    stop("Enter a valid genome object or a character")
  }


  if (!any(GenomeInfoDb::seqlevelsStyle(x) %in%
           GenomeInfoDb::seqlevelsStyle(genome))){
    stop("Seq levels do not match in regions and genome, check the reference")
  }

  tmp <- extend(x,upstream = k,downstream = k)
  seq <-  Biostrings::getSeq(genome,tmp)

  # from here can substitute for simplify_muts function
  alt_seq = VariantAnnotation::alt(x)

  ## if the strand is minus, I need to switch the alt alelle.
  if (any(strand(x) == "-")){
    warning("feature in development")
    alt_seq_dnastringset = Biostrings::DNAStringSet(alt_seq)
    alt_seq_rc = Biostrings::reverseComplement(alt_seq_dnastringset)
    alt_seq = ifelse(test = strand(x) == "-",
                     yes = alt_seq_rc,
                     no = alt_seq)
  }


  #browser()
  muts = paste(seq,alt_seq,sep = sep)

  if (!keep_strand){
    muts = simplify_muts(muts,simplify_set = simplify_set,sep = sep)
  }

  return(muts)
}

#' @describeIn str_reverse_complement Rc of a mutation string
ms_reverse_complement <- function(ms){
  # assumes sep is nchar 1

  nchar_ms = unique(nchar(ms))
  K = nchar_ms - 2

  sep = unique(substr(ms,nchar_ms-1,nchar_ms-1))
  alt = substr(ms,nchar_ms,nchar_ms)
  ctx = substr(ms,1,K)

  ctx_rc = str_reverse_complement(str = ctx)
  alt_rc = str_reverse_complement(str = alt)

  paste(ctx_rc,alt_rc,sep = sep)
}

#' Reverse complement of a string
#'
#' It returns the reverse complement of a context or mutation string.
#'
#' Internally it transforms the string to DNAstringSet from Biostrings.
#'
#' @param str a sequence character vector
#' @param ms a mutation character vector (assumes 1 separation character)
#'
#' @return A string with the reverseComplement of the input
#' @export
#'
#' @examples
#'
#' str_reverse_complement(c("CCC","TTT"))
#'
#' ms_reverse_complement(c("CCC>T","TTT>A"))
#'
str_reverse_complement <- function(str){
  dna_set = Biostrings::DNAStringSet(str)
  dna_set_rc = Biostrings::reverseComplement(dna_set)
  as.character(dna_set_rc)
}


#' Count mutaiton subtypes
#'
#' Given a mutation subtype string vector, automatically indentify the
#' parameters and generate a interger vector with counts.
#'
#' @param ms A mutation subtype string vector
#'
#' @return
#' a integer vector with frequency counts
#' @export
#'
#' @examples
#'
#' count_MS(c("TCA>T","CAC>T"))
#' count_MS(c("TCA_T","CAC_T"))
count_MS <- function(ms) {

  mutaes = identify_mut_aestetics(ms)

  k = mutaes[["k"]]
  simplify_set = mutaes[["mset_ref"]]
  sep = mutaes[["sep"]]

  if (length(simplify_set) > 2) {
    stop("Mutation type input is not simplified, stopping")
  } else if (length(simplify_set) == 2 &&
             simplify_set[1] == str_reverse_complement(simplify_set[2])){
    stop("Mutation type input is not simplified, stopping")
  } else if (length(simplify_set) == 1){
    mssg =
      "A unique central position is presented, returning just those subtypes"
    warning(mssg)
  }

  mut_types = generate_mut_types(k=k,
                                 sep=sep,
                                 simplify_set = simplify_set)

  if (!all(ms %in% mut_types)){
    stop("Some mutations are not standard. Mutations should be simplified")
  }

  ms_fct = factor(ms,levels = mut_types)

  row_counts = (table(ms_fct))
  row_vec = as.vector(row_counts)
  names(row_vec) = names(row_counts)


  # this should be fine if this has passed the previous test
  stopifnot(sum(row_counts) == length(ms))

  return(row_counts)
}


#' Simplify mutations
#'
#' @description
#' Given a set of mutations it removes the strand information and simplifies them.
#' It also translate a mutation set to an alternative refernce set (for instance, C T based to C A based.)
#'
#' @param muts A character vector containing mutations
#' @param ctx A character vector containing contexts
#' @param simplify_set The reference base set. (must be 2)
#' @param sep The string that separates context and alternative
#'
#' @return A set of simplified mutations
#' @export
#'
#' @examples
#'
#' mutset = c("ACT>T","TTA>C")
#' simplify_muts(mutset)
#'
#' ctxset = c("ACT","TTA")
#' simplify_ctx(ctxset)
#'
simplify_muts <- function(muts,simplify_set = c("C","A"),sep=">") {

  if (length(muts) == 0){
    stop("Muts is empty")
  }

  # this is not true from multiple length separators
  sep_nchar = nchar(sep)
  k = unique((nchar(muts) - 2 - sep_nchar)/2)
  K = (k * 2) + 1
  stopifnot(length(k)==1)

  muts_ctx = stringr::str_sub(string = muts,
                              start = 1,
                              end = K
  )
  muts_alt = stringr::str_sub(string = muts,
                              start = K+1+sep_nchar, # this is the end position
                              end = K+1+sep_nchar
  )

  seq = Biostrings::DNAStringSet(muts_ctx)
  alt_seq = Biostrings::DNAStringSet(muts_alt)
  # important to do reverse complementary
  # could add a method to mantain strand
  cp = k + 1

  central_pos = stringr::str_sub(as.character(seq),start = cp,end= cp)

  # first we add the ones in the positive set
  simplify_pattern = glue::glue("[{paste0(simplify_set,collapse='')}]")
  simplify_mask = grepl(pattern = simplify_pattern,
                        x = central_pos)

  # then we add the complementary
  complementary_set = as.character(
    Biostrings::reverseComplement(
      Biostrings::DNAStringSet(
        simplify_set
      )))
  complementary_pattern = glue::glue("[{paste0(complementary_set,collapse='')}]")
  complement_mask = grepl(pattern = complementary_pattern,
                          x = central_pos)

  reverse_mask = !simplify_mask
  seq[reverse_mask] = Biostrings::reverseComplement(seq[reverse_mask])
  alt_seq[reverse_mask] = Biostrings::reverseComplement(alt_seq[reverse_mask])

  seq = as.character(seq)
  alt_res = as.character(alt_seq)

  res = paste(seq,alt_res,sep = sep)

  # apply the NA mask
  if(any(!simplify_mask != complement_mask)){
    warning("Some mutations contain out of set bases. Setting NAs")
    na_mask = !simplify_mask & !complement_mask
    # this means that either the cp is not in the normal set neither the
    # complementary set, then this goes to NA (see issue #25)
    res[na_mask] = NA
  }

  return(res)
}

#' @describeIn simplify_muts Simplify contexts
simplify_ctx <- function(ctx,simplify_set = c("C","A")) {

  if (length(ctx) == 0){
    stop("CTX is empty")
  }

  ###### DANGER!! ###########
  k = unique((nchar(ctx) - 1)/2)
  stopifnot(length(k)==1)

  muts_ctx = ctx

  seq = Biostrings::DNAStringSet(muts_ctx)
  # important to do reverse complementary
  # could add a method to mantain strand
  cp = k + 1

  central_pos = stringr::str_sub(as.character(seq),start = cp,end= cp)

  # first we add the ones in the positive set
  simplify_pattern = glue::glue("[{paste0(simplify_set,collapse='')}]")
  simplify_mask = grepl(pattern = simplify_pattern,x = central_pos)

  # then we add the complementary
  complementary_set = as.character(
    Biostrings::reverseComplement(
      Biostrings::DNAStringSet(
        simplify_set
      )))
  complementary_pattern = glue::glue("[{paste0(complementary_set,collapse='')}]")
  complement_mask = grepl(pattern = complementary_pattern,
                          x = central_pos)

  reverse_mask = !simplify_mask
  seq[reverse_mask] = Biostrings::reverseComplement(seq[reverse_mask])


  seq = as.character(seq)

  res = seq

  # apply the NA mask
  if(any(!simplify_mask != complement_mask)){
    warning("Some mutations contain out of set bases. Setting NAs")
    na_mask = !simplify_mask & !complement_mask
    # this means that either the cp is not in the normal set neither the
    # complementary set, then this goes to NA (see issue #25)
    res[na_mask] = NA
  }

  return(res)
}

#' Make mutation set
#'
#' @param x a simplified version of a mutation set in UPAC format or a vector of non-overlaping sets.
#' @param simplify Boolean to simplify or not the result
#' @param simplify_set Simpliy reference set if simplify is true.
#' @param sep string to separate the context and the alternative allele
#' @param strict removes impossible mutation types from the output
#'
#' @return a mutation set with all the possible mutation types encoded.
#' @export
#'
#' @examples
#'
#' make_set("WAW>T")
#'
make_set <- function(x,
                     simplify = TRUE,
                     simplify_set = c("C","A"),
                     sep=">",
                     strict = TRUE){

  if (length(x)>1){
    lapply(X = x,
           FUN = make_set,
           simplify = simplify,
           simplify_set = simplify_set,
           sep=sep,
           strict = strict) -> sets_list
    return(unlist(sets_list))
  }

  k = (nchar(x) - 3)/2

  a = substr(x,1,k)
  b = substr(x,k+1,k+1) # this is the central position
  e = substr(x,k+2,(k*2)+1)

  to = substr(x,(k*2)+3,(k*2)+3)

  input_to_grid = list()

  sa = 1:nchar(a) %>%
    purrr::map(function(x){
    substr(x = a,start = x,stop = x)
    }) %>%
    purrr::map(function(x){
      helperMut::dna_codes[[x]]
    })

  # yes. but is should be a list.
  sb = list(helperMut::dna_codes[[b]])

  se = 1:nchar(e) %>%
    purrr::map(function(x){
      substr(x = e,start = x,stop = x)
    }) %>%
    purrr::map(function(x){
      helperMut::dna_codes[[x]]
    })

  # this should be a vector because how expand.grid works.
  sto = helperMut::dna_codes[[to]]

  comb_df = expand.grid(c(sa,sb,se))

  from = apply(comb_df, 1, paste0,collapse="")

  comb_df = expand.grid(from,sto)
  if (strict){
    same_mask = comb_df$Var2 == stringr::str_sub(comb_df$Var1,k+1,k+1)
    if (sum(same_mask) > 0){
      warning("Context set contained impossible combinations, removing them from output because strict = TRUE")
      comb_df = comb_df[!same_mask,]
    }
  }

  ms_set = apply(comb_df, 1, paste0,collapse=">")

  if (simplify){
    ms_set = simplify_muts(ms_set,
                           simplify_set = simplify_set,
                           sep = sep )
  }

  return(ms_set)
}



#' Compute test for mutation enrichment in a region
#'
#' This functions computes a fisher test to check for an enrichment
#' or depletion based on the size of the region tested and the genome
#' or genomic mask.
#'
#' It does not control for sequence composition, just for size.
#'
#' Warning: target regions will be intercepted with mask if not NULL
#'
#' @param vr a VRanges object with mutations
#' @param gr a target region to test for enrichment, a GRanges object
#' @param genome a BSgenome object
#' @param genome_mask a regions mask that contains the target regions to test.
#' @param ... arguments for the fisher test
#'
#' @return
#' a tidy dataframe from broom::tidy.
#' @export
#'
#' @examples
#'
#' library(VariantAnnotation)
#' genome = genome_selector("Hsapiens.UCSC.hg19")
#' base_pos = 6e4
#' vr_target = VRanges(seqnames="chr1",
#'                     ranges=IRanges(c(base_pos + 6, base_pos + 16),
#'                                    width = 1),
#'                     ref = "C",alt = "A")
#' gr_target <- GRanges(seqnames="chr1",
#'                      ranges=IRanges(base_pos + 3,base_pos +  10))
#' gr_mask <- GRanges(seqnames=c("chr1", "chr1"),
#'                    ranges=IRanges(c(base_pos + 4,base_pos + 15),
#'                                   width = 10))
#'
#' mutation_enrichment_general(vr = vr_target,
#'                             gr = gr_target,
#'                             genome = genome,
#'                             genome_mask = gr_mask)
mutation_enrichment_general <- function(vr,
                                         gr,
                                         genome = genome_selector(),
                                         genome_mask = NULL,
                                        ...) {

  if (is.null(genome_mask)){
    gr_target = gr
  } else {
    # [^mdfkjs]
    gr_target = GenomicRanges::intersect(gr,genome_mask)
  }

  ovr = GenomicRanges::findOverlaps(query = vr,subject = gr_target)
  ctx_cont = get_k_freq_fromRegion(gr = gr_target, k = 0, genome = genome)

  if (is.null(genome_mask)){
    GenomeInfoDb::seqlevels(gr) = GenomeInfoDb::seqlevels(genome)
    GenomeInfoDb::seqlengths(gr) = GenomeInfoDb::seqlengths(genome)
    GenomeInfoDb::genome(gr) = GenomeInfoDb::genome(genome)
    gr_gaps = gaps(gr)
    gr_gaps = gr_gaps[strand(gr_gaps) == "*"]
  } else {
    # in this scenario gr_target should be always included in genome_mask
    # see [^mdfkjs]
    gr_gaps = GenomicRanges::setdiff(x = genome_mask,y = gr_target)
  }

  ctx_cont_gaps = get_k_freq_fromRegion(gr = gr_gaps,
                                        k = 0,
                                        genome = genome)

  ovr_out = GenomicRanges::findOverlaps(query = vr,subject = gr_gaps)

  muts_in = length(ovr)
  muts_out = length(ovr_out)
  ctx_in = sum(ctx_cont)
  ctx_out = sum(ctx_cont_gaps)

  cnt_table = matrix(c(muts_in,
                       ctx_in,
                       muts_out,
                       ctx_out),
                     byrow = TRUE,
                     nrow = 2,
                     dimnames = list(Mutation = c("Mut","wt"),
                                     Region = c("in","out")))

  res_df = broom::tidy(fisher.test(cnt_table,...))

  res_df
}

#' Compute mutational subtype matrix
#'
#' @description
#' It takes a SNVs VR object and it outputs the matrix with the mutation
#' counts for each mutation subtype. Samples are in rows, subtypes in
#' in columns.
#'
#' see also \link[SomaticSignatures]{motifMatrix}
#'
#' @param vr A VRanges object with SNV calls.
#' @param k The extension of the mutation subtype. This is equivalent to the positions upstream and downstream to be included, k = 1 for trinucleotydes, k = 2 for pentanucleotides.
#' @param sep A string to separate the mutation context to the alternative base.
#' @param genome A BSgenome object installed in the local computer.
#' @param simplify_set Set of base pairs used to simplify the mutation calls.
#'
#' @return A count matrix with samples in rows and mutation subtypes in columns
#' @export
#'
#' @examples
compute_MSM_fast = function(vr,
                            k=1,
                            sep = ">",
                            genome = genome_selector(),
                            simplify_set = c("C","A")) {

  # based on the implementation of somaticSignatures::motifMatrix

  # needs total positive implementation
  stopifnot(requireNamespace("VariantAnnotation",quietly = T))
  vr$ms = get_MS_VR(x = vr,
                     sep = sep,
                     k = k,
                     genome = genome,
                     simplify_set = simplify_set)

  group = sampleNames(vr)
  motif = factor(vr$ms,
                 levels = generate_mut_types(k = k,
                                             sep = sep,
                                             simplify_set = simplify_set))

  y = as(table(motif, group), "matrix")
  dimnames(y) = unname(dimnames(y))
  res = t(y)
  stopifnot( sum(res) == length(vr)) # safety check against issue #22
  return(res)
}


#' Correct MSM per total positives
#'
#' @param msm MSM matrix from compute_MSM_fast
#' @param P vector with the positives per row
#'
#' @export
#'
#' @examples
correct_MSM_P <- function(msm,P){
  freqs = t(apply(msm, 1, function(row){
    row/sum(row)
  }))

  res = freqs * P

  return(res)

}


#' Detect mutation subtype params
#'
#' From a single string of mutations detect which parameters, such as k and
#' sep are being used.
#'
#' It also identifies the reference set of bases used as central position.s
#'
#' @param ms A mutation subtype version (any k)
#' @param force Avoid error if non-standard mutations are used as a central position, such as TNA>C
#'
#' @return A list with the extracted parameters
#' @export
#'
#' @examples
#'
#' helperMut::identify_mut_aestetics(c("TCA>T","CCA>T","TTT>A"))
#'
identify_mut_aestetics = function(ms, force = FALSE){
  l = unique(nchar(ms))
  stopifnot(length(l)==1)

  k = (l - 3 )/ 2
  K = (k * 2) + 1

  mset_ref = stringr::str_sub(ms,k + 1,k + 1)

  cond = mset_ref %in% unique(unlist(helperMut::dna_codes))
  if (any(!cond) & !force){ #shit
    idx = which(!cond)
    idxs = paste0(idx,collapse = ",")
    stop(glue::glue("Non standard mutations included, check index {idxs}"))
  }

  sep = stringr::str_sub(ms,l-1, l-1)

  res = list(
    k = k,
    K = K,
    mset_ref = unique(mset_ref),
    sep = unique(sep)
  )

  return(res)
}


#' Generate Mutation Subtypes
#'
#' @description
#' Generate all the possible mutation subtypes in a mutation set.
#'
#' @param k The number of upstream and downstream bases the motif is generated
#' @param sep The string that separates the motif to the alternative allele.
#' @param simplify_set The reference base set to base the mutation subtypes
#'
#' @return A complete mutation set
#' @export
#'
#' @examples
#' trincl = generate_mut_types(k = 1)
#' length(trincl) # should be 96
#'
#' pentancl = generate_mut_types(k = 2)
#' length(pentancl) # should be 1536
#'
generate_mut_types <- function(k,
                               sep=">",
                               simplify_set = c("C","A")){


  panels = purrr::map(simplify_set,function(bp){
    pe = c("A","C","G","T")
    pc = bp
    tomut = pe[! pe == bp]

    grid_input = list()
    for (i in 1:k){grid_input[[i]] = pe}
    grid_input[[(k+1)]] = pc
    for (i in (k+2):(k*2+1)){grid_input[[i]] = pe}
    grid_input[[(k*2+2)]] = tomut
    lp = expand.grid(grid_input)
    return(lp)
  })

  panel = do.call("rbind",panels)

  step1 = tidyr::unite(data = panel,
                       tidyselect::num_range(prefix = "Var",
                                             range = (k*2+1):1),
                       col = "ctx",
                       sep="")

  step2 = tidyr::unite(data=step1,
                       ctx,
                       tidyselect::num_range(prefix = "Var",
                                             range = (k*2+2)),
                       col = "pos",
                       sep = sep)

  return(as.character(step2$pos))
}


# transforms any mutation level to ms6 or k=1
#' Muts to MS6
#'
#' Transforms any mutation level (any k) to ms6 (a.k.a k=0)
#'
#' @param mutation_vec a string vector with mutation strings
#' @param sep a separator used in the mutation string
#' @return Returns a vector of same length with the mutation types at k=0
#' @export
#' @examples
#' mut2ms6("TCA>T") # returns C>T
#'
mut2ms6 <- function(mutation_vec,sep = ">") {

  K = unique(nchar(stringr::str_extract(
    pattern = glue::glue("[:upper:]+(?={sep})"),
    string = mutation_vec)))

  stopifnot(length(K) == 1)

  k = (K - 1) / 2

  stopifnot((K - 1) %% 2 == 0) # needs to hold

  cp = k + 1
  ap = K + 2

  paste(
    stringr::str_sub(mutation_vec,cp,cp),
    stringr::str_sub(mutation_vec,ap,ap),
    sep = sep
  )
}




