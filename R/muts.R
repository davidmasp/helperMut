### ===================================== ###
#         __  __ _    _ _______ _____
#        |  \/  | |  | |__   __/ ____|
#        | \  / | |  | |  | | | (___
#        | |\/| | |  | |  | |  \___ \
#        | |  | | |__| |  | |  ____) |
#        |_|  |_|\____/   |_| |_____/
#
### ===================================== ###


# function extracted and modified from SomaticSignatures for VRanges.
# Gehring JS, Fischer B, Lawrence M and Huber W (2015).
# Bioinformatics. doi: 10.1093/bioinformatics/btv408.
# see this
# https://github.com/juliangehring/SomaticSignatures/blob/master/R/mutation-context.R#L1-L52
mutationContext <- function(vr,
                            ref,
                            k = 3,
                            strand = FALSE,
                            unify = TRUE,
                            simplify_set = c("C","A"),
                            check = FALSE) {
  #browser()
  ## only SNV substitutions beyond this point
  if(!all(ref(vr) %in% DNA_BASES & alt(vr) %in% DNA_BASES))
    stop("Only SNV substitutions are currently supported.")
  if(k %% 2 != 1)
    stop("'k' must be odd.")
  mid = (k + 1)/2

  gr = granges(vr) ## drop mcols
  ranges = resize(gr, k, fix = "center")
  context = getSeq(ref, ranges)

  ref_base = DNAStringSet(ref(vr))
  alt_base = DNAStringSet(alt(vr))

  ## check against ref
  if(check) {
    ref0 = subseq(context, mid, mid) ## reference base from 'ref', for checking
    idx_invalid = ( ref0 != ref_base )
    if(any(idx_invalid))
      warning(sprintf("References do not match in %d cases",
                      sum(idx_invalid)))
  }

  ## convert to the plus strand, generally not needed
  if(strand) {
    s = strand(gr)
    if(any(s == "*"))
      stop("The strand must be explicit, in order to read the correct strand.")
    idx_minus = (s == "-")
    context[idx_minus] = reverseComplement(context[idx_minus])
    s[idx_minus] = "+"
    strand(gr) = s
  }

  ## convert to alterations starting with "C/T"
  if(unify) {
    idx_complement = !(as.character(ref_base) %in% simplify_set)
    context[idx_complement] = reverseComplement(context[idx_complement])
    ref_base[idx_complement] = reverseComplement(ref_base[idx_complement])
    alt_base[idx_complement] = reverseComplement(alt_base[idx_complement])
  }

  subseq(context, mid, mid) = "."
  alteration = xscat(ref_base, alt_base)

  vr$alteration = alteration
  vr$context = context

  return(vr)
}


#' GRanges extension
#'
#' @description
#' The function extend all the ranges in a \link[GenomicRanges]{GRanges} object by the number
#' of nucleotides defined in the upstream or downstream of each range.
#'
#' It also works for \link[VariantAnnotation]{VRanges} objects.
#'
#' extracted from \link{https://support.bioconductor.org/p/78652/}
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
  if (any(BiocGenerics::strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- BiocGenerics::strand(x) == "+" | BiocGenerics::strand(x) == "*"
  new_start <- BiocGenerics::start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- BiocGenerics::end(x) + ifelse(on_plus, downstream, upstream)
  IRanges::ranges(x) <- IRanges::IRanges(new_start, new_end)
  x = IRanges::trim(x) # this removes the negative Ranges
  return(x)
}


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


  if (!any(GenomeInfoDb::seqlevelsStyle(x) %in% GenomeInfoDb::seqlevelsStyle(genome))){
    stop("Seq levels do not match in regions and genome, check the reference")
  }

  tmp <- extend(x,upstream = k,downstream = k)
  seq <-  Biostrings::getSeq(genome,tmp)

  # from here can substitute for simplify_muts function
  alt_seq = VariantAnnotation::alt(x)
  #browser()
  muts = paste(seq,alt_seq,sep = sep)

  if (!keep_strand){
    muts = simplify_muts(muts,simplify_set = simplify_set)
  }

  return(muts)
}



str_reverse_complement <- function(str){
  as.character(
    Biostrings::reverseComplement(
      Biostrings::DNAString(str)
  ))
}


count_MS <- function(ms) {

  mutaes = identify_mut_aestetics(ms)

  k = mutaes[["k"]]
  simplify_set = mutaes[["mset_ref"]]
  sep = mutaes[["sep"]]

  if (length(simplify_set) > 2) {
    stop("Mutation type input is not simplified, stopping")
  } else if (length(simplify_set) == 2 && simplify_set[1] == str_reverse_complement(simplify_set[2])){
    stop("Mutation type input is not simplified, stopping")
  } else if (length(simplify_set) == 1){
    warning("A unique central position is presented, returning just those subtypes")
  }

  mut_types = generate_mut_types(k=k,
                                 sep=sep,
                                 simplify_set = simplify_set)

  if (!all(ms %in% mut_types)){
    stop("Some mutations are not standard. Mutations should be simplified")
  }

  row_counts = integer((4**(k*2)) * length(simplify_set)*3)
  names(row_counts) = mut_types

  idx = match(ms,table = mut_types)
  freqs = plyr::count(idx)

  row_counts[freqs$x] <- freqs$freq

  # this should be fine if this has passed the previous test
  stopifnot(sum(row_counts) == length(ms))

  return(row_counts)
}


compute_MSR <- function(vr,
                        k=1,
                        sep=">",
                        tp = NULL,
                        genome =  genome_selector()){

  # check class
  class_check =  is(object =vr,class2 ="VRanges")
  if (!class_check){
    stop("VRanges object needed.")
  }

  # check unisample
  stopifnot(length(unique(VariantAnnotation::sampleNames(vr)))==1)

  # get contexts
  ms = get_MS_VR(vr,k=k,simplify_set = c("A","C"), genome =  genome)
  row = count_MS(ms)

  # correct by total positives if needed
  if(!is.null(tp)){
    if (tp > length(vr)){
      row = correct_MSR_TP(x = row,tp=tp) # open issue
    }
  }

  return(row)
}


correct_MSR_TP <- function(x,tp){
  freq = x/sum(x)
  vals = freq*tp
  return(vals)
}


#' Simplify mutations
#'
#' @description
#' Given a set of mutations it removes the strand information and simplifies them.
#' It also translate a mutation set to an alternative refernce set (for instance, C T based to C A based.)
#'
#' @param muts A character vector containing mutations
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
#' simplify_ctx(mutset)
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

  central_pos = stringr::str_sub(seq,start = cp,end= cp)

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

  central_pos = stringr::str_sub(seq,start = cp,end= cp)

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
      dna_codes[[x]]
    })

  # yes. but is should be a list.
  sb = list(dna_codes[[b]])

  se = 1:nchar(e) %>%
    purrr::map(function(x){
      substr(x = e,start = x,stop = x)
    }) %>%
    purrr::map(function(x){
      dna_codes[[x]]
    })

  # this should be a vector because how expand.grid works.
  sto = dna_codes[[to]]

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

dna_codes <- list(
  "A" = c("A"),
  "C" = c("C"),
  "T" = c("T"),
  "G" = c("G"),
  "S" = c("C", "G"), # strong
  "W" = c("A", "T"), # weak
  "M" = c("A", "C"), # amino
  "K" = c("G", "T"), #keto
  "R" = c("A", "G"), #purine
  "Y" = c("C", "T"), # pyrimidine
  "B" = c("C", "T","G"), # not A
  "D" = c("A", "T","G"), # not C
  "H" = c("A", "T","C"), # not G
  "V" = c("A", "G","C"), # not T
  "N" = c("A", "C", "G", "T")
)

get_GR_from_gene_set <- function(GS,txdb){
  g = genes(txdb,
            filter = list("gene_id" = unlist(geneIds(gs) )))

  return(g)
}


get_gene_regions_occurences <- function(GR,k,ref_seq = NULL){
  library(GenomicRanges)
  library(magrittr)

  if (is.null(ref_seq)){
    assembly = unique(seqinfo(GR)@genome)
    if (is.na(assembly)){
      stop("A genome object is needed.")
    } else {
      stop("A genome object is needed.")
    }
  }

  s = Biostrings::getSeq(genome,GR)

  kmercounts = Biostrings::oligonucleotideFrequency(s,k)

  kmer_counts_vec = apply(kmercounts, 2,sum)

  nm = Biostrings::DNAStringSet(names(kmer_counts_vec))
  revComp = substr(nm,3,3) %in% c("T","G")
  nm[revComp] = Biostrings::reverseComplement( nm[revComp])
  nm_df = data.frame(ctx = nm,counts = kmer_counts_vec)
  occ = nm_df %>%
    dplyr::group_by(ctx) %>%
    dplyr::summarise(occurrences = sum(counts))

  return(occ)
}

test_general_enrichment <- function(GR,VR,k) {
  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  tmp <- extend(VR,upstream = k,downstream = k)
  seq <-  Biostrings::getSeq(genome,tmp)

  revComp = substr(seq, 3, 3) %in% c("T", "G")
  seq[revComp] = Biostrings::reverseComplement(seq[revComp])
  seq_df = data.frame(ctx = seq)
  occ = seq_df %>%
    dplyr::group_by(ctx) %>%
    dplyr::summarise(occurrences_obs = n())

  occ_exp = get_gene_regions_occurences(GR,(k*2)+1)
  occ_dat = dplyr::inner_join(occ,occ_exp,by="ctx") # change this for fill = 0
  fisher_df$ctx = occ_dat$ctx
  fisher_df$q.value = p.adjust(fisher_df$p.value,"fdr")
  fisher_df$b.value = p.adjust(fisher_df$p.value,"bonferroni")

  return(fisher_df)
}


fisher_muts <- function(obs,total){
  #browser()
  df = data.frame(obs = obs,
                  total = total)

  df %<>% dplyr::mutate(
    nomut = total - obs,
    rest_mut = sum(obs) - obs,
    rest_nomut = sum(nomut) - nomut
  ) %>% dplyr::select(-total)

  res_df = purrr::map_df(1:nrow(df),function(i){
    # critical step. depends inside the function tho
    tmp_mat = matrix(as.numeric(df[i,]),ncol =2)
    res_df = fisher.test(tmp_mat) %>% broom::tidy()

  })

  return(res_df)
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
#' @param ...
#'
#' @return A count matrix with samples in rows and mutation subtypes in columns
#' @export
#'
#' @examples
compute_MSM_fast = function(vr,
                            k=1,
                            sep = ">",
                            genome = genome_selector(),
                            simplify_set = c("C","A"),
                            ...) {

  # credit to somaticSignatures

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
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
correct_MSM_P <- function(msm,P,...){
  freqs = t(apply(msm, 1, function(row){
    row/sum(row)
  }))

  res = freqs * P

  return(res)

}

# this is a end function (should deprecate!!)
# dependent in clustMUT
compute_MSM <- function(vr,
                        k=1,
                        sep = ">",
                        tp = TRUE,
                        genome = genome_selector(),
                        simplify_set = c("C","A"),
                        ...) {

  ## internal funs
  ## extracted from https://stackoverflow.com/questions/18813526
  identicalValue <- function(x,y) if (identical(x,y)) x else FALSE

  stopifnot(requireNamespace("VariantAnnotation",quietly = T))
  gr_list = GenomicRanges::split(vr,VariantAnnotation::sampleNames(vr))

  # solving issue when a sampl does not have any clusters
  # could be a better way to handle this
  mask = lapply(gr_list, function(x){length(x)>0})
  gr_list = gr_list[unlist(mask)]

  mt = generate_mut_types(k = k,sep = sep,simplify_set = simplify_set)
  #browser()

  if (tp){
    samples = lapply(gr_list, function(x){
      compute_MSR(x,
                  k=k,
                  sep=sep,
                  tp = unique(x$tp),
                  genome=genome)
    })
    res_matrix = matrix(integer(length(mt)*length(samples)),ncol = length(mt))
  } else {
    samples = lapply(gr_list, function(x){
      compute_MSR(x,
                  k=k,
                  sep=sep,
                  genome=genome)
    })
    res_matrix = matrix(integer(length(mt)*length(samples)),ncol = length(mt))
  }

  rownames(res_matrix) <- names(samples)

  # see function in pca_sigs.R
  colnames_list = samples %>% purrr::map(colnames)
  colnames_chr = Reduce(identicalValue,colnames_list)
  colnames(res_matrix) <- mt
  for (i in names(samples)){
    tmp = samples[[i]]
    res_matrix[i,names(tmp)] = tmp
  }
  return(res_matrix)
}


identify_mut_aestetics = function(ms, force = FALSE){
  l = unique(nchar(ms))
  stopifnot(length(l)==1)

  k = (l - 3 )/ 2
  K = (k * 2) + 1

  mset_ref = stringr::str_sub(ms,k + 1,k + 1)

  cond = mset_ref %in% unique(unlist(dna_codes))
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




