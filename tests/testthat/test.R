#$$$$$$$$\ $$$$$$$$\  $$$$$$\ $$$$$$$$\  $$$$$$\
#\__$$  __|$$  _____|$$  __$$\\__$$  __|$$  __$$\
#   $$ |   $$ |      $$ /  \__|  $$ |   $$ /  \__|
#   $$ |   $$$$$\    \$$$$$$\    $$ |   \$$$$$$\
#   $$ |   $$  __|    \____$$\   $$ |    \____$$\
#   $$ |   $$ |      $$\   $$ |  $$ |   $$\   $$ |
#   $$ |   $$$$$$$$\ \$$$$$$  |  $$ |   \$$$$$$  |
#   \__|   \________| \______/   \__|    \______/


## README:
##
## I am not entirely sure about the meaning of context and the test that
## text argument. For now, I will organize this with a section per R/script.R
## and each "group" will have the same context. The name of the test that
## will be diferent.
##



# GENOME ------------------------------------------------------------------

context("genome")
test_that("muts context helpers", {
  expect_error(get_MS_VR(genome = "not a genome"))
  expect_error(compute_MSR(vr = "djf"))
  expect_error(simplify_muts(muts = character()))
  expect_error(simplify_ctx(muts = character()))
})


context("genome")
test_that("genome selector", {

  test = genome_selector()
  expect_equal(length(seqnames(test)), 93)

})

context("genome")
test_that("chunk regions", {


  expect_error(get_region_chunks(gr = "wrong"))
  test = genome_selector()
  res = get_region_chunks(test)

  expect_true(all(width(res) <= 10000))
  expect_equal(sum(width(trim(res))), sum(seqlengths(test)))
})



# MUTS --------------------------------------------------------------------


context("enrichment")
test_that("general enrichment", {
  library(VariantAnnotation)
  genome = genome_selector()
  base_pos = 6e4
  vr_target = VRanges(seqnames="chr1",
                      ranges=IRanges(c(base_pos + 6, base_pos + 16),
                                     width = 1),
                      ref = "C",alt = "A")
  gr_target <- GRanges(seqnames="chr1",
                       ranges=IRanges(base_pos + 3,base_pos +  10))
  gr_mask <- GRanges(seqnames=c("chr1", "chr1"),
                     ranges=IRanges(c(base_pos + 4,base_pos + 15),
                                    width = 10))

  res = mutation_enrichment_general(vr = vr_target,
                              gr = gr_target,
                              genome = genome,
                              genome_mask = gr_mask)

  expect_s3_class(object = res,class = "data.frame")
  expect_equal(round(as.numeric(res$estimate),1), 1.80)
})


context("muts")
test_that("str functions", {
  input = "ACGGT"
  output = "ACCGT"

  expect_equal(str_reverse_complement(input),output)

  input = "ACGGT>T"
  output = "ACCGT>A"

  expect_equal(ms_reverse_complement(input),output)
})


context("muts")
test_that("Simplify", {
  ctx = c("CTC","GGT","GAT")
  res = simplify_ctx(ctx = ctx)
  expect_equal(c("GAG","ACC","GAT"),res)
})


context("muts")
test_that("muts context helpers", {
  # this is to prove that the mutation types are indeed alphabetical
  all_mut_types = generate_mut_types(1,simplify_set = c("C","T"))
  expect_true(all(all_mut_types == order_ms96_cosmicSignatures))

  # general simplifly
  muts = c("AGT>T","ACT>T")
  res = simplify_muts(muts,simplify_set = c("T","G"))
  expect_equal(res,c("AGT>T","AGT>A"))

  # simplyfy muts see issue #25 unique mutation quantifier
  ms = "TCA>T"
  res = simplify_muts(muts = ms,simplify_set = "G",sep = ">")
  expect_equal(res,"TGA>A")
  res = simplify_muts(muts = ms,simplify_set = "C",sep = ">")
  expect_equal(res,"TCA>T")
  # issue #25
  res = expect_warning(simplify_muts(muts = ms,simplify_set = "T",sep = ">"))
  expect_true(is.na(res))

  # indentify muts
  muts = c("TCA>T","TCCAA>T")
  expect_error(identify_mut_aestetics(muts))
  muts = c("TSA>A")
  expect_error(identify_mut_aestetics(muts))

  muts = c("TCA>T","TAC>T")

  # see issue #23
  muts = c("NNN>T")
  expect_error(count_MS(muts))

  # see issue #23
  muts = c("AGT>T","ACT>T")
  expect_error(count_MS(ms = muts),regexp = "simplif")

  # see issue #24 too
  CM = expect_warning(count_MS(simplify_muts(muts,simplify_set = c("T","G"))))
  expect_equal(sum(CM),2)

  # see issue #22 too
  CM = expect_warning(count_MS(simplify_muts(muts,simplify_set = c("C","T"))))
  expect_equal(sum(CM),2)

  # test if the pasting works.
  muts2 = gsub(pattern = ">",replacement = ":",x = muts)
  CM = expect_warning(count_MS(simplify_muts(muts2,simplify_set = c("T","G"),sep=":")))
  expect_equal(unique(stringr::str_sub(names(CM),4,4)) , ":")

  ms = c("CCG>A","CCT>A","GCA>A","GCC>A","CCT>A")
  res = expect_warning(count_MS(ms))

  expect_equal(as.numeric(res["GCC>A"]),1)
  expect_equal(as.numeric( res["CCT>A"]),2)

  expect_equal(mut2ms6("TCA>T"),"C>T")
  expect_equal(mut2ms6("TAA>T"),"A>T")
  expect_equal(mut2ms6("TTCAA>T"),"C>T")
  expect_equal(mut2ms6("C>T"),"C>T")
  expect_equal(mut2ms6("TTCAA~T",sep = "~"),"C~T")
  expect_error(mut2ms6("TTCA~T",sep = "~")) # K needs to be multiple than 2
})


context("muts")
test_that("MS", {
  expect_equal(make_set("ART>T"),c("AAT>T","ACT>A"))
  expect_equal(make_set("YTCAA>T"),c("CTCAA>T","TTCAA>T"))
  expect_equal(make_set("AMT>T"),c("AAT>T","ACT>T"))
})


context("muts")
test_that("MS",{
  library(VariantAnnotation)
  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  ## GGA
  ## GAG
  vr = VRanges(seqnames = "chr1",
               ranges = IRanges::IRanges(start = c(186716278,186716284),
                                         width = c(1)),
               ref = c("A","G"),
               alt = c("T","T"))
  res1 = expect_warning(get_MS_VR(x = vr,genome = genome))

  strand(vr) = "+"
  res2 = get_MS_VR(x = vr,genome = genome)

  strand(vr) = c("-","+")
  res3 = expect_warning(get_MS_VR(x = vr,genome = genome,keep_strand = TRUE))

  expect_equal(res1,res2)
  expect_equal(res1,res3)
  expect_equal(res1,c("TCC>A","GAG>T"))

  strand(vr) = c("+","-")
  res4 = expect_warning(get_MS_VR(x = vr,genome = genome,keep_strand = TRUE))
  expect_equal(res4,c("GGA>T","CTC>A"))

})

context("muts")
test_that("MSM",{
  library(VariantAnnotation)
  vr = VariantAnnotation::VRanges(
    seqnames = c("chr1","chr2"),
    ranges = IRanges(start = c(2000000,200000),
                     end = c(2000000,200000)),
    ref = c("A","G"),
    alt = c("T","T"),
    sampleNames = c("test1","test2")
  )

  vr_test = vr
  seqlevelsStyle(vr_test) = "NCBI"
  # the seqlevels will have different SQlevels
  expect_error(get_MS_VR(x = vr_test,
                         genome = genome_selector(alias="Hsapiens.UCSC.hg19")))

  # contexts need to be simplified
  expect_error(count_MS(c("TCA>T","TCT>T","TGA>T")))
  # this is for non-standard mutations
  expect_error(count_MS(c("TCA>T","TCT>T","TCN>T")))

  genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  m2 = expect_warning(compute_MSM_fast(vr,genome = genome))

  m3 = SomaticSignatures::mutationContext(vr = vr,
                                          ref = genome)
  m3 = glue::glue("{stringr::str_sub(m3$context,1,1)}{stringr::str_sub(m3$alteration,1,1)}{stringr::str_sub(m3$context,3,3)}>{stringr::str_sub(m3$alteration,2,2)}")

  m4 = expect_warning(get_MS_VR(vr,genome = genome,simplify_set = c("C","T")))
  expect_equal(as.character(m3),m4)

  expect_equal(m4,c("TTA>A","CCA>A"))

  expect_equal(sum(correct_MSM_P(m2,c(10,100))), 110)
})


context("muts")
test_that("sets strict",{
  set_test = expect_warning(make_set(x = "NNN>N",simplify = T))
  expect_equal(length(unique(set_test)),96)

  imposible_set = make_set(x = "NNN>N",strict = FALSE)
  expect_equal(length(imposible_set),4^4)

  # this test that the function works with multiple sets.
  set_test2 = make_set(make_set(c("TCW>T","CCW>T")))
  expect_equal(set_test2,c("TCA>T","TCT>T","CCA>T","CCT>T"))

  set_test2 = make_set(c("TCW>T","CCW>T"),
                       simplify_set = c("T","G"))
  expect_equal(set_test2,c("TGA>A","AGA>A","TGG>A","AGG>A"))

})



context("muts")
test_that("extract params", {
  expect_error(helperMut::identify_mut_aestetics(c("TCA>T",
                                                   "CCA>T",
                                                   "TTT>A",
                                                   "TNA>C")))

  val = helperMut::identify_mut_aestetics(c("TCA>T",
                                            "CCA>T",
                                            "TTT>A",
                                            "TNA>C"),force = T)

  expect_true("N" %in% val$mset_ref)

  val = helperMut::identify_mut_aestetics(c("TCA>T",
                                            "CCA>T",
                                            "TTT>A"))

  expect_true(all(c("C","T") %in% val$mset_ref))

  expect_equal(val$k , 1)
  expect_equal(val$K , 3)
  expect_equal(val$sep , ">")

  val = helperMut::identify_mut_aestetics(c("ATCAA~T",
                                            "CCCAC~T",
                                            "CTTTC~A"))

  expect_equal(val$k , 2)
  expect_equal(val$K , 5)
  expect_equal(val$sep , "~")

  val = helperMut::identify_mut_aestetics(c("C~T",
                                            "C~T",
                                            "T~A"))

  expect_equal(val$k , 0)
  expect_equal(val$K , 1)
  expect_equal(val$sep , "~")

})




# REGIONS -----------------------------------------------------------------

context("regions")
test_that("regions", {

  gr_reg_test = GenomicRanges::GRanges(
    seqnames = "chr1",ranges = IRanges::IRanges(start = c(1,20,40,50),
                                                end = c(10,35,47,55))
  )

  mask_gr_test = GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(1,40),end = c(35,100))
  )

  res = shufle_regions2(regions = gr_reg_test,mask = mask_gr_test)


  library(GenomicRanges)
  ovr = findOverlaps(query = res,subject = mask_gr_test)

  # this tests if the result is in the mask
  expect_equal(length(unique(queryHits(ovr))) , length(gr_reg_test))

  # this test if the result is shufled
  expect_false(all(start(res) == start(gr_reg_test)))
})

# PROFILES ----------------------------------------------------------------

context("profiles")
test_that("cosine sim", {
  expect_equal(cos_sim_vect(c("a"=1,"b"=10,"c"=100),c("a"=1,"c"=100,"b"=10)),1)
  expect_warning(cos_sim_vect(x = c(1,2,3),c(2,3,4)))
  expect_error(cos_sim_vect(x = "jshdfs",y = matrix()))
})


context("profiles")
test_that("download sets", {
  expect_error(download_signature_set(version = "v3",synapse_apiKey = NULL))

  test1 = download_signature_set(version = "v3")
  expect_equal(ncol(test1),67)
  expect_equal(nrow(test1),96)
  test2 = download_signature_set() # defaults to v2
  expect_equal(ncol(test2),30)
  expect_equal(nrow(test2),96)

})


# DATA --------------------------------------------------------------------

context("data")
test_that("data is loaded correctly", {
  expect_equal(length(order_ms96_supekCell2017),96)
  expect_equal(length(order_ms96_cosmicSignatures),96)
  expect_equal(length(tr_colors),6)
})

# UTILS -------------------------------------------------------------------

context("utils")
test_that("binom test", {

  resdf = binom_test(x = c(1,2,3),n = c(2,4,6))
  expect_true(all(resdf$estimate == 0.5))
  expect_true(all(resdf$p.value == 1))

  resdf = binom_test(x = c(1,2,3),n = c(2,4,6),p = rep(0.1,3))
  expect_true(all(resdf$p.value < 0.2))

  binom.test(x = 1,n = 2,p = .1)$p.value -> singleP
  resdf$p.value[1] -> multP
  expect_equal(singleP,multP)
})


context("utils")
test_that("random", {
  expect_equal(jaccard(c("a","a","c"),c("a")),.5)


  res = .lenunique(c("A","C","A"))
  expect_equal(res,2)
})



# INDELS ------------------------------------------------------------------

context("indels")
test_that("microhomology", {
  seqvec = c("CTACTCGTTCGA","TTTCCTACTATTTCCGGGTAGGGGGGGGTA")

  res = detect_microhomology(seqvec)
  expect_equal(res[1,1],0)
  expect_equal(res[1,2],3)
  expect_equal(res[2,1],2)
  expect_equal(res[2,2],0)

  seqvec_err = c("CTACTGTTCGA","TTTCCTACTATTTCCGGTAGGGGGGGGTA")
  expect_error(detect_microhomology(seqvec_err))
})



# UNCLASSIFIED ------------------------------------------------------------

