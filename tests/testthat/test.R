
context("Genome object errors")
test_that("muts context helpers", {
  expect_error(get_MS_VR(genome = "not a genome"))
  expect_error(compute_MSR(vr = "djf"))
  expect_error(simplify_muts(muts = character()))
  expect_error(simplify_ctx(muts = character()))
})


context("Simplify ctx")
test_that("Simplify", {
  ctx = c("CTC","GGT","GAT")
  res = simplify_ctx(ctx = ctx)
  expect_equal(c("GAG","ACC","GAT"),res)
})


context("Mutation context helpers")
test_that("muts context helpers", {
  # moved
  expect_true(all(generate_mut_types(1) == pos_ms96))

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


context("t3")
test_that("MS", {
  expect_equal(make_set("ART>T"),c("AAT>T","ACT>A"))
  expect_equal(make_set("YTCAA>T"),c("CTCAA>T","TTCAA>T"))
  expect_equal(make_set("AMT>T"),c("AAT>T","ACT>T"))
})


context("t5")
test_that("cosine sim", {
  expect_equal(cos_sim_vect(c("a"=1,"b"=10,"c"=100),c("a"=1,"c"=100,"b"=10)),1)
})

context("get_MS_VR")
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

context("MSM")
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

  # genome = genome_selector(alias="Hsapiens.UCSC.hg19")
  # m1 = expect_warning(compute_MSM(vr))
  # m2 = expect_warning(compute_MSM_fast(vr))
  # expect_equal(m1,m2)

  # m3 = mutationContext(vr = vr,
  #                      ref = genome_selector(alias="Hsapiens.UCSC.hg19"),
  #                      simplify_set = c("C","A"))
  # m3 = glue::glue("{stringr::str_sub(m3$context,1,1)}{stringr::str_sub(m3$alteration,1,1)}{stringr::str_sub(m3$context,3,3)}>{stringr::str_sub(m3$alteration,2,2)}")
  #m4 = expect_warning(get_MS_VR(vr))
  # expect_equal(as.character(m3),m4)

  #expect_equal(m4,c("TAA>T","CCA>A"))

  #expect_equal(sum(correct_MSM_P(m2,c(10,100))), 110)
})


context("cosine")
test_that("cosine sim", {
  expect_equal(jaccard(c("a","a","c"),c("a")),.5)


  res = .lenunique(c("A","C","A"))
  expect_equal(res,2)
})


context("make sets")
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


# test utils

context("microhomology")
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

#
# # tests
#
# detect_microhomology(seqvec)

# 0 3
# 2 0
#


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

