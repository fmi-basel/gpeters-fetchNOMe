test_that("get_ctable_from_bams works for QuasR bam file", {
  
  quasr_bam <- system.file("extdata",
                           "QuasR_test.bam",
                           package = "fetchNOMe")
  test_genome <- system.file("extdata","random_genome_700bp.fa",package = "fetchNOMe")
  quasR_coocc_matrix <- get_ctable_from_bams(bamfiles = quasr_bam,
                                             samplenames = "quasr_pe",
                                             regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                              strand = "+",
                                                                              IRanges::IRanges(start = 200,end = 500)),
                                             genome = test_genome,
                                             alignerUsed = "QuasR",
                                             min_frag_data_len = 0L,
                                             min_frag_data_dens = 0.0,
                                             max_spacing = 300,
                                             remove_nonunique = F,
                                             clip_until_nbg = 0L,
                                             max_bisC_meth = 1)
  
  quasR_expect_coocc_mat <- readRDS(test_path("testdata","QuasR_expected_coocc_mat.rds")) 
  
  expect_true(all(quasR_coocc_matrix$CoocCountTable[[1]]==quasR_expect_coocc_mat))
  
})

test_that("get_ctable_from_bams works for Bismark bam file", {
  bism_bam <- system.file("extdata",
                          "Bismark_test.bam",
                          package = "fetchNOMe")
  test_genome <- system.file("extdata","random_genome_700bp.fa",package = "fetchNOMe")
  
  bismark_coocc_matrix <- get_ctable_from_bams(bamfiles = bism_bam,
                                               samplenames = "bismark_pe",
                                               regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                strand = "+",
                                                                                IRanges::IRanges(start = 200,end = 500)),
                                               genome = test_genome,
                                               alignerUsed = "Bismark",
                                               min_frag_data_len = 0L,
                                               min_frag_data_dens = 0.0,
                                               max_spacing = 300,
                                               remove_nonunique = F,
                                               clip_until_nbg = 0L,
                                               max_bisC_meth = 1)
  
  expect_coocc_mat <- readRDS(test_path("testdata","QuasR_expected_coocc_mat.rds"))
  
  expect_true(all(bismark_coocc_matrix$CoocCountTable[[1]] == expect_coocc_mat*2))
  
})