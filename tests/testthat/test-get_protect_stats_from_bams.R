test_that("get_protect_stats_from_bams works for QuasR bam file", {
  quasr_bam <- system.file("extdata",
                           "QuasR_test.bam",
                           package = "fetchNOMe")
  test_genome <- system.file("extdata","random_genome_700bp.fa",package = "fetchNOMe")
  
  quasr_protect_stats <- get_protect_stats_from_bams(bamfiles = quasr_bam,
                                                     samplenames = "quasr_pe",
                                                     regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                      strand = "+",
                                                                                      IRanges::IRanges(start = 200,end = 500)),
                                                     genome = test_genome,
                                                     alignerUsed = "QuasR",
                                                     min_frag_data_len = 0L,
                                                     min_frag_data_dens = 0.0,
                                                     remove_nonunique = F,
                                                     clip_until_nbg = 0L,
                                                     max_bisC_meth = 1)
  ## order as currently the order is random due to unordered_map
  quasr_protect_stats$ProtectStats[[1]] <- quasr_protect_stats$ProtectStats[[1]][order(quasr_protect_stats$ProtectStats[[1]][,1],
                                                                                       quasr_protect_stats$ProtectStats[[1]][,2]),]
  
  quasr_expect_protstats <- readRDS(test_path("testdata","QuasR_expected_protectstats.rds"))
  
  expect_true(all(quasr_protect_stats$ProtectStats[[1]] == quasr_expect_protstats))
})


test_that("get_protect_stats_from_bams works for Bismark bam file", {
  
  bism_bam <- system.file("extdata",
                          "Bismark_test.bam",
                          package = "fetchNOMe")
  test_genome <- system.file("extdata","random_genome_700bp.fa",package = "fetchNOMe")
  
  bismark_protect_stats <- get_protect_stats_from_bams(bamfiles = bism_bam,
                                                       samplenames = "bismark_pe",
                                                       regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                        strand = "+",
                                                                                        IRanges::IRanges(start = 200,end = 500)),
                                                       genome = test_genome,
                                                       alignerUsed = "Bismark",
                                                       min_frag_data_len = 0L,
                                                       min_frag_data_dens = 0.0,
                                                       remove_nonunique = F,
                                                       clip_until_nbg = 0L,
                                                       max_bisC_meth = 1)
  ## order as currently the order is random due to unordered_map
  bismark_protect_stats$ProtectStats[[1]] <- bismark_protect_stats$ProtectStats[[1]][order(bismark_protect_stats$ProtectStats[[1]][,1],
                                                                                           bismark_protect_stats$ProtectStats[[1]][,2]),]
  
  expect_protstats <- readRDS(test_path("testdata","QuasR_expected_protectstats.rds"))
  expect_protstats[,3] <- 2*expect_protstats[,3]
  expect_true(all(bismark_protect_stats$ProtectStats[[1]] == expect_protstats))
})