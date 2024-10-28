test_that("get_molecule_data_list_from_bams works for QuasR bam file", {
  quasR_methlist <- get_molecule_data_list_from_bams(bamfiles = "QuasR_test.bam",
                                                     samplenames = "quasr_pe",
                                                     regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                      strand = "+",
                                                                                      IRanges::IRanges(start = 200,end = 500)),
                                                     genome = "random_genome_700bp.fa",
                                                     whichContext ="allC",
                                                     remove_nonunique = F,
                                                     clip_until_nbg = 0L,
                                                     max_bisC_meth = 1,
                                                     min_frag_data_len = 50L,
                                                     min_frag_data_dens = 0.05)
  
  roword <- paste0("f1:test_case_",c("OT","CTOT","OB","CTOB"))
  quasR_methlist <- quasR_methlist[match(roword,quasR_methlist$qname),]
  
  rdat <- lapply(quasR_methlist$allC_fragData_list,
                 function(x){
                   x[is.na(x)] <- 2
                   return(x)
                 })
  
  
  ## load expected data
  dat_reads <- readRDS("QuasR_expected_methylation_reads_list.rds")
  ## test if all data points equal expected values
  expect_true(all(do.call(c,lapply(1:length(rdat),function(i){
    rdat[[i]] == dat_reads[[i]]
  }))))
  ## test if all starts equal expected values
  expect_true(all(quasR_methlist$start == c(207,207,205,205)))
  expect_true(all(quasR_methlist$end == c(498,498,497,497)))
  
  
  ## test if returned Rle encoded data is as expected
  data_vec <- get_molecule_data_list_from_bams(bamfiles = "QuasR_test.bam",
                                               samplenames = "quasr_pe",
                                               regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                strand = "+",
                                                                                IRanges::IRanges(start = 200,end = 500)),
                                               genome = "random_genome_700bp.fa",
                                               whichContext ="GCH",
                                               remove_nonunique = F,
                                               clip_until_nbg = 0L,
                                               max_bisC_meth = 1,
                                               min_frag_data_len = 50L,
                                               min_frag_data_dens = 0.05,
                                               data_as_Rle = FALSE)
  data_rle <- get_molecule_data_list_from_bams(bamfiles = "QuasR_test.bam",
                                                     samplenames = "quasr_pe",
                                                     regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                      strand = "+",
                                                                                      IRanges::IRanges(start = 200,end = 500)),
                                                     genome = "random_genome_700bp.fa",
                                                     whichContext ="GCH",
                                                     remove_nonunique = F,
                                                     clip_until_nbg = 0L,
                                                     max_bisC_meth = 1,
                                                     min_frag_data_len = 50L,
                                                     min_frag_data_dens = 0.05,
                                                     data_as_Rle = TRUE)
  
  data_vec <- data_vec[match(data_rle$qname,data_vec$qname),]
  
  
  
  rdat_vec <- lapply(data_vec$GCH_fragData_list,
                     function(x){
                       x[is.na(x)] <- 2
                       return(x)
                     })
  rdat_rle <- lapply(data_rle$GCH_fragData_list,
                     function(x){
                       x <- as.vector(x)
                       x[is.na(x)] <- 2
                       return(x)
                     })
  expect_true(all(do.call(c,lapply(1:length(rdat_vec),function(i){
    rdat_vec[[i]] == rdat_rle[[i]]
  }))))
  
})


test_that("get_molecule_data_list_from_bams works for Bismark bam file with indels and soft clipping", {
  
  bismark_methlist <- get_molecule_data_list_from_bams(bamfiles = "Bismark_test.bam",
                                                       samplenames = "bismark_pe",
                                                       regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
                                                                                        strand = "+",
                                                                                        IRanges::IRanges(start = 200,end = 500)),
                                                       genome = "random_genome_700bp.fa",
                                                       whichContext ="allC",
                                                       remove_nonunique = F,
                                                       clip_until_nbg = 0L,
                                                       max_bisC_meth = 1)
  ## substitue NA to 2
  roword <- paste0("f1:test_case_",c("OT","OT_I_D","CTOT","CTOT_I_D","OB","OB_I_D", "CTOB","CTOB_I_D"))
  bismark_methlist <- bismark_methlist[match(roword,bismark_methlist$qname),]
  rdat <- lapply(bismark_methlist$allC_fragData_list,
                 function(x){
                   x[is.na(x)] <- 2
                   return(x)
                 })
  
  ## load expected data
  dat_reads <- readRDS("Bismark_expected_methylation_data_read_list.rds")

  ## test if all data points equal expected values
  expect_true(all(do.call(c,lapply(1:length(rdat),function(i){
    rdat[[i]] == dat_reads[[i]]
  }))))
  
  ## test if all starts and ends equal expected values
  expect_true(all(bismark_methlist$start == c(207,207,207,207,205,205,205,205)))
  expect_true(all(bismark_methlist$end == c(498,498,498,498,497,497,497,497)))
  
  
  
})


