test_that("get_data_matrix_from_bams works for QuasR bam file", {
  #library(GenomicRanges)
  #browser()
	quasR_methmat <- get_data_matrix_from_bams(bamfiles = "QuasR_test.bam",
																						 samplenames = "quasr_pe",
																						 regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
																						 																 strand = "+",
																						 																 IRanges::IRanges(start = 200,end = 500)),
																						 genome = "random_genome_700bp.fa",
																						 whichContext ="allC",
																						 remove_nonunique = F,
																						 clip_until_nbg = 0L,
																						 max_bisC_meth = 1)
	## substitue NA to 2
	roword <- paste0("test_case_",c("OT","CTOT","OB","CTOB"))
	quasR_methmat$allC_DataMatrix[[1]][is.na(quasR_methmat$allC_DataMatrix[[1]])] <- 2
	row.names(quasR_methmat$allC_DataMatrix[[1]]) <- gsub("f1:","",row.names(quasR_methmat$allC_DataMatrix[[1]]))
	quasR_methmat$allC_DataMatrix[[1]] <- quasR_methmat$allC_DataMatrix[[1]][roword,]
	
	## load expected matrix
	quasr_expected_methmat <- readRDS("QuasR_expected_methylation_data_matrix.rds")[roword,]
	
	## test that identical
	expect_true(all(quasR_methmat$allC_DataMatrix[[1]] == quasr_expected_methmat))
	
	
	# #### check if filtering by data vector length works ####
	quasR_methmat <- get_data_matrix_from_bams(bamfiles = "QuasR_test.bam",
	                                           samplenames = "quasr_pe",
	                                           regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
	                                                                            strand = "+",
	                                                                            IRanges::IRanges(start = 200,end = 500)),
	                                           genome = "random_genome_700bp.fa",
	                                           whichContext ="allC",
	                                           remove_nonunique = F,
	                                           clip_until_nbg = 0L,
	                                           max_bisC_meth = 1,
	                                           min_frag_data_len = 500L)
	
	# ## test whether 0 reads are fetched
	expect_true(nrow(quasR_methmat$allC_DataMatrix[[1]]) == 0)
	quasR_methmat <- get_data_matrix_from_bams(bamfiles = "QuasR_test.bam",
	                                           samplenames = "quasr_pe",
	                                           regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
	                                                                            strand = "+",
	                                                                            IRanges::IRanges(start = 200,end = 500)),
	                                           genome = "random_genome_700bp.fa",
	                                           whichContext ="allC",
	                                           remove_nonunique = F,
	                                           clip_until_nbg = 0L,
	                                           max_bisC_meth = 1,
	                                           min_frag_data_dens  = 1.0)
	expect_true(nrow(quasR_methmat$allC_DataMatrix[[1]]) == 0)
	
	
})


test_that("get_data_matrix_from_bams works for Bismark bam file with indels and soft clipping", {
#   library(GenomicRanges)
# 	## coordinates of the test fragment
# 	frag_coord_OT_gr <- GRanges(seqnames = "random_genome_700bp",
# 															strand = "+",
# 															IRanges(start = 200,end = 500))
	
	bismark_methmat <- get_data_matrix_from_bams(bamfiles = "Bismark_test.bam",
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
	roword <- paste0("test_case_",c("OT","OT_I_D","CTOT","CTOT_I_D","OB","OB_I_D", "CTOB","CTOB_I_D"))
	bismark_methmat$allC_DataMatrix[[1]][is.na(bismark_methmat$allC_DataMatrix[[1]])] <- 2
	row.names(bismark_methmat$allC_DataMatrix[[1]]) <- gsub("f1:","",row.names(bismark_methmat$allC_DataMatrix[[1]]))
	bismark_methmat$allC_DataMatrix[[1]] <- bismark_methmat$allC_DataMatrix[[1]][roword,]
	
	## load expected matrix
	bismark_expected_methmat <- readRDS("Bismark_expected_methylation_data_matrix.rds")[roword,]
	
	
	## test that identical
	expect_true(all(bismark_methmat$allC_DataMatrix[[1]] == bismark_expected_methmat))
	
	# #### test if filtering by data vector length works ####
	bismark_methmat <- get_data_matrix_from_bams(bamfiles = "Bismark_test.bam",
	                                             samplenames = "bismark_pe",
	                                             regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
	                                                                              strand = "+",
	                                                                              IRanges::IRanges(start = 200,end = 500)),
	                                             genome = "random_genome_700bp.fa",
	                                             whichContext ="allC",
	                                             remove_nonunique = F,
	                                             clip_until_nbg = 0L,
	                                             min_frag_data_len = 500L)
	expect_true(nrow(bismark_methmat$allC_DataMatrix[[1]]) == 0)
	bismark_methmat <- get_data_matrix_from_bams(bamfiles = "Bismark_test.bam",
	                                             samplenames = "bismark_pe",
	                                             regions = GenomicRanges::GRanges(seqnames = "random_genome_700bp",
	                                                                              strand = "+",
	                                                                              IRanges::IRanges(start = 200,end = 500)),
	                                             genome = "random_genome_700bp.fa",
	                                             whichContext ="allC",
	                                             remove_nonunique = F,
	                                             clip_until_nbg = 0L,
	                                             min_frag_data_dens = 1.0
	                                             )
	expect_true(nrow(bismark_methmat$allC_DataMatrix[[1]]) == 0)
	
	
	
})