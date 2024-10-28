


#' Fetch data for each molecule from NOMe-seq BAM files for a set of regions
#'
#'
#'
#' @param bamfiles \code{character} paths to bam files
#' @param samplenames \code{character} names of samples. If \code{collapseBySample} is \code{TRUE} counts are aggregated across bam files with same sample name.
#' @param regions \linkS4class{GRanges} object for regions of interest.
#' @param genome \code{character} path to fasta file which was used as reference for alignment.
#' @param whichContext \code{character} which can be one of "GCH","WCG","bisC", "otherC" or "allC" and sets which data to retrieve.
#' \describe{
#' \item{"GCH}{protection data;}
#' \item{"WCG"}{endogenous CpG methylation data;}
#' \item{"bisC}{methylation states of non-GCH and non-WCG cytosins used for estimating efficiency of bisulfite conversion;}
#' \item{"otherC"}{methylation states of all other cytosines which do not fall into above categories.}
#' \item{"allC"}{methylation states of all cytosines combined.}
#' }
#' @param alignerUsed \code{character} that defines which aligner was used to generate BAM files. Currently supported aligners are QuasR, Bismark and BISCUIT.
#' @param min_frag_data_len \code{integer} ignore fragments that have genomic lengths from most-left to most-right data points less than \code{min_frag_data_len}.
#' @param min_frag_data_dens \code{numeric} ignore fragments that have density of data-containing positions lower than \code{min_frag_data_dens}.
#' @param collapseBySample \code{logical} indicating whether to collapse counts for bam files with same \code{samplenames}. If \code{FALSE} prefix \code{s} followed by index is added to \code{samplenames}.
#' @param remove_nonunique \code{logical} if \code{TRUE} only unique fragments are analyzed. Uniqueness is defined by fragments start, end and methylation states of all cytosines.
#' @param clip_until_nbg \code{integer} controlling clipping partial footprints at both ends of fragments. Namely, protected states are omitted from each end until \code{clip_until_nbg} unprotected states are met.
#' Set \code{clip_until_nbg=0L} to skip clipping pre-processing.
#' @param noclip_protect_frac_above \code{numeric} controlling fragment clipping, i.e. fragments that have fraction of protected states higher than \code{max_protect_frac} escape clipping. This parameter allows to avoid removal of fully protected fragments during clipping process.
#' @param max_bisC_meth \code{numeric} controlling removal of fragments with failed bisulfite conversion, i.e. fragments with fraction of methylated non-GCH, and non-WCG cytosines higher then \code{max_bisC_meth} are removed from the analysis to get rid of fragments for which bisulfite conversion failed.
#' @param min_bisC_size \code{integer} setting minimum number non-GCH, and non-WCG cytosines for filtering according to bisulfite conversion efficiency,
#' i.e. only fragments with number of non-GCH, and non-WCG cytosines higher then \code{min_bisC_size} are subject to filtering based on bisulfite conversion efficiency.
#' @param mapqMin \code{integer} setting minimum MAPQ of reads to be included into analysis.
#' @param mapqMax \code{integer} setting maximum MAPQ of reads to be included into analysis.
#' @param max_read_size \code{integer}. maximum read length which is expected in the data. This parameter controls only extension of regions for extracting reference sequence.
#' @param ncores number of cores to use.
#' @param data_as_Rle \code{logical} if \code{TRUE} data vectors are \code{Rle} compressed.
#' @return \code{tibble} where each row corresponds to sample - molecule combination.
#' Please note that if \code{regions} contain overlapping ranges it is possible that same molecule to appear several times.
#' 
#' Columns of the returned \code{tibble} represent:
#' \describe{
#'   \item{SampleName}{Name of a sample as defined by \code{samplenames}.}
#'   \item{seqnames, start, end}{seqnames (or chromosomes) and coordinates of molecules. Please note, that start and end are positions of the first and last data points and not positions of alignments.}
#'   \item{R1strand}{strands of alignments for R1 reads}
#'   \item{qname}{names of a molecules}
#'   \item{nDataPoints}{numbers of informative data points within molecules}
#'   \item{startDataPoint}{most-left data points}
#'   \item{endDataPoint}{most-right data points}
#'   \item{fragData_list}{vectors with molecule data
#'   }
#' }
#' @export
#'

get_molecule_data_list_from_bams <- function(bamfiles,
                                      samplenames,
                                      regions,
                                      genome,
                                      whichContext = c("GCH","WCG","bisC","otherC", "allC"),
                                      alignerUsed = c("QuasR","Bismark","BISCUIT"),
                                      min_frag_data_len = 50L,
                                      min_frag_data_dens = 0.05,
                                      collapseBySample = TRUE,
                                      remove_nonunique = TRUE,
                                      clip_until_nbg  = 1L,
                                      noclip_protect_frac_above = 0.90,
                                      max_bisC_meth = 0.1,
                                      min_bisC_size = 10,
                                      mapqMin = 0,
                                      mapqMax = 255,
                                      max_read_size = 1000L,
                                      ncores = 1L,
                                      data_as_Rle = FALSE){
  if(!inherits(bamfiles,"character"))
    stop("bamfiles must be vector of character")
  if(!inherits(samplenames,"character"))
    stop("samplenames must be vector of character")
  if(any(!file.exists(bamfiles))){
    nonex <- which(!file.exists(bamfiles))
    nonex <- paste0(bamfiles[nonex],collapse = "\n")
    stop(paste0("Could not find following BAM files:\n",nonex))
  }
  
  if(!inherits(regions,"GRanges"))
    stop("regions must be Granges")
  
  
  if(!inherits(genome,"character") & length(genome) !=1)
    stop("genome must be character and contain only one element")
  
  if(!file.exists(genome))
    stop("Could not find genome:", genome)
  
  whichContext <- match.arg(whichContext)
  
  alignerUsed <- match.arg(alignerUsed)
  
  message(paste0("Fetching data for every molecule from BAM files generated by ",alignerUsed,"."))
  
  
  ## remove metadata from regions
  GenomicRanges::mcols(regions) <- NULL
  
  ## load refsequences for regions
  
  ## extend regions from left and right
  fa_index <- Rsamtools::scanFaIndex(genome)
  GenomeInfoDb::seqlevels(regions) <- GenomeInfoDb::seqlevels(fa_index)
  GenomeInfoDb::seqinfo(regions) <- suppressMessages(GenomeInfoDb::Seqinfo(seqnames = as.character(GenomeInfoDb::seqnames(fa_index)),
                                                                           seqlengths = GenomicRanges::width(fa_index)))
  refseq_gr <- suppressWarnings(GenomicRanges::resize(regions,width = GenomicRanges::width(regions) + 2*max_read_size,fix="center"))
  
  ## trim to remove parts out of chromosomes
  refseq_gr <- GenomicRanges::trim(refseq_gr)
  refseqs <- Rsamtools::scanFa(genome,refseq_gr)
  
  # 2. split vector by sample name. If collapseBySample is FALSE add to samplenames some suffixes
  if(!collapseBySample){
    samplenames <- paste0("s",seq_along(samplenames),"_",samplenames)
  }
  
  bamfiles <- split(bamfiles, f=samplenames)
  
  
  molec_data_list <- do.call(rbind,lapply(names(bamfiles),
                                        function(bamgroupname){
                                          ctbl <- do.call(rbind,parallel::mclapply(seq_along(regions),
                                                                                   function(regi){
                                                                                     
                                                                                     regdf <- GenomicRanges::as.data.frame(regions[regi])
                                                                                     regdf$seqnames <- as.character(regdf[,"seqnames"])
                                                                                     seqgr <- as.data.frame(refseq_gr[regi])
                                                                                     
                                                                                     outlist <- fetch_molec_data_list_from_bams_cpp(whichContext = whichContext,
                                                                                                                                    infiles = bamfiles[[bamgroupname]],
                                                                                                                                    regionChr = regdf[1,"seqnames"],
                                                                                                                                    regionStart = regdf[1,"start"] - 1, # convert to 0-based
                                                                                                                                    regionEnd = regdf[1,"end"],
                                                                                                                                    remove_nonunique = remove_nonunique,
                                                                                                                                    seqstring = as.character(refseqs[regi]),
                                                                                                                                    seqStart = seqgr[1,"start"] - 1,
                                                                                                                                    seqEnd = seqgr[1,"end"],
                                                                                                                                    clip_until_nbg  = clip_until_nbg,
                                                                                                                                    max_protect_frac = noclip_protect_frac_above,
                                                                                                                                    max_bisC_meth = max_bisC_meth,
                                                                                                                                    min_bisC_size = min_bisC_size,
                                                                                                                                    mapqMin = mapqMin,
                                                                                                                                    mapqMax = mapqMax,
                                                                                                                                    alignerUsed = alignerUsed,
                                                                                                                                    min_frag_data_len = min_frag_data_len,
                                                                                                                                    min_frag_data_dens = min_frag_data_dens,
                                                                                                                                    data_as_rle = data_as_Rle)
                                                                                     
                                                                                     
                                                                                     outtbl <- tibble::tibble("SampleName" = bamgroupname,
                                                                                                              "seqnames" = regdf[1,"seqnames"],
                                                                                                              "start" = outlist[["start"]],
                                                                                                              "end" = outlist[["end"]],
                                                                                                              "R1strand" = outlist[["R1strand"]],
                                                                                                              "fragConfig" = outlist[["fragConfigs"]],
                                                                                                              "qname" = outlist[["qnames"]],
                                                                                                              "nDataPoints" = outlist[["nDataPoints"]],
                                                                                                              "startDataPoint" = outlist[["startDataPoint"]],
                                                                                                              "endDataPoint" = outlist[["endDataPoint"]],
                                                                                                              "fragData_list" = outlist[["fragData_list"]])
                                                                                     colnames(outtbl)[which(colnames(outtbl) == "fragData_list")] <- paste0(whichContext,"_fragData_list")
                                                                                     
                                                                                     return(outtbl)
                                                                                   },mc.cores = ncores))
                                          
                                          return(ctbl)
                                        }))
}