


#' Fetch co-occurrence count table (ctable) of GCH protection states from BAM files
#' 
#' Given a set of bamfiles and regions of interest this function searches for C to T and G to A mismatches in sequencing reads for several tri-nucleotide contexts
#' and collects statistics of observing 0-0, 0-1, 1-0 and 1-1 at different distances.
#'
#' @param bamfiles \code{character} vector with paths to bam files
#' @param samplenames \code{character} names of samples. If \code{collapseBySample} is \code{TRUE} counts are aggregated across bam files with same sample name.
#' @param regions \linkS4class{GRanges} object for regions of interest.
#' @param genome \code{character} path to fasta file which was used as reference for alignment.
#' @param alignerUsed \code{character} that defines which aligner was used to generate BAM files. Currently supported aligners are QuasR, Bismark and BISCUIT.
#' @param SMFenzymeUsed \code{character} informs the function about what enzyme was used for marking accessible positions. 
#' If \code{M.CviPI} is chosen, then only cytosines in GCH are considered informative.
#' If \code{DddB} is chosen, then all cytosines are considered informative.
#' @param min_frag_data_len \code{integer} ignore fragments that have genomic lengths from most-left to most-right data points less than \code{min_frag_data_len}.
#' @param min_frag_data_dens \code{numeric} ignore fragments that have density of data-containing positions lower than \code{min_frag_data_dens}.
#' @param collapseBySample \code{logical} indicating whether to collapse counts for bam files with same \code{samplenames}. If \code{FALSE} prefix \code{s} followed by index is added to \code{samplenames}.
#' @param max_spacing maximum spacing length for which to collect co-occurrence counts.
#' @param remove_nonunique \code{logical} if \code{TRUE} only unique fragments are analyzed. Uniqueness is defined by fragments start, end and methylation states of all cytosines. 
#' @param clip_until_nbg \code{integer} controlling clipping partial footprints at both ends of fragments. Namely, protected states are erased from each end until \code{clip_until_nbg} unprotected states are met.
#' Set \code{clip_until_nbg=0L} to skip clipping preprocessing.
#' @param noclip_protect_frac_above \code{numeric} controlling fragment clipping, i.e. fragments which have fraction of protected states higher than \code{max_protect_frac} escape clipping. This parameter allows to avoid removal of fully protected fragments during clipping process.
#' @param max_bisC_meth \code{numeric} controlling removal of fragments with failed bisulfite conversion, i.e. fragments with fraction of methylated non-GCH, and non-WCG cytosines higher then \code{max_bisC_meth} are removed from the analysis to get rid of fragments for which bisulfite conversion failed.
#' @param min_bisC_size \code{integer} setting minimum number non-GCH, and non-WCG cytosines for filtering according to bisulfite conversion efficiency, 
#' i.e. only fragments with number of non-GCH, and non-WCG cytosines higher then \code{min_bisC_size} are subject to filtering based on bisulfite conversion efficiency.
#' @param mapqMin \code{integer} setting minimum MAPQ of reads to be included into analysis.
#' @param mapqMax \code{integer} setting maximum MAPQ of reads to be included into analysis.

#' @param max_read_size \code{integer}. maximum read length which is expected in the data. This parameter controls only extension of regions for extracting reference sequence.
#' @param ncores number of cores to use.
#'
#' @return \code{tibble} where each row corresponds to sample - region combination. 
#' Columns of the returned \code{tibble} represent:
#' \describe{
#'   \item{SampleName}{Name of a sample as defined by \code{samplenames}}.
#'   \item{seqnames, start, end}{seqnames (or chromosomes) and coordinates of corresponding regions.}
#'   \item{names}{the same as \code{names(regions)} if it is not \code{NULL} or index of a corresponding region in \code{regions}.}
#'   \item{nFragsFetched}{number of fetched fragments before filtering.}
#'   \item{nFragsNonUnique}{number of non-unique fragments as defined by coordinates and methylation profiles.}
#'   \item{nFragsBisFailed}{number of fragments for which bisulfite conversion may have failed, namely non-GCH and non-WCG methylation level is higher than \code{max_bisC_meth}.}
#'   \item{nFragsAnalyzed}{number of fragments analysed for co-occurrence table after filtering non-unique and bisulfite conversion failed fragments.}
#'   \item{ClippedUntilNbg}{parameter \code{clip_until_nbg} used for clipping partial footprints from both ends.}
#'   \item{NotClippedProtectAbove}{parameter \code{noclip_protect_frac_above} used for clipping partial footprints from both ends.}
#'   \item{nFragsClipped}{number of fragments clipped from both ends to get rid of paetial footprints.}
#'   \item{CoocCountTable}{\code{matrix} of dimension \code{max_spacing}x4 for number of observed combination of protection states at different window sizes (distances). 
#'   Row names correspond to window size and columns contain number of observed 0-0, 0-1, 1-0 and 1-1. Note that first row contains counts for window of size 1 (or distance 0) and equal to total number of 0s and 1s in the data.}
#' }
#' 
#' @export
#'

get_ctable_from_bams <- function(bamfiles,
                                 samplenames,
                                 regions,
                                 genome,
                                 alignerUsed = c("QuasR","Bismark","BISCUIT"),
                                 SMFenzymeUsed = c("M.CviPI","DddB"),
                                 min_frag_data_len = 50L,
                                 min_frag_data_dens = 0.05,
                                 collapseBySample = TRUE,
                                 max_spacing = 300,
                                 remove_nonunique = TRUE,
                                 clip_until_nbg  = 1L,
                                 noclip_protect_frac_above = 0.90,
                                 max_bisC_meth = 0.1,
                                 min_bisC_size = 10,
                                 mapqMin = 0,
                                 mapqMax = 255,
                                 max_read_size = 1000L,
                                 ncores = 1L){
  
  # 1. check arguments
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
  
  alignerUsed <- match.arg(alignerUsed)
  SMFenzymeUsed <- match.arg(SMFenzymeUsed)
  
  message(paste0("Collecting co-occurrence frequencies from BAM files generated by ",alignerUsed,". Assumed enzyme ",SMFenzymeUsed,"."))
  
  
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
  # 3. for each set of bam files call C++ function
  ctable_list <- do.call(rbind,lapply(names(bamfiles),
                                      function(bamgroupname){
                                        ctbl <- do.call(rbind,parallel::mclapply(seq_along(regions),
                                                                                 function(regi){
                                                                                   
                                                                                   regdf <- GenomicRanges::as.data.frame(regions[regi])
                                                                                   regdf$seqnames <- as.character(regdf[,"seqnames"])
                                                                                   seqgr <- as.data.frame(refseq_gr[regi])
                                                                                   outlist <- fetch_cooc_ctable_from_bams_cpp(infiles = bamfiles[[bamgroupname]],
                                                                                                                              regionChr = regdf[1,"seqnames"],
                                                                                                                              regionStart = regdf[1,"start"] - 1, # convert to 0-based
                                                                                                                              regionEnd = regdf[1,"end"],
                                                                                                                              max_spac = max_spacing,
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
                                                                                                                              SMFenzymeUsed = SMFenzymeUsed,
                                                                                                                              min_frag_data_len = min_frag_data_len,
                                                                                                                              min_frag_data_dens = min_frag_data_dens
                                                                                                                              )
                                                                                   
                                                                                   
                                                                                   #coocmatr <- outlist[["CoocCountTable"]]
                                                                                   
                                                                                   
                                                                                   outtbl <- tibble::tibble_row("SampleName" = bamgroupname,
                                                                                                                "seqnames" = regdf[1,"seqnames"],
                                                                                                                "start" = regdf[1,"start"],
                                                                                                                "end" = regdf[1,"end"],
                                                                                                                "strand" = regdf[1,"strand"],
                                                                                                                "names" = ifelse(!is.null(names(regions)), names(regions)[regi],regi),
                                                                                                                "nFragsFetched" = outlist[["nFragsFetched"]],
                                                                                                                "nFragsNonUnique" = outlist[["nFragsNonUnique"]],
                                                                                                                "nFragsBisFailed" = outlist[["nFragsBisFailed"]],
                                                                                                                "nFragsAnalyzed" = outlist[["nFragsAnalyzed"]],
                                                                                                                "ClippedUntilNbg" = outlist[["ClippedUntilNbg"]],
                                                                                                                "NotClippedProtectAbove" = outlist[["NotClippedProtectAbove"]],
                                                                                                                "nFragsClipped" = outlist[["nFragsClipped"]],
                                                                                                                "CoocCountTable" = list(outlist[["CoocCountTable"]]))
                                                                                   return(outtbl)
                                                                                 },mc.cores = ncores))
                                        
                                        return(ctbl)
                                      }))
  
  
}