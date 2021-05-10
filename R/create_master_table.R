#' Initiate master table
#'
#' This function create the first columns of the master table. Those columns indicate
#' whether a method could call of not a variant. A master table is used to get
#' information form VCF files generated from different methods and compare to
#' ground-truth VCF files.
#'
#' @param ... VCF file addresses.
#' @param method_names Vector of strings of same length of number of input VCF files.
#'
#' @return A data.frame.
#' @export
initiate_master_table <- function(..., method_names) {
  ### load PASS variants
  vcf_file_list <- list(...)

  vcfs <- lapply(vcf_file_list, function(vcf_file) {
    vcf <- read.table(vcf_file)
    vcf[ vcf[[7]] == "PASS", ]
  })

  ### split VCFs by chromosomes
  present_chromosomes <- bind_rows(vcfs) %>%
    pull(1) %>%
    unique
  vcfs <- sapply(vcfs, function(vcf) {
    lapply(present_chromosomes, function(pchrm) {
      vcf[ vcf[,1]==pchrm, ]
    })
  })

  ### which variants are contained in each vcf and their DP tag (from VCF file)
  master_table <- apply(vcfs, 1, function(vcfs_chrmI) {
    all_positions <- bind_rows(vcfs_chrmI) %>%
      pull(2) %>%
      unique %>%
      sort
    in_chrmI_methodJ <- sapply(vcfs_chrmI, function(vcf_chrmI_methodJ) {
      1 * ( all_positions %in% vcf_chrmI_methodJ[[2]] )
    })

    dv_chrmI_methodJ <- sapply(vcfs_chrmI, function(vcf_chrmI_methodJ) {
      res <- strsplit(vcf_chrmI_methodJ[[10]], ":") %>%
        sapply("[", 3) %>%
        as.integer() %>%
        setNames(vcf_chrmI_methodJ[[2]])
      res <- unname( res[ as.character(all_positions) ] )
    })

    data.frame(chrom=vcfs_chrmI[[1]] [[1]] [1],
               in_chrmI_methodJ,
               dv_chrmI_methodJ)
  })
  master_table <- bind_rows(master_table)

  method_names <- c( paste0("in_", method_names),
                     paste0("dp_", method_names) )
  names(master_table)[-1] <- method_names

  master_table
}






# Helpers ------------------------------------------------------------------------



#' Get all splice sites positions of a BAM files
#'
#' Get the position of all splice sites in a BAM file and the number of reads that
#' support each one of them. Besides that, indicate whether the splice site is
#' acceptor or donor site.
#'
#' This function may take several hours to run.
#'
#' @param input_bam The input BAM file to extract splice site positions from.
#' @param threads Number of threads.
#'
#' @return A data.frame in which each row is a splice-site position.
#' @export
get_splice_sites_information <- function(input_bam, threads){
  ss <- extractAlignmentRangesOnReference( cigar(input_bam), start(input_bam),  )

  # add chromossome name
  chrm_name <- as.vector( seqnames(input_bam) )

  multicoreParam <- MulticoreParam(workers = threads)
  ss <- bpmapply(function(ss_i, chrm_name_i){
    if(length(ss_i) > 1){
      ss_start <- start(ss_i)[-1]
      ss_end <- head( end(ss_i), -1 )
      is_acceptor_site <- rep( c(1L:0L), c( length(ss_start), length(ss_end) ) )
      data.frame(chrm= chrm_name_i, pos= c(ss_start, ss_end),
                 is_acceptor_site)
    }else{
      NULL
    }
  }, ss, chrm_name, BPPARAM=multicoreParam)

  ss <- do.call(rbind, ss)
  if( is.null(ss) )
    return(NULL)
  ss <- group_by(ss, chrm, pos, is_acceptor_site) %>%
    tally

  ss
}
