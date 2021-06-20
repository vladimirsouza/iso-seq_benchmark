#' Make the flags of a BAM be the same as they were before SplitNCigarReads
#' 
#' This function is specially important to be used combined with DeepVariant
#'   to call variants from Iso-Seq data.
#' 
#' @param input_bam A 1-length string. Path of the BAM file before
#'   SplitNCigarReads.
#' @param input_sncr_bam A 1-length string. Path of the BAM file after
#'   SplitNCigarReads.
#' @param output_bam A 1-length string. File path to save the 
#'   output BAM file with corrected flags.
#' @param threads A 1-length integer. Number of threads.
#'
#' @return NULL.
#' 
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools ScanBamParam
#' 
#' @export
flagCorrection <- function(input_bam, input_sncr_bam, output_bam, threads){
  ### load, subset and sort bam's flags according sncr_bam's qnames.
  ### and save the flags to a temp file.
  sncr_bam_qname_order <- scanBam(input_sncr_bam,
                                  param=ScanBamParam(what="qname")) [[1]]$qname
  output_dir <- dirname(output_bam)
  ordered_flag_file <- tempfile(tmpdir=output_dir, fileext=".txt")
  k <- scanBam(input_bam,
               param=ScanBamParam(what=c("qname", "flag"))) [[1]]
  k <- setNames(k$flag, k$qname)
  k <- k[sncr_bam_qname_order]
  k <- unname(k)
  k <- data.frame(k)
  write.table(k, ordered_flag_file, col.names=FALSE, row.names=FALSE, quote=FALSE,
              sep="\t")
  
  ### replace the flag
  # write SAM body only (no headers)
  sncrBam_sam_body_file <- tempfile(tmpdir=output_dir, fileext=".sam")
  cmds <- gettextf("samtools view -@ %i -o %s %s", threads, sncrBam_sam_body_file,
                   input_sncr_bam)
  system(cmds)
  # replace the flag column
  # (https://stackoverflow.com/questions/7846476/replace-column-in-one-file-with-column-from-another-using-awk)
  sncrBam_sam_body_flagCorrection_file <- tempfile(tmpdir=output_dir, fileext=".sam")
  cmds <- gettextf("awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' OFS='\t' %s %s > %s",
                   ordered_flag_file, sncrBam_sam_body_file,
                   sncrBam_sam_body_flagCorrection_file)
  system(cmds)
  
  ### write the headers
  sncrBam_sam_flagCorrection_file <- tempfile(tmpdir=output_dir, fileext=".sam")
  cmds <- gettextf("samtools view -H -@ %i -o %s %s", threads,
                   sncrBam_sam_flagCorrection_file, input_sncr_bam)
  system(cmds)
  
  ### add body to headers
  cmds <- gettextf("cat %s >> %s", sncrBam_sam_body_flagCorrection_file,
                   sncrBam_sam_flagCorrection_file)
  system(cmds)
  
  ### convert sam to bam
  cmds <- gettextf("samtools view -S -b -@ %i -o %s %s", threads, output_bam,
                   sncrBam_sam_flagCorrection_file)
  system(cmds)
  
  unlink( c(ordered_flag_file,
            sncrBam_sam_body_file,
            sncrBam_sam_body_flagCorrection_file,
            sncrBam_sam_flagCorrection_file) )
}
