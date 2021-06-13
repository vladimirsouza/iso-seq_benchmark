#' Generate igv batch screenshots script
#'
#' @param chrm A vector of string. Chromossome name, e.g. "chrm1".
#' @param pos A vector of integers. Positions to center the screenshot.
#' @param output_dir A 1-length string. Directory where the IGV-batch-screenshot 
#'   script is saved.
#' @param prefix A 1-length string.
#' @param snapshot_path A 1-length string. Directory where the IGV-batch-screenshot 
#'   script will save the screenshots.
#' @param windows_size 1-length integer. Window size of the screenshot.
#' @param screenshot_number 1-length integer. Number of screenshot to generate.
#' @param output_positions If TRUE, output a data.frame with variant positions and 
#'   the reagions to visualize.
#'
#' @return NULL
#' 
#' @importFrom IRanges IRanges resize start end
#' @importFrom rlang .data
#' @importFrom utils write.table
#' @import dplyr
#' 
#' @export
igv_batch_screenshots <- function(chrm, pos, output_dir, prefix, snapshot_path, windows_size=1501, screenshot_number=100, output_positions=FALSE) {
  ### create the bed file that store the area of interest to visualize
  bed <- IRanges(pos, width=1) %>% 
    resize(windows_size, "center") 
  bed <- data.frame(chrm=chrm, start=start(bed), end=end(bed))
  
  ### subset from bed
  if( nrow(bed) < screenshot_number ) screenshot_number <- nrow(bed)
  k <- {1:nrow(bed)} %>% 
    sample(screenshot_number)
  bed <- bed[k,]
  
  ### output
  if(output_positions) {
    output_dataframe_positions <- data.frame( variants= paste0(chrm[k], ":", pos[k]),
                                              regions= paste0(bed$chrm, ":", bed$start, "-", bed$end) )
  }
  
  ### write bed file
  bed_file <- file.path( output_dir, paste0(prefix, ".bed") )
  write.table(bed,
              bed_file,
              sep="\t",
              quote=FALSE,
              row.names=FALSE,
              col.names=FALSE)
  
  ### convert the bed file into an IGV batch script
  batch_file <- file.path(output_dir, paste0(prefix, ".batch"))
  cmds <- gettextf("bedtools igv -path %s -i %s > %s",
                   file.path(snapshot_path, prefix), 
                   bed_file, batch_file)
  system(cmds)
  unlink(bed_file)
  
  if(output_positions) {
    output_dataframe_positions
  }
}







#' Calculate precision, sensitivity and F1-score of a method in a master table
#'
#' @param input_table A data.frame. The master table.
#' @param method_name A 1-length string. The name of the method to calculate the
#'   accuracy measures.
#' @param truth_name A 1-length string. The name of the ground-truth to validate
#'   the method.
#'
#' @return A named vector with the accuracy measures.
#' @export
calc_accuracy_measures <- function(input_table, method_name, truth_name) {
  in_method <- paste0("in_", method_name)
  in_truth <- paste0("in_", truth_name)
  method_classification <- paste0(method_name, "_classification")
  
  tp_count <- sum(input_table[,method_classification] == "TP")
  positive_method <- sum(input_table[,in_method])
  positive_truth <- sum(input_table[,in_truth])
  
  precision <- tp_count/positive_method
  sensitivity <- tp_count/positive_truth
  f1Score <- 2*(sensitivity * precision) / (sensitivity + precision)
  
  c(precision=precision, sensitivity=sensitivity, f1Score=f1Score)
}







#' Barplot of accuracy measures for different methods and coverage
#' 
#' Return a ggplot object for visualization.
#'
#' @param master_table A data.frame. The input master table.
#' @param method_names A vector of strings. The names of the methods to be compared
#'   in the master table.
#' @param output_method_names A vector of strings. How to output the names of the
#'   methods to be compared. The defaut is NULL.
#' @param data_name A 1-length string. The name of the dataset used with the methods 
#'   to be compared.
#' @param truth_name A 1-length string. The name of the ground-truth.
#' @param interval_start A vector of integers. The start position of all intervals.
#' @param interval_end A vector of integers. The end position of all intervals.
#'
#' @return A ggplot object.
#' 
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @import ggplot2
#' @import dplyr
#' 
#' @export
check_accuracy_per_coverage <- function(master_table, 
                                        method_names,
                                        output_method_names=NULL,
                                        data_name,
                                        truth_name, 
                                        interval_start, interval_end) {
  intervals <- matrix( c(interval_start, interval_end), byrow=TRUE, nrow=2) %>% 
    data.frame
  
  mt_intervalI_methodJ <- lapply(intervals, function(interv) {
    k <- master_table[ ,paste0(data_name, "_coverage") ]
    mt_intervalI <- master_table[ k>=interv[1] & k<=interv[2], ]
    
    accur_intervalI_methodJ <- sapply(method_names, function(method_names_i) {
      calc_accuracy_measures(mt_intervalI, method_names_i, truth_name)
    })
    
    k <- data.frame(accur_intervalI_methodJ)
    k <- rownames_to_column(k, "measure")
    cbind( k, interval=paste(interv, collapse="-") )
  })
  
  k <- bind_rows(mt_intervalI_methodJ) 
  mt_intervalI_methodJ <- gather(k, "method", "score", -.data$measure, -.data$interval, factor_key=TRUE)
  
  if( !is.null(output_method_names) ){
    if( length(method_names) != length(output_method_names) ){
      stop("The lengths of method_names and output_method_names must be equal.")
    }
    stopifnot( identical(method_names, levels(mt_intervalI_methodJ$method)) )
    levels(mt_intervalI_methodJ$method) <- output_method_names
  }
  
  mt_intervalI_methodJ$measure <- factor(mt_intervalI_methodJ$measure)
  
  k <- sapply(unname(intervals), paste, collapse="-")
  mt_intervalI_methodJ$interval <- factor(mt_intervalI_methodJ$interval, levels=k, ordered=TRUE)
  
  ggplot(mt_intervalI_methodJ, aes(x=.data$method, y=.data$score, fill=.data$method)) +
    facet_grid(.data$measure~.data$interval) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 270)) +
    theme(legend.position="bottom") +
    scale_x_discrete(labels=method_names)
  
}





#' Add a INFO tag from a VCF file to a master table
#' 
#' If a variant in the VCF file doesn'r contain the specified tag, the function
#' return -1 for that variant.
#' 
#' @param input_table A data.frame. The input master table.
#' @param col_name A 1-length string. the name of the column to add. If NULL (default),
#'   `col_name` is equal to `tag`.
#' @param tag A 1-length string. The name of the tag to add. It must be exactly how
#'   the tag is written in the VCF file (case sensitive).
#' @param vcf_file A 1-length string. The path to the VCF file from which the tag
#'   is extracted.
#'
#' @return A data.frame.
#' 
#' @importFrom vcfR read.vcfR
#' @importFrom dplyr left_join
#' 
#' @export
add_info_tag_from_vcf <- function(input_table, col_name=NULL, tag, vcf_file){
  if(is.null(col_name)){
    col_name <- tag
  }
  
  vcf <- read.vcfR(vcf_file)
  
  k <- gettextf(".*%s=(.+?);.*", tag)
  k <- sub( k, "\\1", 
            paste0(vcf@fix[,8], ";") )
  k <- as.numeric(k)
  k[is.na(k)] <- -1
  tag_values <- data.frame(chrm=vcf@fix[,1],
                           pos=as.integer(vcf@fix[,2]),
                           tag=k)
  names(tag_values)[3] <- col_name
  
  left_join(input_table, tag_values)
}



#' Add some statistic of a tag about reads of a BAM file
#'
#' @param input_table A data.frame. The input master table.
#' @param col_name A 1-length string. The name of the new column to add.
#' @param galn A GAlignment object. The alingnment from where the tag is 
#'   extracted.
#' @param tag A 1-length string. The name of tag that is wanted to calculate
#'   the statistic `stat`. The name must be exactly the same of one column
#'   name from `mcols(galn)`.
#' @param stat A function used to calculate the desired statistic.
#'
#' @return A data.frame.
#' 
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom GenomicAlignments findOverlaps
#' @importFrom S4Vectors mcols
#' 
#' @export
add_some_read_statistic_from_bam_metacolumn <- function(input_table, col_name, galn, tag, stat){
  dat_pos <- GRanges(input_table$chrm, IRanges(input_table$pos, width=1))
  ovl <- findOverlaps(dat_pos, galn)
  output_stat <- tapply( subjectHits(ovl), queryHits(ovl), function(i){
    stat( mcols(galn) [,tag] [i] )
  })
  input_table[,col_name] <- NA
  k <- as.integer( names(output_stat) )
  input_table[k ,col_name] <- unname(output_stat)
  
  input_table
}





#' Add column of density of variants around each variant called by a method
#'
#' @param input_table A data.frame. The input master table.
#' @param window_size A 1-length integer. The size of the window where variants
#'   are used to calculate the density of variants around a each variant called
#'   by the specified method.
#' @param method A 1-length string. The name of the method to be used.
#'
#' @return A data.frame.
#' 
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom IRanges resize
#' @importFrom IRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom snakecase to_lower_camel_case
#' 
#' @export
add_variant_density_of_a_method <- function(input_table, window_size, method){
  # input_table <- input_table_bk
  # window_size <- 101
  # method <- "tr_dna_merged"
  
  
  method_calls <- paste0("in_", method)
  in_method <- input_table[,method_calls] == 1
  in_method_range <- GRanges( input_table$chrm[in_method], 
                              IRanges(input_table$pos[in_method], width=1) )
  
  all_variants <- GRanges( input_table$chrm,
                           IRanges(input_table$pos, width=1) )
  all_windows <- resize(all_variants, window_size, "center")
  
  ovl <- findOverlaps(in_method_range, all_windows)
  variant_density <- tapply(queryHits(ovl), subjectHits(ovl), length)
  
  variant_density_column <- paste0( "variantDensity_", to_lower_camel_case(method) )
  input_table[,variant_density_column] <- 0
  k <- as.numeric( names(variant_density) )
  input_table[k, variant_density_column] <- variant_density
  
  input_table
}
