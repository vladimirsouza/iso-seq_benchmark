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
#' @param method_names A vector of strings. The name of the methods to be compared.
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
#' @import ggplot
#' @import dplyr
#' 
#' @export
check_accuracy_per_coverage <- function(master_table, 
                                        method_names, data_name,
                                        truth_name, 
                                        interval_start, interval_end) {
  intervals <- matrix( c(interval_start, interval_end), byrow=TRUE, nrow=2) %>% 
    data.frame
  
  mt_intervalI_methodJ <- lapply(intervals, function(interv) {
    mt_intervalI <- paste0(data_name, "_coverage") %>% 
      master_table[,.data] %>% 
      { .data>=interv[1] & .data<=interv[2] } %>% 
      master_table[.data,]
    accur_intervalI_methodJ <- sapply(method_names, function(method_names_i) {
      calc_accuracy_measures(mt_intervalI, method_names_i, truth_name)
    })
    data.frame(accur_intervalI_methodJ) %>% 
      rownames_to_column(.data, "measure") %>% 
      cbind( interval=paste(.data, interv, collapse="-") )
  })
  
  mt_intervalI_methodJ <- bind_rows(mt_intervalI_methodJ) %>% 
    gather(.data, "method", "score", -measure, -interval, factor_key=TRUE)
  mt_intervalI_methodJ$measure <- factor(mt_intervalI_methodJ$measure)
  mt_intervalI_methodJ$interval <- sapply(unname(intervals), paste, collapse="-") %>% 
    factor(mt_intervalI_methodJ$interval, levels=.data, ordered=TRUE)
  
  ggplot(mt_intervalI_methodJ, aes(x=method, y=score, fill=method, colour=method)) +
    facet_grid(measure~interval) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 270))
}
