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
#' @param coverage_threshold A vector of integers. The minimum thresholds to filer 
#'   by read coverage.
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
                                        coverage_threshold){
  
  mt_thresholdI_methodJ <- lapply(coverage_threshold, function(coverage_threshold_i) {
    k <- master_table[ ,paste0(data_name, "_coverage") ]
    mt_thresholdI <- master_table[ k>=coverage_threshold_i, ]
    
    accur_thresholdI_methodJ <- sapply(method_names, function(method_names_i) {
      calc_accuracy_measures(mt_thresholdI, method_names_i, truth_name)
    })
    
    k <- data.frame(accur_thresholdI_methodJ)
    k <- rownames_to_column(k, "measure")
    cbind( k, threshold=paste(coverage_threshold_i, collapse="-") )
  })
  
  k <- bind_rows(mt_thresholdI_methodJ) 
  mt_thresholdI_methodJ <- gather(k, "method", "score", -.data$measure, -.data$threshold, factor_key=TRUE)
  
  if( !is.null(output_method_names) ){
    if( length(method_names) != length(output_method_names) ){
      stop("The lengths of method_names and output_method_names must be equal.")
    }
    stopifnot( identical(method_names, levels(mt_thresholdI_methodJ$method)) )
    levels(mt_thresholdI_methodJ$method) <- output_method_names
  }
  
  mt_thresholdI_methodJ$measure <- factor(mt_thresholdI_methodJ$measure)
  
  mt_thresholdI_methodJ$threshold <- factor(mt_thresholdI_methodJ$threshold, 
                                            levels=sort(coverage_threshold),
                                            ordered=TRUE)
  
  ggplot(mt_thresholdI_methodJ, aes(x=.data$method, y=.data$score, fill=.data$method)) +
    facet_grid(.data$measure~.data$threshold) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 270)) +
    theme(legend.position="bottom") +
    scale_x_discrete(labels=method_names)
  
}










#' Precision-recall curves to compare methods, using different filtering on read
#'   coverage
#'
#' @param master_table A data.frame. The input master table.
#' @param method_names A vector of strings. The names of the methods to be compared.
#' @param output_method_names A vector of strings. How to output the names of the
#'   methods to be compared. The default is NULL.
#' @param data_name A 1-length string. The name of the dataset used with the methods 
#'   to be compared.
#' @param truth_name A 1-length string. The name of the ground-truth.
#' @param coverage_thresholds A vector of integers. The minimum thresholds to filer 
#'   by read coverage. Each element defines a point in the precision-recall curves.
#' @param what A 1-length string. Possible values are: "snps_indels" (default), to 
#'   make curves for SNPs and indels separately; "snps", to make curves only for
#'   SNPs; "indels", to make curves only for indels; "overall", to make curves 
#'   without distinguishing SNPs and indels.
#' 
#' @return A ggplot object.
#' 
#' @import ggplot2
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @importFrom dplyr bind_rows
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' 
#' @export
precision_recall_curve_per_coverage <- function(master_table,
                                                method_names,
                                                output_method_names=NULL,
                                                data_name,
                                                truth_name,
                                                coverage_thresholds,
                                                what){
  
  if( !any(what %in% c("snps_indels", "snps", "indels", "overall")) ){
    stop("`what` argument must be either \"snps_indels\", \"snps\", \"indels\", or \"overall\"")
  }
  
  mt_thresholdI_methodJ <- lapply(coverage_thresholds, function(threshold_i) {
    k <- master_table[ ,paste0(data_name, "_coverage") ]
    mt_thresholdI <- master_table[ k>=threshold_i, ]
    
    if( what %in% c("snps_indels", "snps") ){
      snp_accur_thresholdI <- lapply(method_names, function(method_names_i){
        k <- paste0("is_indel_", method_names_i)
        k <- which( mt_thresholdI[,k] == 0)
        k <- mt_thresholdI[k,]
        calc_accuracy_measures(k, method_names_i, truth_name)
      })
      snp_accur_thresholdI <- do.call(rbind, snp_accur_thresholdI)
      snp_accur_thresholdI <- data.frame(snp_accur_thresholdI)
      snp_accur_thresholdI <- cbind(snp_accur_thresholdI,
                                    variant="snps",
                                    method=method_names,
                                    coverage_thresholds=threshold_i)
    }
    if( what %in% c("snps_indels", "indels") ){
      indel_accur_thresholdI <- lapply(method_names, function(method_names_i){
        k <- paste0("is_indel_", method_names_i)
        k <- which( mt_thresholdI[,k] == 1)
        k <- mt_thresholdI[k,]
        calc_accuracy_measures(k, method_names_i, truth_name)
      })
      indel_accur_thresholdI <- do.call(rbind, indel_accur_thresholdI)
      indel_accur_thresholdI <- data.frame(indel_accur_thresholdI)
      indel_accur_thresholdI <- cbind(indel_accur_thresholdI,
                                      variant="indels",
                                      method=method_names,
                                      coverage_thresholds=threshold_i)
    }
    if( what == "overall" ){
      overall_accur_thresholdI <- lapply(method_names, function(method_names_i){
        k <- mt_thresholdI
        calc_accuracy_measures(k, method_names_i, truth_name)
      })
      overall_accur_thresholdI <- do.call(rbind, overall_accur_thresholdI)
      overall_accur_thresholdI <- data.frame(overall_accur_thresholdI)
      overall_accur_thresholdI <- cbind(overall_accur_thresholdI,
                                        variant="overall",
                                        method=method_names,
                                        coverage_thresholds=threshold_i)
    }
    
    accur_thresholdI <- switch(
      what, 
      "snps_indels"={
        rbind(snp_accur_thresholdI, indel_accur_thresholdI)
      },
      "snps"={
        snp_accur_thresholdI
      },
      "indels"={
        indel_accur_thresholdI
      },
      "overall"={
        overall_accur_thresholdI
      }
    )
  })
  
  dat <- do.call(rbind, mt_thresholdI_methodJ) 
  names(dat) [names(dat) == "sensitivity"] <- "recall"
  dat$variant <- factor(dat$variant)
  dat$method <- factor(dat$method)
  dat$coverage_thresholds <- factor(dat$coverage_thresholds,
                                    levels=sort(coverage_thresholds),
                                    ordered=TRUE)
  if( !is.null(output_method_names) ){
    if( length(method_names) != length(output_method_names) ){
      stop("The lengths of `method_names` and `output_method_names` must be equal.")
    }
    stopifnot( identical(method_names, levels(dat$method)) )
    levels(dat$method) <- output_method_names
  }
  names(dat) [names(dat) == "coverage_thresholds"] <- "coverage >= n"
  
  
  p <- ggplot(dat, aes(.data$recall, .data$precision, fill=.data$method, colour=.data$method))
  if(what == "snps_indels"){
    p <- p +
      geom_point(aes(shape=.data$variant, size=.data$`coverage >= n`), alpha=.5) +
      geom_path(aes(linetype=.data$variant))
  }else{
    p <- p +
      geom_point(aes(size=.data$`coverage >= n`), alpha=.5) +
      geom_path()
  }
  p <- p +
    theme(text = element_text(size = 20))
  p
}








#' Calculate accuracy measures for multiple master tables
#' 
#' Create a data.frame with accuracy measures for multiple master tables.
#'
#' @param ... Master tables.
#' @param experiment_names Vector of strings. Name of the experiments in the same 
#'   order of the master tables.
#' @param method_names Vector or list of strings. Names of the methods to compare.
#'   If all experiments compare the same methods, `method_names` can be a single
#'   vector. Otherwise, `method_names` must be a list, in which each element 
#'   specifies the names of the methods to compare for each experiment.
#' @param output_method_names NULL (default), or a vector or a list of strings.
#'   Names of the methods to output. 
#' @param data_names Vector or list of strings.
#' @param truth_names Vector or list of strings.
#' @param coverage_thresholds Vector or list of integers.
#' @param what Vector or list of strings.
#'
#' @return A data.frame.
#' @export
calculate_precision_recall_for_multi_master_tables <- function(
  ...,
  experiment_names,
  method_names,
  output_method_names=NULL,
  data_names,
  truth_names,
  coverage_thresholds=c(3,15,40,100),
  what
){
  
  master_tables <- list(...)
  
  stopifnot( length(experiment_names) == length(master_tables) )
  
  if( length(master_tables) == 1 ){
    stopifnot( is.vector(method_names) & !is.list(method_names) )
    
    if( !is.null(output_method_names) ){
      stopifnot( is.vector(output_method_names) & !is.list(output_method_names) )
    }
    
    stopifnot( is.vector(data_names) & !is.list(data_names))
    stopifnot( length(data_names)==1 )
    
    stopifnot( is.vector(truth_names) & !is.list(truth_names) )
    stopifnot( length(truth_names)==1 )
    
    stopifnot( is.vector(coverage_thresholds) & !is.list(coverage_thresholds) )
    
    stopifnot( is.vector(what) & !is.list(what) )
    stopifnot( length(what)==1 )
    
    method_names <- list(method_names)
    output_method_names <- list(output_method_names)
    coverage_thresholds <- list(coverage_thresholds)
    
    mt_len <- 1
  }else{
    mt_len <- length(master_tables)
    
    if( is.null(output_method_names) ){
      output_method_names <- rep( list(output_method_names), mt_len )
    }else{
      if( is.vector(output_method_names) & !is.list(output_method_names) ){
        stopifnot( length(output_method_names) == length(method_names) )
        output_method_names <- rep( list(output_method_names), mt_len )
      }else{
        if( is.list(output_method_names) ){
          stopifnot( length(output_method_names) == length(master_tables) )
        }else{
          stop("Check argument output_method_names")
        }
      }
    }
    
    if( is.vector(method_names) & !is.list(method_names) ){
      method_names <- rep( list(method_names), mt_len )
    }else{
      if( is.list(method_names) ){
        stopifnot( length(method_names) == length(master_tables) )
      }else{
        stop("Check argument method_names")
      }
    }
    
    if( is.vector(data_names) & !is.list(data_names) ){
      if( length(data_names) == 1 ){
        data_names <- rep(data_names, mt_len)
      }else{
        stopifnot( length(data_names) == length(data_names) )
      }
    }else{
      stop("Check argument data_names")
    }
    
    if( is.vector(truth_names) & !is.list(truth_names) ){
      if( length(truth_names) == 1 ){
        truth_names <- rep(truth_names, mt_len)
      }else{
        stopifnot( length(truth_names) == length(truth_names) )
      }
    }else{
      stop("Check argument truth_names")
    }
    
    if( is.vector(coverage_thresholds) & !is.list(coverage_thresholds) ){
      coverage_thresholds <- rep( list(coverage_thresholds), mt_len )
    }else{
      if( is.list(coverage_thresholds) ){
        stopifnot( length(coverage_thresholds) == length(master_tables) )
      }else{
        stop("Check argument coverage_thresholds")
      }
    }
  }
  
  if( is.vector(what) & !is.list(what) ){
    if( length(what) == 1 ){
      what <- rep(what, mt_len)
    }else{
      stopifnot( length(what) == length(what) )
    }
  }else{
    stop("Check argument what")
  }
  
  
  
  dat_mts <- mapply(
    function(master_tables_l,
             method_names_l,
             output_method_names_l,
             data_names_l,
             truth_names_l,
             coverage_thresholds_l,
             what_l,
             experiment_names_l){
      
      # master_tables_l <- master_tables[[1]]
      # method_names_l <- method_names[[1]]
      # output_method_names_l <- output_method_names[[1]]
      # data_names_l <- data_names[[1]]
      # truth_names_l <- truth_names[[1]]
      # coverage_thresholds_l <- coverage_thresholds[[1]]
      # what_l <- what[[1]]
      # experiment_names_l <- experiment_names[[1]]
      
      if( !any(what_l %in% c("snps_indels", "snps", "indels", "overall")) ){
        stop("`what` argument must be either \"snps_indels\", \"snps\", \"indels\", or \"overall\"")
      }
      
      mt_thresholdI_methodJ <- lapply(coverage_thresholds_l, function(threshold_i) {
        k <- master_tables_l[ ,paste0(data_names_l, "_coverage") ]
        mt_thresholdI <- master_tables_l[ k>=threshold_i, ]
        
        if( what_l %in% c("snps_indels", "snps") ){
          snp_accur_thresholdI <- lapply(method_names_l, function(method_names_i){
            k <- paste0("is_indel_", method_names_i)
            k <- which(mt_thresholdI[,k] == 0)
            k <- mt_thresholdI[k,]
            calc_accuracy_measures(k, method_names_i, truth_names_l)
          })
          snp_accur_thresholdI <- do.call(rbind, snp_accur_thresholdI)
          snp_accur_thresholdI <- data.frame(snp_accur_thresholdI)
          snp_accur_thresholdI <- cbind(snp_accur_thresholdI,
                                        variant="snps",
                                        method=method_names_l,
                                        coverage_thresholds=threshold_i)
        }
        if( what_l %in% c("snps_indels", "indels") ){
          indel_accur_thresholdI <- lapply(method_names_l, function(method_names_i){
            k <- paste0("is_indel_", method_names_i)
            k <- which(mt_thresholdI[,k] == 1)
            k <- mt_thresholdI[k,]
            calc_accuracy_measures(k, method_names_i, truth_names_l)
          })
          indel_accur_thresholdI <- do.call(rbind, indel_accur_thresholdI)
          indel_accur_thresholdI <- data.frame(indel_accur_thresholdI)
          indel_accur_thresholdI <- cbind(indel_accur_thresholdI,
                                          variant="indels",
                                          method=method_names_l,
                                          coverage_thresholds=threshold_i)
        }
        if( what_l == "overall" ){
          overall_accur_thresholdI <- lapply(method_names_l, function(method_names_i){
            k <- mt_thresholdI
            calc_accuracy_measures(k, method_names_i, truth_names_l)
          })
          overall_accur_thresholdI <- do.call(rbind, overall_accur_thresholdI)
          overall_accur_thresholdI <- data.frame(overall_accur_thresholdI)
          overall_accur_thresholdI <- cbind(overall_accur_thresholdI,
                                            variant="overall",
                                            method=method_names_l,
                                            coverage_thresholds=threshold_i)
        }
        
        accur_thresholdI <- switch(
          what_l, 
          "snps_indels"={
            rbind(snp_accur_thresholdI, indel_accur_thresholdI)
          },
          "snps"={
            snp_accur_thresholdI
          },
          "indels"={
            indel_accur_thresholdI
          },
          "overall"={
            overall_accur_thresholdI
          }
        )
      })
      
      dat <- do.call(rbind, mt_thresholdI_methodJ) 
      names(dat) [names(dat) == "sensitivity"] <- "recall"
      dat$variant <- factor(dat$variant)
      dat$method <- factor(dat$method, levels=method_names_l)
      dat$coverage_thresholds <- factor(dat$coverage_thresholds,
                                        levels=sort(coverage_thresholds_l),
                                        ordered=TRUE)
      if( !is.null(output_method_names_l) ){
        if( length(method_names_l) != length(output_method_names_l) ){
          stop("The lengths of `method_names` and `output_method_names` must be equal.")
        }
        stopifnot( identical(method_names_l, levels(dat$method)) )
        levels(dat$method) <- output_method_names_l
      }
      names(dat) [names(dat) == "coverage_thresholds"] <- "coverage >= n"
      
      dat <- cbind(dat, experiment=experiment_names_l)
      
    },
    master_tables,
    method_names,
    output_method_names,
    data_names,
    truth_names,
    coverage_thresholds,
    what,
    experiment_names,
    SIMPLIFY=FALSE
  )
  dat_mts <- do.call(rbind, dat_mts)
  
  dat_mts
  
  # p <- ggplot(dat_mts, aes(.data$recall, .data$precision,
  #                          fill=.data$method, colour=.data$method)) +
  #   facet_grid(experiment~.)    ####### <<<<<--------=== not using .data$
  # if(what == "snps_indels"){
  #   p <- p +
  #     geom_point(aes(shape=.data$variant, size=.data$`coverage >= n`), alpha=.5) +
  #     geom_path(aes(linetype=.data$variant))
  # }else{
  #   p <- p +
  #     geom_point(aes(size=.data$`coverage >= n`), alpha=.5) +
  #     geom_path()
  # }
  # p <- p +
  #   coord_fixed(ratio=1) +
  #   theme(text = element_text(size = 20))
  # 
  # p
}




#' For each N-cigar-read-count interval, calculate accuracy measures for 
#'   multiple master tables
#' 
#' Create a data.frame with accuracy measures for multiple master tables and 
#'   different intervals for N-cigar-read counts.
#'
#' @param ... Master tables.
#' @param experiment_names Vector of strings. Name of the experiments in the same 
#'   order of the master tables.
#' @param method_names Vector or list of strings. Names of the methods to compare.
#'   If all experiments compare the same methods, `method_names` can be a single
#'   vector. Otherwise, `method_names` must be a list, in which each element 
#'   specifies the names of the methods to compare for each experiment.
#' @param output_method_names NULL (default), or a vector or a list of strings.
#'   Names of the methods to output. 
#' @param data_names Vector or list of strings.
#' @param truth_names Vector or list of strings.
#' @param start_n_cigar_read_percent_intervals vector of doubles. The start of
#'   each interval of N-cigar-read counts.
#' @param end_n_cigar_read_percent_intervals vector of doubles. The end of
#'   each interval of N-cigar-read counts.
#' @param what Vector or list of strings.
#'
#' @return A data.frame.
#' @export
calculate_precision_recall_for_n_cigar_read_count_intervarls <- function(
  ...,
  experiment_names,
  method_names,
  output_method_names=NULL,
  data_names,
  truth_names,
  start_n_cigar_read_percent_intervals=c(0, 0.05, 0.9),
  end_n_cigar_read_percent_intervals=c(0.05, 0.9, 1),
  what
){
  
  
  # master_tables <- list(dat_jurkat_cover10to20, dat_wtc11_cover10to20)
  # experiment_names <- c("jurkat", "wtc11")
  # method_names <- c("dv", "dv_s", "dv_s_fc", "gatk_s")
  # output_method_names = c("DV", "SNCR+DV", "SNCR+FC+DV", "SNCR+GATK")
  # data_names <- "isoSeq"
  # truth_names <- c("jurkat_dna_merged", "allen")
  # start_n_cigar_read_percent_intervals <- c(0, 0.05, 0.9)
  # end_n_cigar_read_percent_intervals <- c(0.05, 0.9, 1)
  # what <- "snps_indels"
  
  
  n_cigar_read_percent_intervals <- matrix(
    c(start_n_cigar_read_percent_intervals,
      end_n_cigar_read_percent_intervals),
    ncol=2
  )
  
  
  master_tables <- list(...)
  
  stopifnot( length(experiment_names) == length(master_tables) )
  
  if( length(master_tables) == 1 ){
    stopifnot( is.vector(method_names) & !is.list(method_names) )
    
    if( !is.null(output_method_names) ){
      stopifnot( is.vector(output_method_names) & !is.list(output_method_names) )
    }
    
    stopifnot( is.vector(data_names) & !is.list(data_names))
    stopifnot( length(data_names)==1 )
    
    stopifnot( is.vector(truth_names) & !is.list(truth_names) )
    stopifnot( length(truth_names)==1 )
    
    # stopifnot( is.vector(n_cigar_read_percent_intervals) & !is.list(n_cigar_read_percent_intervals) )
    
    stopifnot( is.vector(what) & !is.list(what) )
    stopifnot( length(what)==1 )
  }else{
    mt_len <- length(master_tables)
    
    if( is.null(output_method_names) ){
      output_method_names <- rep( list(output_method_names), mt_len )
    }else{
      if( is.vector(output_method_names) & !is.list(output_method_names) ){
        stopifnot( length(output_method_names) == length(method_names) )
        output_method_names <- rep( list(output_method_names), mt_len )
      }else{
        if( is.list(output_method_names) ){
          stopifnot( length(output_method_names) == length(master_tables) )
        }else{
          stop("Check argument output_method_names")
        }
      }
    }
    
    if( is.vector(method_names) & !is.list(method_names) ){
      method_names <- rep( list(method_names), mt_len )
    }else{
      if( is.list(method_names) ){
        stopifnot( length(method_names) == length(master_tables) )
      }else{
        stop("Check argument method_names")
      }
    }
    
    if( is.vector(data_names) & !is.list(data_names) ){
      if( length(data_names) == 1 ){
        data_names <- rep(data_names, mt_len)
      }else{
        stopifnot( length(data_names) == length(data_names) )
      }
    }else{
      stop("Check argument data_names")
    }
    
    if( is.vector(truth_names) & !is.list(truth_names) ){
      if( length(truth_names) == 1 ){
        truth_names <- rep(truth_names, mt_len)
      }else{
        stopifnot( length(truth_names) == length(truth_names) )
      }
    }else{
      stop("Check argument truth_names")
    }
    
    # if( is.vector(n_cigar_read_percent_intervals) & !is.list(n_cigar_read_percent_intervals) ){
    #   n_cigar_read_percent_intervals <- rep( list(n_cigar_read_percent_intervals), mt_len )
    # }else{
    #   if( is.list(n_cigar_read_percent_intervals) ){
    #     stopifnot( length(n_cigar_read_percent_intervals) == length(master_tables) )
    #   }else{
    #     stop("Check argument n_cigar_read_percent_intervals")
    #   }
    # }
    n_cigar_read_percent_intervals <- rep(
      list(n_cigar_read_percent_intervals),
      mt_len
    )
  }
  
  if( is.vector(what) & !is.list(what) ){
    if( length(what) == 1 ){
      what <- rep(what, mt_len)
    }else{
      stopifnot( length(what) == length(what) )
    }
  }else{
    stop("Check argument what")
  }
  
  
  
  dat_mts <- mapply(
    function(master_tables_l,
             method_names_l,
             output_method_names_l,
             data_names_l,
             truth_names_l,
             n_cigar_read_percent_intervals_l,
             what_l,
             experiment_names_l){
      
      # master_tables_l <- master_tables[[2]]
      # method_names_l <- method_names[[2]]
      # output_method_names_l <- output_method_names[[2]]
      # data_names_l <- data_names[[2]]
      # truth_names_l <- truth_names[[2]]
      # n_cigar_read_percent_intervals_l <- n_cigar_read_percent_intervals[[2]]
      # what_l <- what[[2]]
      # experiment_names_l <- experiment_names[[2]]
      
      
      
      
      if( !any(what_l %in% c("snps_indels", "snps", "indels", "overall")) ){
        stop("`what` argument must be either \"snps_indels\", \"snps\", \"indels\", or \"overall\"")
      }
      
      mt_intervalI_methodJ <- apply(n_cigar_read_percent_intervals_l, 1, function(interval_i) {
        # interval_i <- n_cigar_read_percent_intervals_l[1,]
        
        k <- master_tables_l[ ,paste0("percent_n_cigar_reads") ]
        mt_intervalI <- master_tables_l[ k>=interval_i[1] & k<interval_i[2], ]
        
        if( what_l %in% c("snps_indels", "snps") ){
          snp_accur_thresholdI <- lapply(method_names_l, function(method_names_i){
            # method_names_i <- method_names_l[[1]]
            
            k <- which( mt_intervalI$is_indel_dv_s_fc==0 )
            k <- mt_intervalI[k,]
            calc_accuracy_measures(k, method_names_i, truth_names_l)
          })
          snp_accur_thresholdI <- do.call(rbind, snp_accur_thresholdI)
          snp_accur_thresholdI <- data.frame(snp_accur_thresholdI)
          interval_i <- paste0("[", interval_i[1], "-", interval_i[2], ")")
          snp_accur_thresholdI <- cbind(snp_accur_thresholdI,
                                        variant="snps",
                                        method=method_names_l,
                                        n_cigar_read_percent_intervals=interval_i)
        }
        if( what_l %in% c("snps_indels", "indels") ){
          indel_accur_thresholdI <- lapply(method_names_l, function(method_names_i){
            k <- which( mt_intervalI$is_indel_dv_s_fc==1 )
            k <- mt_intervalI[k,]
            calc_accuracy_measures(k, method_names_i, truth_names_l)
          })
          indel_accur_thresholdI <- do.call(rbind, indel_accur_thresholdI)
          indel_accur_thresholdI <- data.frame(indel_accur_thresholdI)
          indel_accur_thresholdI <- cbind(indel_accur_thresholdI,
                                          variant="indels",
                                          method=method_names_l,
                                          n_cigar_read_percent_intervals=interval_i)
        }
        if( what_l == "overall" ){
          overall_accur_thresholdI <- lapply(method_names_l, function(method_names_i){
            k <- mt_intervalI
            calc_accuracy_measures(k, method_names_i, truth_names_l)
          })
          overall_accur_thresholdI <- do.call(rbind, overall_accur_thresholdI)
          overall_accur_thresholdI <- data.frame(overall_accur_thresholdI)
          overall_accur_thresholdI <- cbind(overall_accur_thresholdI,
                                            variant="overall",
                                            method=method_names_l,
                                            n_cigar_read_percent_intervals=interval_i)
        }
        
        accur_thresholdI <- switch(
          what_l, 
          "snps_indels"={
            rbind(snp_accur_thresholdI, indel_accur_thresholdI)
          },
          "snps"={
            snp_accur_thresholdI
          },
          "indels"={
            indel_accur_thresholdI
          },
          "overall"={
            overall_accur_thresholdI
          }
        )
      })
      
      dat <- do.call(rbind, mt_intervalI_methodJ) 
      names(dat) [names(dat) == "sensitivity"] <- "recall"
      dat$variant <- factor(dat$variant)
      dat$method <- factor(dat$method)
      
      
      k <- apply(n_cigar_read_percent_intervals_l, 1, function(x){
        paste0()
      })
      dat$n_cigar_read_percent_intervals <- factor(dat$n_cigar_read_percent_intervals,
                                                   levels=sort(unique(dat$n_cigar_read_percent_intervals)),
                                                   ordered=TRUE)
      if( !is.null(output_method_names_l) ){
        if( length(method_names_l) != length(output_method_names_l) ){
          stop("The lengths of `method_names` and `output_method_names` must be equal.")
        }
        stopifnot( identical(method_names_l, levels(dat$method)) )
        levels(dat$method) <- output_method_names_l
      }
      names(dat) [names(dat) == "n_cigar_read_percent_intervals"] <- "% N-cigar reads"
      
      dat <- cbind(dat, experiment=experiment_names_l)
      
    },
    master_tables,
    method_names,
    output_method_names,
    data_names,
    truth_names,
    n_cigar_read_percent_intervals,
    what,
    experiment_names,
    SIMPLIFY=FALSE
  )
  dat_mts <- do.call(rbind, dat_mts)
  
  dat_mts
  
  # p <- ggplot(dat_mts, aes(.data$recall, .data$precision,
  #                          fill=.data$method, colour=.data$method)) +
  #   facet_grid(experiment~.)    ####### <<<<<--------=== not using .data$
  # if(what == "snps_indels"){
  #   p <- p +
  #     geom_point(aes(shape=.data$variant, size=.data$`coverage >= n`), alpha=.5) +
  #     geom_path(aes(linetype=.data$variant))
  # }else{
  #   p <- p +
  #     geom_point(aes(size=.data$`coverage >= n`), alpha=.5) +
  #     geom_path()
  # }
  # p <- p +
  #   coord_fixed(ratio=1) +
  #   theme(text = element_text(size = 20))
  # 
  # p
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







#' Add a FORMAT tag from a VCF file to a master table
#' 
#' The function adds to `input_table` the new column <tagName>_<methodName>,
#'   where <tagName> is `tag` in lower case and <methodName> is `method_name`.
#' 
#' @param input_table A data.frame. The input master table.
#' @param method_name A 1-length string. The name of the method.
#' @param tag A 1-length string. The name of the tag to add. It must be exactly
#'   how the tag is written in the VCF file (case sensitive).
#' @param vcf_file A 1-length string. The path to the VCF file from which the tag
#'   is extracted.
#'
#' @return A data.frame.
#' 
#' @importFrom vcfR read.vcfR
#' @importFrom dplyr left_join
#' 
#' @export
add_format_tag_from_vcf <- function(input_table, method_name, tag, vcf_file){
  
  vcf <- read.vcfR(vcf_file)
  
  cmds <- gettextf("bcftools query -f '%%CHROM %%POS[ %%%s]\n' %s",
                   tag, vcf_file)
  x <- system(cmds, intern=TRUE)
  x <- strsplit(x, " ")
  x <- do.call(rbind, x)
  stopifnot( all(vcf@fix[,1:2] == x[,-3]) )
  x <- data.frame(x)
  col_name <- paste( tolower(tag), method_name, sep="_")
  names(x) <- c("chrm", "pos", col_name)
  x$pos <- as.integer(x$pos)
  res <- left_join(input_table, x)
  stopifnot( identical(input_table[,1:2], res[,1:2]) )
  
  res
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
#' @param used_methods A vector of strings. The names of the methods to be used.
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
add_variant_density_of_a_method <- function(input_table, window_size, used_methods){
  # input_table <- dat
  # window_size <- 101
  # used_methods <- "dv_s_fc"
  # used_methods <- c("dv_s_fc", "gatk_s")
  
  
  
  # method_calls <- paste0("in_", method)
  # in_method <- input_table[,method_calls] == 1
  # in_method_range <- GRanges( input_table$chrm[in_method], 
  #                             IRanges(input_table$pos[in_method], width=1) )
  methods_calls <- paste0("in_", used_methods)
  k <- input_table[,methods_calls, drop=FALSE]
  
  k <- rowSums(k)
  in_methods <- k > 0
  in_methods_range <- GRanges( input_table$chrm[in_methods],
                               IRanges(input_table$pos[in_methods], width=1) )
  
  
  
  
  all_variants <- GRanges( input_table$chrm,
                           IRanges(input_table$pos, width=1) )
  all_windows <- resize(all_variants, window_size, "center")
  
  ovl <- findOverlaps(in_methods_range, all_windows)
  variant_density <- tapply(queryHits(ovl), subjectHits(ovl), length)
  
  variant_density_column <- paste(
    c( "variantDensity", to_lower_camel_case(used_methods) ),
    collapse="_"
  )
  input_table[,variant_density_column] <- 0
  k <- as.numeric( names(variant_density) )
  input_table[k, variant_density_column] <- variant_density
  
  input_table
}






#' Standardize genotypes
#' 
#' This function rewrite genotype like:
#' * 0|1 to 0/1;
#' * 0/2, 0/3, ..., to 0/1;
#' * 2/2, 3/3, ..., to 1/1;
#' * 1/3, 2/3, ..., to 1/2
#'
#' @param gt A vector of strings. The input genotypes.
#' 
#' @importFrom stringr str_split
#'
#' @return A vector of strings.
#' @export
standardize_genotype <- function(gt){
  gt_sd <- sub("\\|", "/", gt)
  gt_sd_split <- str_split(gt_sd, "/", simplify=TRUE)
  # heterozigous
  het <- apply(gt_sd_split=="0", 1, sum)
  het <- het==1
  gt_sd[het] <- "0/1"
  gt_sd_split[het, 2] <- "1"
  # homozigous alternative
  homAlt <- gt_sd_split[,1] != "0" & gt_sd_split[,1]==gt_sd_split[,2]
  gt_sd[homAlt] <- "1/1"
  gt_sd_split[homAlt, 1] <- "1"
  gt_sd_split[homAlt, 2] <- "1"
  # heterozygous alternative
  hetRef <- gt_sd=="0/0"
  hetAlt <- !hetRef & !het & !homAlt
  gt_sd[hetAlt] <- "1/2"
  gt_sd
}




#' get the genotype and the variant type of the calls of a primary and 
#' secundary method
#' 
#' Create a new data.frame from a input master table. The rows correspond to
#'   the same variants in `input_table` and in the same order.
#'
#' @param input_table A data.frame. The input master table.
#' @param gt_first A 1-length string. The name of the column that stores the
#'   genotypes called by a method. This method is called the first method.
#' @param gt_second A 1-length string. If the first method does not call some 
#'   variants, extract their genotype from column `gt_second`.
#' @param vt_first A 1-length string. The name of the column that stores the
#'   variant type called by a method. This method is called the first method.
#' @param vt_second A 1-length string. If the first method does not call some 
#'   variants, extract their variant type from column `vt_second`.
#'
#' @return A data.frame.
#' 
#' @export
gt_vt_method <- function(input_table, gt_first, gt_second, vt_first, vt_second){
  
  ### genotype (gt)
  gt <- input_table[,gt_first]
  k <- is.na(gt) | gt=="0/0"
  gt[k] <- input_table[,gt_second] [k]
  
  ### gt -- standardized 
  gt_sd <- standardize_genotype(gt)
  
  ### variant type (vt)
  vt <- input_table[,vt_first]
  k <- is.na(vt)
  vt[k] <- input_table[,vt_second] [k]
  
  data.frame(gt, gt_sd, vt)
}











#' Get indels in homopolymers
#' 
#' Create a new data.frame that stores deletions of inside homopolymers, and
#' insertions of a same nucleotide type that happen inside homopolymers of
#' that same nucleotide type. Homopolymer length equal to 1 means
#' non-homopolymers.
#'
#' @param input_table A data.frame. The input master table.
#' @param first_method_name A 1-length string of the name of a method. First,
#'   the variant information that is taken is related to this method. If a
#'   variant is not called by the method, its information is taken relating
#'   to the method specifiec in `second_method_name`.
#' @param second_method_name A 1-length string. The name of a method.
#' @param vcf_first A 1-length string. The path of the VCF file of the first
#'   method.
#' @param vcf_second A 1-length string. The path of the VCF file of the
#'   second method.
#' @param method_dataset_name A 1-length string. Name of the dataset used to
#'   call variants by the methods to be compared.
#' @param homopolymers A CompressedIRangesList object. It should store all 
#'   homopolymers, it's nucleotive types and lengths, of the genome used as
#'   the reference to call the variants. It is gerated by the function
#'   `sarlacc::homopolymerFinder`.
#' @param ref_fasta_seqs A DNAStringSet object. The sequences of the genome
#'   used as the reference to call the variants. It's names must be in the
#'   form like "chr1", "chr2", ..., "chrX", "chrY".
#' @param min_isoseq_coverage min iso-seq read coverage to filter.
#' @param genotyped_alt One of the strings "find" or "same". If "find", the
#'   function finds the correct alternative allele based on the genotype. This
#'   option must be chosen if there are multiple values for column ALT, which
#'   is the case of VCF files output by DeepVariant and Clair3. Since GATK's
#'   VCF files keep only the genotyped alternative allele in column ALT, this
#'   argument should be "same".
#'
#' @return A data.frame.
#' 
#' @importFrom vcfR read.vcfR
#' @importFrom dplyr left_join
#' @importFrom Biostrings readDNAStringSet
#' @importFrom IRanges IRanges
#' @importFrom XVector subseq compact
#' 
#' @export
method_homopolymer_indels <- function(input_table, first_method_name, second_method_name,
                                      vcf_first, vcf_second, method_dataset_name, homopolymers,
                                      ref_fasta_seqs, min_isoseq_coverage, genotyped_alt){
  
  stopifnot( genotyped_alt %in% c("find", "same") )
  k <- grep("dv|c3|clair3", second_method_name, ignore.case=TRUE)
  if( identical(k,1L) ){
    if(genotyped_alt!="find"){
      message("It looks like the second method is DeepVariant or Clair3, but genotyped_alt is not 'find'.")
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }else{
    k <- grep("gatk", second_method_name, ignore.case=TRUE)
    if( identical(k,1L) ){
      if(genotyped_alt!="same"){
        message("It looks like the second method is GATK, but genotyped_alt is not 'same'.")
        invisible(readline(prompt="Press [enter] to continue"))
      }
    }
  }
  
  ### genotype (gt)
  gt_first <- paste0("gt_", first_method_name)
  gt_second <- paste0("gt_", second_method_name)
  gt <- input_table[,gt_first]
  het <- is.na(gt) | gt=="0/0"
  gt[het] <- input_table[,gt_second] [het]
  
  ### gt -- standardized 
  gt_sd <- standardize_genotype(gt)
  
  ### variant type (vt)
  vt_first <- paste0("variantType_", first_method_name)
  vt_second <- paste0("variantType_", second_method_name)
  vt <- input_table[,vt_first]
  vt[het] <- input_table[,vt_second] [het]
  
  ### add ref/alt alleles from `vcf_first` to `input_table`
  vcf_fst <- read.vcfR(vcf_first)
  k <- vcf_fst@fix[ ,c("CHROM", "POS", "REF", "ALT") ]
  k <- data.frame(k)
  k <- rename(k, "chrm"="CHROM", "pos"="POS", "ref_fst"="REF", "alt_fst"="ALT")
  k$pos <- as.integer(k$pos)
  input_table <- left_join(input_table, k)
  
  ### add ref/alt alleles from `vcf_second` to `input_table`
  vcf_snd <- read.vcfR(vcf_second)
  k <- vcf_snd@fix[ ,c("CHROM", "POS", "REF", "ALT") ]
  k <- data.frame(k)
  k <- rename(k, "chrm"="CHROM", "pos"="POS", "ref_snd"="REF", "alt_snd"="ALT")
  k$pos <- as.integer(k$pos)
  input_table <- left_join(input_table, k)
  
  ### ref/alt alleles
  ref_alt_alleles <- input_table[, c("ref_fst", "alt_fst")]
  ref_alt_alleles[het,] <- input_table[het, c("ref_snd", "alt_snd")]
  names(ref_alt_alleles) <- c("ref_allele", "alt_allele")
  
  ### make a data frame and filter
  k <- paste0( "in_",  c(first_method_name, second_method_name) )
  k <- c( k, paste0(second_method_name, "_classification") )
  coverage_dataset_column <- paste0(method_dataset_name, "_coverage")
  k <- c("chrm", "pos", "homopolymer_length_indel", coverage_dataset_column, k)
  res <- cbind( input_table[,k], gt, gt_sd, vt, ref_alt_alleles )
  k <- (res$gt_sd %in% c("0/1", "1/1")) & (res$vt %in% c("insertion", "deletion"))
  res <- res[k,]
  
  ### filter by iso-seq coverage
  res <- res[ res[,coverage_dataset_column] >= min_isoseq_coverage, ]
  
  ### genotyped aternative allele
  if(genotyped_alt=="find"){
    k <- sub("^[0-9]+/", "", res$gt_sd)
    k <- as.integer(k)
    alt_split <- strsplit(res$alt_allele, ",")
    alt_genotyped <- mapply( function(a, ki){
      a[ki]
    }, alt_split, k)
    res <- cbind(res, alt_genotyped)
  }else{
    res$alt_genotyped <- res$alt_allele
  }
  
  ### there might be a mix between insertions and deletions. to make it simple,
  ### removed all those cases
  res <- res[ !(res$vt=="deletion" & nchar(res$alt_genotyped)>1), ]
  res <- res[ !(res$vt=="insertion" & nchar(res$ref_allele)>1), ]
  
  ### deleted nts must be the same
  res$ref_nts <- sub(".", "", res$ref_allele)
  res$alt_nts <- sub(".", "", res$alt_genotyped)
  
  k <- strsplit(res$ref_nts, NULL)
  k <- sapply(k, function(x){
    unique(x)
  })
  k <- lengths(k)
  stopifnot( all(k[res$vt=="insertion"]==0) )
  res <- res[ !(res$vt=="deletion" & k!=1), ]
  
  ### inserted nts must be the same
  k <- strsplit(res$alt_nts, NULL)
  k <- sapply(k, function(x){
    unique(x)
  })
  k <- lengths(k)
  stopifnot( all(k[res$vt=="deletion"]==0) )
  res <- res[ !(res$vt=="insertion" & k!=1), ]
  
  ### nts of indels and homopolymers (from genome of reference) must be the same
  # add nt type of homopolymers
  add_homopolymer_nt_when_indels <- add_homopolymer_length_when_indels
  res <- add_homopolymer_nt_when_indels(res, homopolymers, ouput_what="nts")
  # add nt type of non-homopolymers
  res_split <- split(res, res$chrm)
  ref_fasta_seqs <- ref_fasta_seqs[names(res_split)]
  res_split <- lapply( seq_along(res_split), function(i){
    r <- res_split[[i]]
    s <- ref_fasta_seqs[i]
    is_not_hom <- is.na(r$homopolymer_nt_indel)
    k_pos <- r$pos[is_not_hom] +1
    s <- rep(s, length(k_pos))
    nts <- compact(
      subseq(
        s,
        k_pos,
        k_pos
      )
    )
    r$homopolymer_nt_indel[is_not_hom] <- unname( as.vector(nts) )
    r
  })
  res <- do.call(rbind, res_split)
  
  ### filter out insertions that the nts of the alternative allele are not the
  ### same of the nts of the homopolymer in the reference fasta
  k <- substring(res$alt_nts, 1, 1)
  k <- res$vt=="insertion" & res$homopolymer_nt_indel!=k
  res <- res[!k,]
  
  res
}











#' Make plot to compare accuracy per homopolymer length
#' 
#' The function uses the output from function `method_homopolymer_indels`
#'   to make a plot to compare how the accuracy vary according to the
#'   homopolymer length or not within homopolymers. Those accuracy measures
#'   can be either rates of TPs, FNs and FPs, or the precision, recall
#'   and F1-score.
#'
#' @param input_hom_table A data.frame. The output of the function 
#'   `method_homopolymer_indels`.
#' @param variant_type A 1-length string. The variant type. Possivle
#'   values are: "snp", "deletion", or "insertion".
#' @param method_name A 1-length string. The name of the method to analyse.
#' @param truth_name A 1-length string. The name of the ground-truth.
#' @param hom_length_intervals A vector of integers. Must be the same
#'   length of `interval_names`. The inferior limit for each interval for
#'   homopolymer length.
#' @param interval_names A vector of strings. Must be the same length of
#'   `hom_length_intervals`. The name of each interval of homopolymer 
#'   length.
#' @param to_calculate A 1-length string. Possible values are: "rates" or
#'   "pre_rec_f1".
#'
#' @return A list, in which the first element, named `p`, is a ggplot
#'   object, and the second element, named `interval_counts`, is a named
#'   vector of integers that stores the counts of variants in each 
#'   interval specified by the x-axis of the plot in `p`.
#' 
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr rename
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr summarise
#' @importFrom tidyr gather
#' @import ggplot2
#' @importFrom rlang .data
#' 
#' @export
make_homopolymer_plot <- function(input_hom_table, variant_type,
                                  method_name, truth_name, hom_length_intervals,
                                  interval_names, to_calculate){
  
  k <- to_calculate %in% c("rates", "pre_rec_f1")
  if(!k){
    stop("'to_calculate' argument must be either 'rates' or 'pre_rec_f1'.")
  }
  
  input_hom_table <- input_hom_table[input_hom_table$vt==variant_type,]
  
  k <- findInterval(input_hom_table$homopolymer_length_indel,
                    hom_length_intervals)
  stopifnot( !any(k==0) )
  k <- interval_names[k]
  k <- factor(k, levels=interval_names, ordered=TRUE)
  input_hom_table$homopolymer_length_intervals <- k
  
  method_class_name <- paste0(method_name, "_classification")
  k <- input_hom_table[,method_class_name] %in% c("FN", "TP", "FP")
  input_hom_table <- input_hom_table[k,]
  
  if(to_calculate=="rates"){
    k <- rename(input_hom_table, "Classification"=all_of(method_class_name))
    k <- mutate(k, Classification=droplevels(.data$Classification))
    k <- group_by(k, .data$homopolymer_length_intervals, .data$Classification)
    class_counts <- summarise(k, count=n())
    k <- group_by(class_counts, .data$homopolymer_length_intervals)
    total_counts <- summarise(k, total_count=sum(.data$count))
    k <- left_join(class_counts, total_counts)
    class_counts <- mutate(k, percent=.data$count/.data$total_count)
  }else{
    input_hom_table_split <- split(input_hom_table,
                                   input_hom_table$homopolymer_length_intervals)
    k <- lapply(input_hom_table_split, function(x){
      y <- calc_accuracy_measures(x, method_name, truth_name)
      c(y, total_count=nrow(x))
    })
    k <- do.call(rbind, k)
    k <- data.frame(k)
    k <- rownames_to_column(k, "homopolymer_length_intervals")
    k$homopolymer_length_intervals <- factor(k$homopolymer_length_intervals,
                                             levels=interval_names, ordered=TRUE)
    class_counts <- gather(k, "Classification", "percent", .data$precision:.data$f1Score,
                           factor_key=TRUE)
  }
  
  k <- split(class_counts$total_count, class_counts$homopolymer_length_intervals)
  interval_counts <- sapply(k, unique)
  stopifnot( is.atomic(interval_counts) )
  dat_text <- data.frame(
    label=paste0("n=", interval_counts),
    x=names(interval_counts),
    y=1.05*max(class_counts$percent),
    Classification=NA
  )
  
  p <- ggplot(class_counts, aes(x=.data$homopolymer_length_intervals, y=.data$percent,
                                group=.data$Classification, colour=.data$Classification)) +
    geom_point() +
    geom_line() +
    xlab("Variant is in a homopolymer of length n") +
    ylab("Proportion of each classification") +
    # ggtitle("Classifications from SNCR+FC+DeepVariant") +
    theme(text = element_text(size=18)) +
    # ylim(0,1) +
    geom_text( data=dat_text, mapping= aes(x=x, y=y, label=label) ) +
    NULL
  p
}







#' Organize the data used to make plots about homopolymer analysis
#' 
#' The user may want to use this function several times to pull information about
#'   different combinations between `variant_type` and `method_name` from
#'   `input_hom_table`. In the future, the function should make the job automatically
#'   using loops.
#'
#' @param input_hom_table A data.frame generated by function `method_homopolymer_indels`.
#' @param variant_type A 1-length string. Possible values are "insertion" or "deletion".
#' @param method_name A 1-length string. The name of the method from which is desired
#'   to extract information.
#' @param truth_name A 1-length string. The name of the ground-truth.
#' @param hom_length_intervals A vector of integers. The minimum for each interval of
#'   homopolymer length. Each interval `i` ranges from `hom_length_intervals[i]` to
#'   `hom_length_intervals[i+1]`, except the last interval which upper limit is `Inf`.
#' @param interval_names A vector of characters with the same length of
#'   `hom_length_intervals`. The name for each interval of homopolymer length.
#' @param to_calculate A 1-length string. Possible values are "rates" or "pre_rec_f1".
#'   If "rates", the functin calculates the rates of TPs, FNs and FPs. If "pre_rec_f1",
#'   it calculates the precision, the recall and the F1-score.
#' @param output_method_name A 1-length string. The label of the method specified in
#'   `method_name` to be output.
#'
#' @return A 2-length list (`class_counts` and `dat_text`).
#' 
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' 
#' @export
make_homopolymer_table_to_plot <- function(input_hom_table, variant_type,
                                           method_name, truth_name,
                                           hom_length_intervals, interval_names,
                                           to_calculate, output_method_name){
  
  k <- to_calculate %in% c("rates", "pre_rec_f1")
  if(!k){
    stop("'to_calculate' argument must be either 'rates' or 'pre_rec_f1'.")
  }
  
  input_hom_table <- input_hom_table[input_hom_table$vt==variant_type,]
  
  k <- findInterval(input_hom_table$homopolymer_length_indel,
                    hom_length_intervals)
  stopifnot( !any(k==0) )
  k <- interval_names[k]
  k <- factor(k, levels=interval_names, ordered=TRUE)
  input_hom_table$homopolymer_length_intervals <- k
  
  method_class_name <- paste0(method_name, "_classification")
  k <- input_hom_table[,method_class_name] %in% c("FN", "TP", "FP")
  input_hom_table <- input_hom_table[k,]
  
  input_hom_table_split <- split(input_hom_table,
                                 input_hom_table$homopolymer_length_intervals)
  k <- lapply(input_hom_table_split, function(x){
    y <- calc_accuracy_measures(x, method_name, truth_name)
    c(y, total_count=nrow(x))
  })
  k <- do.call(rbind, k)
  k <- data.frame(k)
  k <- rownames_to_column(k, "homopolymer_length_intervals")
  k$homopolymer_length_intervals <- factor(k$homopolymer_length_intervals,
                                           levels=interval_names, ordered=TRUE)
  class_counts <- gather(k, "Measures", "percent",
                         .data$precision:.data$f1Score, factor_key=TRUE)
  
  k <- split(class_counts$total_count, class_counts$homopolymer_length_intervals)
  interval_counts <- sapply(k, unique)
  stopifnot( is.atomic(interval_counts) )
  dat_text <- data.frame(
    label=paste0("n=", interval_counts),
    x=names(interval_counts),
    y=1.05*max(class_counts$percent),
    Measures=NA
  )
  
  class_counts$method <- dat_text$method <- output_method_name
  class_counts$variant_type <- dat_text$variant_type <- variant_type
  
  res <- list(class_counts=class_counts, dat_text=dat_text)
  res
}







#' Generate data for splice-junction analysis using multiple master tables
#' 
#' Generate data to plot variant performance of sites near to and far from splice junctions
#'   comparisons using multiple master table as facets.
#'
#' @param ... Data.frames. Each data.frame is a master table.
#' @param experiment_names A vector of strings with length equal to the number of master
#'   tables input in `...`. Names of the datasets for each data table.
#' @param truth_names A vector of strings. Names of the ground truth in each data table
#'   especified in `...`.
#' @param method_dataset_name A vector of strings. Name of the datasets used to call
#'   variants, by the methods to be compared, in each input master table.
#' @param method_names A vector of strings. The name of the methods to be compared defined in
#'   the input master tables. All master table must have the same `method_names`.
#' @param output_method_names A vector of strings. The output name of the methods to be
#'   compared defined in the input master tables. The output name of the methods will be
#'   the same for all master table.
#' @param variant_type A 1-length string. Possible values are "snp" or "indel".
#' @param min_isoseq_coverage A 1-length integer. The threshod value for the minimum
#'   Iso-Seq read coverage.
#'
#' @return A list of two elements (`acc_sj` and `n_test`) to be used to draw a chart to
#'  analyse the variant-calling performance of sites near to and far from splice junctions.
#'  `acc_sj` stores the calculated performance measures, and `n_test` stores the number of
#'  observed variants in each case.
#' 
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr recode
#' @importFrom rlang .data
#' 
#' @export
splice_junction_analysis_table <- function(..., experiment_names, truth_names, method_dataset_name,
                                           method_names, output_method_names, variant_type,
                                           min_isoseq_coverage) {
  
  data_tables <- list(...)
  
  if( !all( sapply(data_tables, class) == "data.frame" ) ){
    stop("All objects in `...` must be data.frames. Make sure it's a master table.")
  }
  
  if( !(variant_type %in% c("snp", "indel")) ){
    stop("`variant_type` argument must be either 'snp' or 'indel'.")
  }else{
    variant_type <- ifelse(variant_type=="snp", 0, 1)
  }
  
  acc_n_experiments <- mapply(function(data_tables_i, experiment_names_i, truth_names_i, method_dataset_name_i){
    ### filter master table by iso-seq read coverage
    dataset_coverage_column <- paste0(method_dataset_name_i, "_coverage")
    k <- data_tables_i[,dataset_coverage_column] >= min_isoseq_coverage
    data_tables_i <- data_tables_i[k,]
    
    ### take variants near and far from splice junctions
    ### if near, many reads must contain the splice junction (at least 50% of them)
    ### if far, no one read that contain a splice junction
    ### if we do not consider n-cigar reads, percent_ss may be higher than 1.
    ###   but this is ok, because we remove variants that overlap n-cigar reads
    k <- data_tables_i$ss_highest_num / data_tables_i[,dataset_coverage_column]
    data_tables_i <- mutate(data_tables_i, "percent_ss" = k)
    data_tables_i <- filter(data_tables_i, .data$percent_ss>=.5 | .data$is_near_ss==0)
    
    ### remove variants that overlap with any intronic region
    dataset_ncrNum_column <- paste0(method_dataset_name_i, "_ncr_num")
    k <- data_tables_i[,dataset_ncrNum_column] == 0
    data_tables_i <- data_tables_i[k,]
    
    ### get the distance of the most frequent variant
    data_tables_i <- mutate(data_tables_i, ss_dist_of_the_most_freq={
      mapply(function(num, dis){
        if( any(is.na(dis)) ){
          stopifnot( length(dis)==1 )
          NA
        }else{
          dis[ which.max(num) ]
        }
      }, .data$ss_num, .data$ss_dist)
    })
    
    ### get the is_acceptor_site of the most frequent variant
    data_tables_i <- mutate(data_tables_i, is_acceptor_site_of_the_most_freq={
      mapply(function(num, acc){
        if( any(is.na(acc)) ){
          stopifnot( length(acc)==1 )
          NA
        }else{
          acc[ which.max(num) ]
        }
      }, .data$ss_num, .data$is_acceptor_site)
    })
    
    ### for each method to compare, calculate performance measures of
    ### variant calling from sites near to and far from splice junctions
    data_tables_i_ori <- data_tables_i
    env <- environment()
    env[["n_methods"]] <- NULL
    acc_methods <- mapply(function(method_names_i, output_method_names_i){
      k <- paste0("is_indel_", method_names_i)
      k <- data_tables_i_ori[,k] == variant_type
      k <- which(k)
      data_tables_i <- data_tables_i_ori[k,]
      
      env[["n_methods"]] <- c( env[["n_methods"]], paste0("n=", table(data_tables_i$is_near_ss)) )
      
      k <- split(data_tables_i, data_tables_i$is_near_ss)
      k <- sapply(k, calc_accuracy_measures, method_name=method_names_i, truth_name=truth_names_i)
      k <- as.table(k)
      acc <- as.data.frame(k)
      names(acc) <- c("Measures", "is_near", "Score")
      cbind(acc, Method=output_method_names_i)
    }, method_names, output_method_names, SIMPLIFY=FALSE)
    acc_methods <- do.call(rbind, acc_methods)
    
    acc_methods$is_near <- recode(acc_methods$is_near, "0"="No", "1"="Yes")
    acc_methods$Measures <- recode(acc_methods$Measures, "precision"="Precision",
                                          "sensitivity"="Recall", "f1Score"="F1-score")
    acc_methods$Method <- factor(acc_methods$Method, levels=output_method_names, ordered=TRUE)
    acc_methods$experiment <- experiment_names_i
    
    list(acc_methods=acc_methods, n_methods=n_methods)
  }, data_tables, experiment_names, truth_names, method_dataset_name, SIMPLIFY=FALSE)
  
  ### table of performace measures
  k <- lapply(acc_n_experiments, function(acc_n_experiments_i){
    acc_n_experiments_i$acc_methods
  })
  k <- unname(k)
  acc_sj <- do.call(rbind, k)
  rownames(acc_sj) <- NULL
  
  ### table of text annotation. usefull to inclute text into future charts to create
  k <- lapply(acc_n_experiments, function(acc_n_experiments_i){
    acc_n_experiments_i$n_methods
  })
  k <- unname(k)
  k <- do.call(c, k)
  n_text <- data.frame(
    label=k,
    experiment=rep(
      experiment_names,
      rep(length(method_names)*2,
          length(experiment_names))
    ),
    Method=gl(n=length(method_names),
              k=2,
              length=length(method_names)*2*length(experiment_names),
              labels=output_method_names,
              ordered=TRUE
    ),
    x=rep(
      c("No", "Yes"),
      length(method_names)*length(experiment_names)
    ),
    y=1.10*max(acc_sj$Score),
    Measures=NA
  )
  
  list(acc_sj=acc_sj, n_text=n_text)
}
