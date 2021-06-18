#' The jurkat_master_table_v5_notFiltered master table
#'
#' A master table generated from Jurkat human cells. Variant callers to compare
#' are DeepVariant and GATK on Iso-Seq data. The ground-truth was generated with
#' GATK from short-read data.
#' 
#' Methods to avaluate included: dv, dv_s, dv_s_fc, gatk_s.
#' 
#' Ground-truth: tr_dna_merged.
#'
#' @docType data
#'
#' @usage data(jurkat_master_table_v5_notFiltered)
#'
#' @format A data.frame with dimensions equal to (1020553, 35).
#'
#' Each row is a variant called by at least one of the methods. For column
#' descriptions, see \url{https://docs.google.com/document/d/1REIXRS6_oGEjTtlokwZ0pBQmoZcEE7L6Nkq-nCMoI5c/edit}
"jurkat_master_table_v5_notFiltered"




#' The jurkat_master_table_v6_notFiltered master table
#'
#' A master table generated from Jurkat human cells. Variant callers to compare
#' are DeepVariant and GATK on Iso-Seq data. The ground-truth was generated with
#' GATK from short-read data.
#' 
#' Updates from `jurkat_master_table_v6_notFiltered`:
#'   * DeepVariant calls filtered by QUAL>=20;
#'   * VCF tag QD of the gound-truth was added as column `qd_trDnaMerged`.
#' 
#' Methods to avaluate included: dv_s_fc_q20, gatk_s.
#' 
#' Ground-truth: tr_dna_merged.
#'
#' @docType data
#'
#' @usage data(jurkat_master_table_v6_notFiltered)
#'
#' @format A data.frame with dimensions equal to (1003676, 26).
#'
#' Each row is a variant called by at least one of the methods. For column
#' descriptions, see \url{https://docs.google.com/document/d/1REIXRS6_oGEjTtlokwZ0pBQmoZcEE7L6Nkq-nCMoI5c/edit}
"jurkat_master_table_v6_notFiltered"




#' The jurkat_master_table_v7_notFiltered master table
#'
#' A master table generated from Jurkat human cells. Variant callers to compare
#' are DeepVariant and GATK on Iso-Seq data. The ground-truth was generated with
#' GATK from short-read data.
#' 
#' Each row is a variant called by at least one of the methods. For column
#' descriptions, see \url{https://docs.google.com/document/d/1REIXRS6_oGEjTtlokwZ0pBQmoZcEE7L6Nkq-nCMoI5c/edit}
#' 
#' Methods to evaluate: dv_s_fc, gatk_s.
#' 
#' Ground-truth: tr_dna_merged.
#' 
#' No filtering applied.
#' 
#' @docType data
#'
#' @usage data(jurkat_master_table_v7_notFiltered)
#'
#' @format A data.frame with dimensions equal to (1018480, 32).
"jurkat_master_table_v7_notFiltered"
