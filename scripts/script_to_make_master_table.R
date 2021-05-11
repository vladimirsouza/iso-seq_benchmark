##### this script uses functions from the package to compare and validate VCF files.



### packages
library(variantCallingFromIsoSeq)



### some input variables
THREADS <- 40
MIN_DIST_FROM_SPLICE_SITE <- 20


### load the coveraga of BAM files generated from methods desired to compare
load("/home/vbarbo/load_later/dv_cover.RData")
load("/home/vbarbo/load_later/dv_nc_cover.RData")
load("/home/vbarbo/load_later/dv_nc_wh_cover.RData")
load("/home/vbarbo/load_later/dv_u_cover.RData")
load("/home/vbarbo/load_later/splice_sites.RData")
load("/home/vbarbo/load_later/tr_100_cover.RData")
load("/home/vbarbo/load_later/tr_150_cover.RData")
load("/home/vbarbo/load_later/tr_dna_merged_cover.RData")
load("/home/vbarbo/load_later/tr_rna_split_cover.RData")

### load BAM files
BAM_FILE_NOSPLIT <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/aln_rg_dedupped.bam"
bam_noSplit <- readGAlignments(BAM_FILE_NOSPLIT)



DEEPVARIANT <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noSplitBam/deepvariant_calls_pass.vcf.gz"
DEEPVARIANT_SNCR <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/deepvariant_calls_pass.vcf.gz"
DEEPVARIANT_SNCR_RC <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noClipping/deepvariant_noClipping.vcf.gz"
DEEPVARIANT_SNCR_RC_WH <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/whatshap/deepvariant2_noClipping_haplotagged.vcf"
ISOSEQ_GATK <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/gatk_calls_isoSeq/notSplit/isoSeq_jurkat.recal_pass.vcf.gz"
ISOSEQ_GATK_SNCR <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/gatk_calls_isoSeq/isoSeq_split_gatk.recal_pass.vcf.gz"
ISOSEQ_GATK_SNCR_RC <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/gatk_calls_isoSeq/split_removeClipping/isoSeq_split_removeClipping_gatk.recal_pass.vcf.gz"
ISOSEQ_GATK_SNCR_RC_WH <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/gatk_calls_isoSeq/split_removeClipping_whatsHap/isoSeq_jurkat_split_removeClipping_whatsHap_haplotagged_pass.vcf.gz"
TRUTH_DNA_100 <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_100bp/remove_intronic_regions/jurkat100bp.recal_exons.vcf.gz"
TRUTH_DNA_150 <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/my_ground_truth_jurkat_wgs_pe_150bp/remove_intronic_regions/jurkat150bp.recal_exons.vcf.gz"
TRUTH_RNA <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/rna_ground_truth/remove_intronic_regions/gloria_jurkat_shortreads.recal_exons_sampleName.vcf.gz"
TRUTH_DNA_MERGED <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/truth_merged_bams/merged.recal_exons_pass.vcf.gz"


method_names <- c("dv", "dv_s", "dv_s_rc", "dv_s_rc_wh",
                  "gatk", "gatk_s", "gatk_s_rc", "gatk_s_rc_wh",
                  "tr_100", "tr_150", "tr_rna",
                  "tr_dna_merged")
master_table <- initiate_master_table(DEEPVARIANT,
                                      DEEPVARIANT_SNCR,
                                      DEEPVARIANT_SNCR_RC,
                                      DEEPVARIANT_SNCR_RC_WH,
                                      ISOSEQ_GATK,
                                      ISOSEQ_GATK_SNCR,
                                      ISOSEQ_GATK_SNCR_RC,
                                      ISOSEQ_GATK_SNCR_RC_WH,
                                      TRUTH_DNA_100,
                                      TRUTH_DNA_150,
                                      TRUTH_RNA,
                                      TRUTH_DNA_MERGED,
                                      method_names=method_names)

