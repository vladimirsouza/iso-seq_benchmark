##### In this script, we add the variant calls of a new method to a existing master
##### table.
##### To do this, first we need to generate the master table of the new method.
##### Once we get the master table of the new method, we can merge it with a
##### existing master table.



### input variables
VCF_NEW_METHOD_FILE <- "/home/vbarbo/project_2021/datasets/gloria_data/analysis/rna_ground_truth/no_splitNCigarReads/jurkat_rna_notSplit.recal_pass.vcf.gz"



### libraries
library(variantCallingFromIsoSeq)

### load a existing master table
data("vlad_master_table_v3")
dat <- vlad_master_table_v3

###
new_dat <- initiate_master_table(VCF_NEW_METHOD_FILE, method_names="shortRead_rna_gatk")
