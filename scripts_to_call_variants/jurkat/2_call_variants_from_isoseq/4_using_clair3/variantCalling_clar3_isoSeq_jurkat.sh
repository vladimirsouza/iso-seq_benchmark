#!/bin/bash



###################################
###### run clair3 on iso-seq ######
###################################

conda deactivate
conda activate clair3


REF="/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta"
THREADS=10


### Clair3 (alone)
INPUT_BAM="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/aln.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/alone"
# call variants
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/pileup.vcf.gz \
  > $OUTPUT_DIR/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/pileup_pass.vcf \
  > $OUTPUT_DIR/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/pileup_pass.vcf.gz
rm $OUTPUT_DIR/pileup_pass.vcf


### SNCR + Clair3
INPUT_BAM="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate/aln_split.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/sncr"
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/pileup.vcf.gz \
  > $OUTPUT_DIR/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/pileup_pass.vcf \
  > $OUTPUT_DIR/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/pileup_pass.vcf.gz
rm $OUTPUT_DIR/pileup_pass.vcf


### SNCR + FC + Clair3
INPUT_BAM="/home/vbarbo/project_2021/datasets/gloria_data/analysis/dv_calls/noMarkDuplicate/aln_split_flagCorrection.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/gloria_data/analysis/clair3_isoSeq/sncr_fc"
/home/vbarbo/programs/Clair3/run_clair3.sh \
  --bam_fn=$INPUT_BAM \
  --ref_fn=$REF \
  --threads=$THREADS \
  --platform="hifi" \
  --model_path="/home/vbarbo/programs/Clair3/models/hifi" \
  --output=$OUTPUT_DIR
# filter only PASS
bcftools view -f PASS \
  $OUTPUT_DIR/pileup.vcf.gz \
  > $OUTPUT_DIR/pileup_pass.vcf
# compress and index
bgzip -c $OUTPUT_DIR/pileup_pass.vcf \
  > $OUTPUT_DIR/pileup_pass.vcf.gz
bcftools index $OUTPUT_DIR/pileup_pass.vcf.gz
tabix -p vcf $OUTPUT_DIR/pileup_pass.vcf.gz
rm $OUTPUT_DIR/pileup_pass.vcf



