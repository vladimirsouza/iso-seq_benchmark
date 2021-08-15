#!/bin/bash



############################
###### install clair3 ######
############################


conda deactivate

# create and activate an environment named clair3
conda create -n clair3 python=3.6.10 -y
source activate clair3

# install pypy and packages in the environemnt
conda install -c conda-forge pypy3.6 -y
pypy3 -m ensurepip
pypy3 -m pip install mpmath==1.2.1

# install python packages in environment
pip3 install tensorflow==2.2.0
pip3 install tensorflow-addons==0.11.2 tables==3.6.1
conda install -c anaconda pigz==2.4 -y
conda install -c conda-forge parallel=20191122 zstd=1.4.4 -y
conda install -c conda-forge -c bioconda samtools=1.10 -y
conda install -c conda-forge -c bioconda whatshap=1.0 -y

# clone Clair3
git clone https://github.com/HKU-BAL/Clair3.git
cd Clair3

# download pre-trained models
mkdir models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar -zxvf clair3_models.tar.gz -C ./models



###################################
###### run clair3 on iso-seq ######
###################################

conda deactivate
conda activate clair3


REF="/home/vbarbo/project_2021/datasets/reference/GRCh38.p13_genome_only_chrm/GRCh38.p13_all_chr.fasta"
THREADS=30


### Clair3 (alone)
INPUT_BAM="/home/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_s.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/wtc11/methods_to_comp/clair3/alone"
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
INPUT_BAM="/home/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_sncr.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/wtc11/methods_to_comp/clair3/sncr"
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
INPUT_BAM="/home/vbarbo/project_2021/datasets/wtc11/manipulate_data/aln_sncr_fc.bam"
OUTPUT_DIR="/home/vbarbo/project_2021/datasets/wtc11/methods_to_comp/clair3/sncr_fc"
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



