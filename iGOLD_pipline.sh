#!/bin/bash
fn_DEG=$1
fn_SNP=$2
#enrich DEG and SNPs
./enrich_DEG_SNP_forline.py $fn_DEG $fn_SNP
#LDSC enrich, generate sumstat
python ~/softwares/ldsc/munge_sumstats.py \
    --sumstats ./input_sumstat.txt \
    --merge-alleles ~/database_fan/1000G/w_hm3.snplist \
    --out input_sumstat \
    --a1-inc
#calculate LDSC enrichment
python ~/softwares/ldsc/ldsc.py \
    --h2 input_sumstat.sumstats.gz \
    --intercept-h2 1 \
    --ref-ld-chr ./enrich_LDSC_lds/chr@ \
    --w-ld-chr ~/database_fan/1000G/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.\
    --overlap-annot \
    --frqfile-chr ~/database_fan/1000G/1000G_Phase3_frq/1000G.EUR.QC.\
    --out enrich_LDSC
#summarize results
./summary_result.py
#initiate
#rm input+sumstat.txt
#rm enrich_DEG_SNP.csv
#rm enrich_LDSC.result
#rm enrich_result.csv
