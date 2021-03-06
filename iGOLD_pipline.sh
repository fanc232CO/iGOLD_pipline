#!/bin/bash
fn_DEG=$1
fn_SNP=$2
#change by users:
run_anaconda3= ~/anaconda3/bin/python3
run_munge_sumstats= ~/softwares/ldsc/munge_sumstats.py
run_ldsc= ~/softwares/ldsc/ldsc.py
fn_snplist= ~/database_fan/1000G/w_hm3.snplist
path_weights_hm3= ~/database_fan/1000G/1000G_Phase3_weights_hm3_no_MHC
path_frq= ~/database_fan/1000G/1000G_Phase3_frq


#enrich DEG and SNPs
$run_anaconda3 ./enrich_DEG_SNP_forline.py $fn_DEG $fn_SNP
#LDSC enrich, generate sumstat
python $run_munge_sumstats \
    --sumstats ./input_sumstat.txt \
    --merge-alleles $fn_snplist \
    --out input_sumstat \
    --a1-inc
#calculate LDSC enrichment
python $run_ldsc \
    --h2 input_sumstat.sumstats.gz \
    --intercept-h2 1 \
    --ref-ld-chr ./enrich_LDSC_lds/chr@ \
    --w-ld-chr $path_weights_hm3/weights.hm3_noMHC.\
    --overlap-annot \
    --frqfile-chr $path_frq/1000G.EUR.QC.\
    --out enrich_LDSC
#summarize results
./summary_result.py
#initiate
#rm input+sumstat.txt
#rm enrich_DEG_SNP.csv
#rm enrich_LDSC.result
#rm enrich_result.csv
