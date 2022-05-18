#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:19:47 2022

@author: fanc232
"""

import pandas as pd
import sys
from scipy import stats
import numpy as np

#read modules
fn_mod1='./ori_files/PD_module_gene_unsigned.map'
df_mod1=pd.read_csv(fn_mod1,sep=' ')
fn_mod2='./ori_files/PD_BR_mod_unsigned.map'
df_mod2=pd.read_csv(fn_mod2,sep=' ')
df_mod2.rename(columns={'trans':'probe'},inplace=True)
df_mod_T=pd.concat([df_mod1,df_mod2],axis=0)

#read DEGs
fn_DEG=sys.argv[1]
#fn_DEG='./enrich_DEG/test_input.xlsx' #check code!!!
df_DEG=pd.read_excel(fn_DEG,header=None)
l_DEG=set(df_DEG.iloc[:,0]) #113 DEGs in the example

#must contain information for LDSC: rs, N,Z,A1,A2,P
#must contain information for SNP: CHR,POS
fn_SNP=sys.argv[2]
#fn_SNP='./enrich_SNP/example_input_lancet.xlsx' #check code!!!
df_SNP=pd.read_excel(fn_SNP)
#process SNP files for LDSC enrichment
df_SNP1=df_SNP.iloc[:,0:6]
df_SNP1.to_csv('input_sumstat.txt',index=False,sep=' ')
#process SNP files for SNP enrichment, cutoff as 1E-5
df_SNP2=df_SNP.loc[df_SNP.P<1E-5,['CHR','LOC']]
#fine trans-regulational SNPs
df_SNP2.CHR=df_SNP2['CHR'].astype('str')
df_SNP2.CHR=[x.strip('chr') for x in df_SNP2.CHR]
df_SNP2.CHR=df_SNP2['CHR'].astype('int')
df_SNP2.LOC=df_SNP2['LOC'].astype('int')
pd_gene=pd.read_csv('./ori_files/probe_location_noXY.info',sep=' ')
pd_gene.eval('S1=START-10000',inplace=True)
pd_gene.eval('S2=END+10000',inplace=True)
def map_SNP(nrow):
    chr=pd_gene.iloc[nrow,1]
    st=pd_gene.iloc[nrow,4]
    ed=pd_gene.iloc[nrow,5]
    df_SNP3=df_SNP2[df_SNP2.CHR==chr]
    l_SNP_pos=list(df_SNP3.LOC)
    l_SNP=[x for x in l_SNP_pos if x <= ed and x >= st]
    return(len(l_SNP))
pd_gene['nSNP']=pd_gene.index.map(map_SNP) #time costing
pd_gene1=pd_gene[pd_gene.nSNP>0] #only four left in the example.
l_SNP_probe=set(pd_gene.GENE)

#enrich DEG. 
def enrich_DEG(modd:str): #do not use "mod"!
    l_mod_gene=set(df_mod_T.gene[df_mod_T.Mod==modd])
    not_mod_gene=set(df_mod_T.gene[df_mod_T.Mod!=modd])
    c1=len(l_mod_gene & l_DEG)
    c2=len(l_mod_gene - l_DEG)
    c3=len(not_mod_gene & l_DEG)
    c4=len(not_mod_gene - l_DEG)
    #print([c1,c2,c3,c4])
    res=stats.fisher_exact([[c1,c2],[c3,c4]],alternative='greater')
    p=res[1]
    return(-np.log10(p))
l_mod=list(set(df_mod_T.Mod))
l_deg_p=list(map(enrich_DEG,l_mod))
#enrich SNP. 
def enrich_SNP(modd:str): #do not use "mod"!
    l_mod_probe=set(df_mod_T.probe[df_mod_T.Mod==modd])
    not_mod_probe=set(df_mod_T.probe[df_mod_T.Mod!=modd])
    c1=len(l_mod_probe & l_SNP_probe)
    c2=len(l_mod_probe - l_SNP_probe)
    c3=len(not_mod_probe & l_SNP_probe)
    c4=len(not_mod_probe - l_SNP_probe)
    #print([c1,c2,c3,c4])
    res=stats.fisher_exact([[c1,c2],[c3,c4]],alternative='greater')
    p=res[1]
    return(-np.log10(p))
l_snp_p=list(map(enrich_SNP,l_mod))

#output
df_res=pd.DataFrame(data={'Mod':l_mod,'logP_DEG':l_deg_p,'logP_SNP':l_snp_p})
df_res.to_csv('enrich_DEG_SNP.csv',index=False)
