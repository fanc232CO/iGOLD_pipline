#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 15:29:30 2022

@author: fanc232
"""
import pandas as pd
import numpy as np

#read results of enriching DEG_SNP
df_res1=pd.read_csv('enrich_DEG_SNP.csv')
#read result of sLDSC
def read_LDSC_res(fn):
    df_res=pd.read_csv(fn,sep='\t',usecols=[0,6])
    df_res['Mod']=[str(x).split('_')[0].split('L2')[0] for x in df_res.Category]
    df_res['logP_ldsc']=df_res.Enrichment_p.apply(lambda x:-np.log10(x))
    return(df_res)
#df_res2=pd.read_csv('./enrich_LDSC.results')
df_res2=read_LDSC_res('./enrich_LDSC.results')
#check code!
#df_res2=read_LDSC_res('./ori_files/GCST007780.results')

#df_res2=pd.concat([df_res2,df_res3],axis=0)

df_res_T=pd.merge(df_res1,df_res2.loc[:,['Mod','logP_ldsc']], on='Mod')
#ldsc files should be re-generated to include both the CCM andn SCM modules.
df_res_T['logP_sum']=df_res_T.iloc[:,1:4].sum(1)
df_res_T.sort_values(by='logP_sum',inplace=True,ascending=False)
df_res_T.to_csv('enrich_result.csv',index=False)

