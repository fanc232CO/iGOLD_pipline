#!/usr/bin/python3
import pandas as pd
import sys

Mod_name=sys.argv[1]
input_DEG=sys.argv[2] #excel,differential expressed genes, columns should be gene, FC and Pval
fn_dis_gene=sys.argv[3]

fn_mod_gene='./ori_files/PD_module_gene_unsigned.map'
if 'BR' in Mod_name:
    fn_mod_gene='./ori_files/PD_BR_mod_unsigned.map'
pd_mod_gene=pd.read_csv(fn_mod_gene,sep=' ',usecols=[1,2])
mod_gene=pd_mod_gene['gene'][pd_mod_gene.Mod==Mod_name].to_list()

pd_dis_gene=pd.read_excel(fn_dis_gene,header=None)
dis_gene=pd_dis_gene.iloc[:,0].to_list()
dis_mod_gene=set(mod_gene) & set(dis_gene)
pd1=pd.DataFrame(data={'dis_mod_gene':list(dis_mod_gene)})
pd1.to_excel('disease_associated_genes_in_module.xlsx')

pd_DEG=pd.read_excel(input_DEG)
print(pd_DEG)
pd2=pd_DEG[pd_DEG['gene'].isin(mod_gene)]
pd2.to_excel('DEG_in_module.xlsx')

