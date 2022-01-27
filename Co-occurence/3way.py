import itertools


import pandas as pd
import os
import numpy as np
import datetime
import dask
import dask.dataframe as dd
from tqdm import tqdm

cnv = pd.read_table('~/data/Cancer_ML/cnv/cnv_merged.tsv',sep='\t',index_col=0)
mut = pd.read_table('~/data/Cancer_ML/mutation/mut_merged.tsv',sep='\t',index_col=0)
genes = (set(mut.columns) & set(cnv.columns))
samples = (set(mut.index) & set(cnv.index))
cnv = cnv.loc[samples]
mut = mut.loc[samples]
cnv = cnv[genes]
mut = mut[genes]
cnv = cnv.replace(-1,'Loss').replace(1,'Gain').replace(0,'WT').replace(-2,'Loss').replace(2,'Gain')
mut = mut.replace(1,'Mut').replace(0,'NoMut')

co=[]
for c in tqdm(co_occur.columns.tolist()):
    df_tmp = pd.DataFrame(co_occur[[c]].value_counts())
    df_tmp.columns = [c]
    co.append(df_tmp.T)
co = pd.concat(co)
co = co.fillna(0)

pools = itertools.permutations(co.index, 2)


def get_3way_interaction(p):
    without2ndmt=[]
    with2ndmt=[]
    for i in range(0,len(co_occur)):
        if (co_occur[p[1]][i].count('NoMut')==0):
            with2ndmt.append(co_occur[p[0]][i])
        elif (co_occur[p[1]][i].count('NoMut')!=0):
            without2ndmt.append(co_occur[p[0]][i])
    
    a = pd.value_counts(np.array(with2ndmt))
    a['Tissue'] = 'HNSC'
    a['Target'] = 'WithSecondM'
    b = pd.value_counts(np.array(without2ndmt))
    b['Tissue'] = 'HNSC'
    b['Target'] = 'WOSecondM' 
    column_name = ['Mut_Loss','NoMut_Loss','Mut_Gain','Mut_Loss','Mut_WT','NoMut_Loss']
    for c_name in column_name:
        try:
            a[c_name]
        except KeyError:
            a[c_name]=0
    for c_name in column_name:
        try:
            b[c_name]
        except KeyError:
            b[c_name]=0
    
    c = pd.concat([a,b],axis=1).T.fillna(0)
    c.index =[(p[0]+'_'+p[1]),(p[0]+'_'+p[1])]

    c = c[['Tissue','Target','Mut_Loss','NoMut_Loss','Mut_Gain','Mut_Loss','Mut_WT','NoMut_Loss']]
    return(c)


from multiprocessing import Process, Manager


if __name__ == "__main__":
    pool=multiprocessing.Pool(processes=60)
    manager=multiprocessing.Manager()
    sharedlist=manager.list()
    tasks = [(x,p) for x in list]
    pool.starmap(get_3way_interactiondo_stuff, tasks)
    pool.close()
    
    
L.to_csv('/home/suhan/data/Cancer_ML/3way_result.tsv',sep='\t')
        
        
