#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 14:36:56 2022

@author: joannafernandez
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

#Derine a function to make dataframes and normalise

def Collect_data(file1,file2,fact):
    forward = pd.read_csv(file1)
    reverse = pd.read_csv(file2)
    forward['score'] = forward['count'] + reverse['count']
    print(forward)
    print (forward['score'].sum())
    forward['norm']= (forward['score'] / forward['score'].sum())*fact
    print(forward)

    chr1= forward[(forward["chro"] == "chr1")]
#    chr1.plot.line('pos','norm')
    chr2= forward[(forward["chro"] == "chr2")]
    chr3= forward[(forward["chro"] == "chr3")]
    return chr1,chr2,chr3, forward

wtchr1, wtchr2, wtchr3, wt = Collect_data("2-RNApol2-ChIP-RZ316-w-thia.cen.f-w50.count.csv",
                                      '2-RNApol2-ChIP-RZ316-w-thia.cen.r-w50.count.csv', 1000000)

sen1chr1,sen1chr2,sen1chr3, sen1 = Collect_data("4-RNApol2-ChIP-RZ317-w-thia.cen.f-w50.count.csv",
                                          '4-RNApol2-ChIP-RZ317-w-thia.cen.r-w50.count.csv', 1000000)

dbl8chr1,dbl8chr2,dbl8chr3, dbl8 = Collect_data("6-RNApol2-ChIP-RZ318-w-thia.cen.f-w50.count.csv",
                                          '6-RNApol2-ChIP-RZ318-w-thia.cen.r-w50.count.csv', 1000000)

dschr1,dschr2,dschr3, ds = Collect_data("8-RNApol2-ChIP-RZ319-w-thia.cen.f-w50.count.csv",
                                    '8-RNApol2-ChIP-RZ319-w-thia.cen.f-w50.count.csv', 1000000)
#nmt




wtchr1, wtchr2, wtchr3, wt = Collect_data("1-RNApol2-ChIP-RZ316.cen.f-w50.count.csv",
                                      '1-RNApol2-ChIP-RZ316.cen.r-w50.count.csv', 1000000)

sen1chr1,sen1chr2,sen1chr3, sen1 = Collect_data("3-RNApol2-ChIP-RZ317.cen.f-w50.count.csv",
                                          '3-RNApol2-ChIP-RZ317.cen.r-w50.count.csv', 1000000)

dbl8chr1,dbl8chr2,dbl8chr3, dbl8 = Collect_data("5-RNApol2-ChIP-RZ318.cen.f-w50.count.csv",
                                          '5-RNApol2-ChIP-RZ318.cen.r-w50.count.csv', 1000000)

dschr1,dschr2,dschr3, ds = Collect_data("7-RNApol2-ChIP-RZ319.cen.f-w50.count.csv",
                                    '7-RNApol2-ChIP-RZ319.cen.r-w50.count.csv', 1000000)

def Find(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]


    genesfor['Sbinpos'] = genesfor['start']/50
    genesfor['Sbinpos'] = genesfor['Sbinpos'].astype(int)
    genesfor['Sbinpos'] = genesfor['Sbinpos']*50 +25
    genesfor['Ebinpos'] = genesfor['end']/50
    genesfor['Ebinpos'] = genesfor['Ebinpos'].astype(int)
    genesfor['Ebinpos'] = genesfor['Ebinpos']*50 +25


    genesrev['Sbinposr'] = genesrev['end']/50
    genesrev['Sbinposr'] = genesrev['Sbinposr'].astype(int)
    genesrev['Sbinposr'] = genesrev['Sbinposr']*50 +25
    genesrev['Ebinposr'] = genesrev['start']/50
    genesrev['Ebinposr'] = genesrev['Ebinposr'].astype(int)
    genesrev['Ebinposr'] = genesrev['Ebinposr']*50 +25

    return genesfor, genesrev, genes

genesfor, genesrev, ggenes = Find("dbl8_stall_sites_direction.txt")
controlfor, controlrev, ccontrol = Find('control_genes.txt')

sen1stall = ggenes[(ggenes['genotype'] == 'sen1D')]
dbl8stall = ggenes[(ggenes['genotype'] == 'dbl8D')]
doublestall = ggenes[(ggenes['genotype'] == 'sen1dbl8DD_unique')]

sen1stalls = sen1stall['ID'].to_list()
dbl8stalls = dbl8stall['ID'].to_list()
dsstalls = doublestall['ID'].to_list()

def Findfeat(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genes.loc[genes['Chromosome'] == "I", 'chro'] = 'chr1'
    genes.loc[genes['Chromosome'] == "II", 'chro'] = 'chr2'
    genes.loc[genes['Chromosome'] == "III", 'chro'] = 'chr3'
    
    featfor = genes[(genes['Strand'] == 'forward')]
    featrev = genes[(genes['Strand'] == 'reverse')]
    
    featrev.loc[featrev['Start position'] < 25, 'Start position'] = 25
     
    featfor['Sbinpos'] = featfor['Start position']/50
    featfor['Sbinpos'] = featfor['Sbinpos'].astype(int)
    featfor['Sbinpos'] = featfor['Sbinpos']*50 + 25
    featfor['Ebinpos'] = featfor['End position']/50
    featfor['Ebinpos'] = featfor['Ebinpos'].astype(int)
    featfor['Ebinpos'] = featfor['Ebinpos']*50 +25
    
    featrev['Sbinpos'] = featrev['End position']/50 
    featrev['Sbinpos'] = featrev['Sbinpos'].astype(int)
    featrev['Sbinpos'] = featrev['Sbinpos']*50 +25
    featrev['Ebinpos'] = featrev['Start position']/50
    featrev['Ebinpos'] = featrev['Ebinpos'].astype(int)
    featrev['Ebinpos'] = featrev['Ebinpos']*50 +25
    
    return featfor, featrev, genes

featfor, featrev, ffeat = Findfeat('protein_coding_gene_list.tsv')

feat1 = ffeat[(ffeat['chro'] == 'chr1')]
feat2 = ffeat[(ffeat['chro'] == 'chr2')]
feat3 = ffeat[(ffeat['chro'] == 'chr3')]



gene1 = ggenes[(ggenes['chro'] == 'chr1')]
gene2 = ggenes[(ggenes['chro'] == 'chr2')]
gene3 = ggenes[(ggenes['chro'] == 'chr3')]

#%%

selectionf = featfor.sample(200)
selectionr = featrev.sample(200)

#highsecf = controlfor.sample(9)
#highsecr = controlrev.sample(9)

#%%



def RStart(genesfor, chr1, chr2, chr3, p, k):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[7] == k:
            if g[8] == p:
                if g[5] == 'chr1':
                    x = chr1.index[chr1['pos'] == g[9]].tolist()
                    xx.append(x) 

                if g[5] == 'chr2':
                    x2 = chr2.index[chr2['pos'] == g[9]].tolist()
                    xx.append(x2)

                if g[5] == 'chr3':
                    x3 = chr3.index[chr3['pos'] == g[9]].tolist()
                    xx.append(x3)
                
    return xx
dsforxxR = RStart(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsrevxxR = RStart(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsforxxL = RStart(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')
dsrevxxL = RStart(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')


def REnd(genesfor, chr1, chr2, chr3, p, k):
    xxe =[]
    for ge in genesfor.itertuples():
        if ge[7] == k:
            if ge[8] == p:
                if ge[5] == 'chr1':
                    xe = chr1.index[chr1['pos'] == ge[10]].tolist()
                    xxe.append(xe) 

                if ge[5] == 'chr2':
                    xe2 = chr2.index[chr2['pos'] == ge[10]].tolist()
                    xxe.append(xe2)

                if ge[5] == 'chr3':
                    xe3 = chr3.index[chr3['pos'] == ge[10]].tolist()
                    xxe.append(xe3)
    return xxe

dsforxxeR = REnd(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsrevxxeR = REnd(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward')
dsforxxeL = REnd(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')
dsrevxxeL = REnd(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward')

#lets make two list, one is ds hh start 
# rev right and 
# for left 




#for wt , forward coding strand genes only 
#TSS
def Start(genesfor, chr1, chr2, chr3, p):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[8] == p:
            if g[5] == 'chr1':
                x = chr1.index[chr1['pos'] == g[9]].tolist()
                xx.append(x) 

            if g[5] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[9]].tolist()
                xx.append(x2)

            if g[5] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[9]].tolist()
                xx.append(x3)
                
    return xx

senforxx = Start(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
senrevxx = Start(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxx = Start(genesfor, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
dbl8revxx = Start(genesrev, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
#dsforxx = Start(genesfor, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')
#dsrevxx = Start(genesrev, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')


def ControlStart(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[2] == 'chr1':
                x = chr1.index[chr1['pos'] == g[8]].tolist()
                xx.append(x) 

            if g[2] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[8]].tolist()
                xx.append(x2)

            if g[2] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[8]].tolist()
                xx.append(x3)
                
    return xx
controlforxx = ControlStart(controlfor, wtchr1, wtchr2, wtchr3)
controlrevxx = ControlStart(controlrev, wtchr1, wtchr2, wtchr3)


def ControlEnd(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[2] == 'chr1':
                x = chr1.index[chr1['pos'] == g[9]].tolist()
                xx.append(x) 

            if g[2] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[9]].tolist()
                xx.append(x2)

            if g[2] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[9]].tolist()
                xx.append(x3)
                
    return xx
controlforxxe = ControlEnd(controlfor, wtchr1, wtchr2, wtchr3)
controlrevxxe = ControlEnd(controlrev, wtchr1, wtchr2, wtchr3)

#TES
def End(genesfor, chr1, chr2, chr3, p):
    xxe =[]
    for ge in genesfor.itertuples():
        if ge[8] == p:
            if ge[5] == 'chr1':
                xe = chr1.index[chr1['pos'] == ge[10]].tolist()
                xxe.append(xe) 

            if ge[5] == 'chr2':
                xe2 = chr2.index[chr2['pos'] == ge[10]].tolist()
                xxe.append(xe2)

            if ge[5] == 'chr3':
                xe3 = chr3.index[chr3['pos'] == ge[10]].tolist()
                xxe.append(xe3)
    return xxe

senforxxe = End(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
senrevxxe = End(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxxe = End(genesfor, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
dbl8revxxe = End(genesrev, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
#dsforxxe = End(genesfor, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')
#dsrevxxe = End(genesrev, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')


#selection of random 200 genes 

#####PLEASE NOTE RIGHT NOW ITS JUST ALL GENES
def selectionlStart(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[9] == 'chr1':
                x = chr1.index[chr1['pos'] == g[10]].tolist()
                xx.append(x) 

            if g[9] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[10]].tolist()
                xx.append(x2)

            if g[9] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[10]].tolist()
                xx.append(x3)
                
    return xx
selforxx = selectionlStart(featfor, wtchr1, wtchr2, wtchr3)
selrevxx = selectionlStart(featrev, wtchr1, wtchr2, wtchr3)


def selectionEnd(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[9] == 'chr1':
                x = chr1.index[chr1['pos'] == g[11]].tolist()
                xx.append(x) 

            if g[9] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[11]].tolist()
                xx.append(x2)

            if g[9] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[11]].tolist()
                xx.append(x3)
                
    return xx
selforxxe = selectionEnd(featfor, wtchr1, wtchr2, wtchr3)
selrevxxe = selectionEnd(featrev, wtchr1, wtchr2, wtchr3)


def Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]:ze[0],'norm'].tolist()
        yucky.append(yy)
    stalls = pd.DataFrame(yucky)   

    return stalls

sen1stallswt = Gene_bits(senforxx, senforxxe, wt)
sen1stallssen1 = Gene_bits(senforxx, senforxxe, sen1)
sen1stallsdbl8 = Gene_bits(senforxx, senforxxe, dbl8)
sen1stallsds = Gene_bits(senforxx, senforxxe, ds)

dbl8stallswt = Gene_bits(dbl8forxx, dbl8forxxe, wt)
dbl8stallssen1 = Gene_bits(dbl8forxx, dbl8forxxe, sen1)
dbl8stallsdbl8 = Gene_bits(dbl8forxx, dbl8forxxe, dbl8)
dbl8stallsds = Gene_bits(dbl8forxx, dbl8forxxe, ds)

#for left double
dsstallswtL = Gene_bits(dsforxxL, dsforxxeL, wt)
dsstallssen1L = Gene_bits(dsforxxL, dsforxxeL, sen1)
dsstallsdbl8L = Gene_bits(dsforxxL, dsforxxeL, dbl8)
dsstallsdsL = Gene_bits(dsforxxL, dsforxxeL, ds)

#for right double 
dsstallswtR = Gene_bits(dsforxxR, dsforxxeR, wt)
dsstallssen1R = Gene_bits(dsforxxR, dsforxxeR, sen1)
dsstallsdbl8R = Gene_bits(dsforxxR, dsforxxeR, dbl8)
dsstallsdsR = Gene_bits(dsforxxR, dsforxxeR, ds)

#for control 
constallswt = Gene_bits(controlforxx, controlforxxe, wt)
constallssen1 = Gene_bits(controlforxx, controlforxxe, sen1)
constallsdbl8 = Gene_bits(controlforxx, controlforxxe, dbl8)
constallsds = Gene_bits(controlforxx, controlforxxe, ds)

#for random selection
selbodywt = Gene_bits(selforxx, selforxxe, wt)
selbodysen1 = Gene_bits(selforxx, selforxxe, sen1)
selbodydbl8 = Gene_bits(selforxx, selforxxe, dbl8)
selbodyds = Gene_bits(selforxx, selforxxe, ds)



def Rev_gene_bits (list1, list2, frame):
    alldatar =[]
    for y, ye  in zip(list1, list2):
        yyr = frame.loc[ye[0]:y[0],'norm'].tolist()
        alldatar.append(yyr)
    stalls = pd.DataFrame(alldatar)
    return stalls

sen1stallsrwt = Rev_gene_bits(senrevxx, senrevxxe, wt)
sen1stallsrsen1 = Rev_gene_bits(senrevxx, senrevxxe, sen1)
sen1stallsrdbl8 = Rev_gene_bits(senrevxx, senrevxxe, dbl8)
sen1stallsrds = Rev_gene_bits(senrevxx, senrevxxe, ds)

dbl8stallsrwt = Rev_gene_bits(dbl8revxx, dbl8revxxe, wt)
dbl8stallsrsen1 = Rev_gene_bits(dbl8revxx, dbl8revxxe, sen1)
dbl8stallsrdbl8 = Rev_gene_bits(dbl8revxx, dbl8revxxe, dbl8)
dbl8stallsrds = Rev_gene_bits(dbl8revxx, dbl8revxxe, ds)

#rev right double
dsstallsrwtR = Rev_gene_bits(dsrevxxR, dsrevxxeR, wt)
dsstallsrsen1R = Rev_gene_bits(dsrevxxR, dsrevxxeR, sen1)
dsstallsrdbl8R = Rev_gene_bits(dsrevxxR, dsrevxxeR, dbl8)
dsstallsrdsR = Rev_gene_bits(dsrevxxR, dsrevxxeR, ds)

#rev left double 
dsstallsrwtL = Rev_gene_bits(dsrevxxL, dsrevxxeL, wt)
dsstallsrsen1L = Rev_gene_bits(dsrevxxL, dsrevxxeL, sen1)
dsstallsrdbl8L = Rev_gene_bits(dsrevxxL, dsrevxxeL, dbl8)
dsstallsrdsL = Rev_gene_bits(dsrevxxL, dsrevxxeL, ds)


#rev control 
constallsrwt = Rev_gene_bits(controlrevxx, controlrevxxe, wt)
constallsrsen1 = Rev_gene_bits(controlrevxx, controlrevxxe, sen1)
constallsrdbl8 = Rev_gene_bits(controlrevxx, controlrevxxe, dbl8)
constallsrds = Rev_gene_bits(controlrevxx, controlrevxxe, ds)

#rev random


selbodyrwt = Rev_gene_bits(selrevxx, selrevxxe, wt)
selbodyrsen1 = Rev_gene_bits(selrevxx, selrevxxe, sen1)
selbodyrdbl8 = Rev_gene_bits(selrevxx, selrevxxe, dbl8)
selbodyrds = Rev_gene_bits(selrevxx, selrevxxe, ds)

#%%

def Expand_genes(stall_df):
    alldata = []
    for row in stall_df.iterrows():
        xx = row[1]
        xx = xx.dropna()
        length = len(xx)
        num = 1000/length
        print (num*length)

        te = np.arange(0,1000,num, dtype = int)

        expand = np.zeros(1000)

        for itter, val in enumerate(te):
            if itter<=length-2:
                expand[val:te[itter+1]] = xx[itter]
            else: continue
        expand[te[length-1]:] = xx[length-1]
        alldata.append(expand)
    return alldata

#control
cw = Expand_genes(constallswt)
cs = Expand_genes(constallssen1)
cb = Expand_genes(constallsdbl8)
cds = Expand_genes(constallsds)

cwr = Expand_genes(constallsrwt)
csr = Expand_genes(constallsrsen1)
cbr = Expand_genes(constallsrdbl8)
cdsr = Expand_genes(constallsrds)



#sen1
aw = Expand_genes(sen1stallswt)
ass = Expand_genes(sen1stallssen1)
ad = Expand_genes(sen1stallsdbl8)
ads = Expand_genes(sen1stallsds)

raw= Expand_genes(sen1stallsrwt)
ras = Expand_genes(sen1stallsrsen1)
rad = Expand_genes(sen1stallsrdbl8)
rads = Expand_genes(sen1stallsrds)



#dbl8
dbw = Expand_genes(dbl8stallswt)
dbss = Expand_genes(dbl8stallssen1)
dbd = Expand_genes(dbl8stallsdbl8)
dbds = Expand_genes(dbl8stallsds)

rdbw= Expand_genes(dbl8stallsrwt)
rdbs = Expand_genes(dbl8stallsrsen1)
rdbd = Expand_genes(dbl8stallsrdbl8)
rdbds = Expand_genes(dbl8stallsrds)


#ds right 
dswR = Expand_genes(dsstallswtR)
dssR = Expand_genes(dsstallssen1R)
dsbR = Expand_genes(dsstallsdbl8R)
dsdsR = Expand_genes(dsstallsdsR)

rdswR = Expand_genes(dsstallsrwtR)
rdssR = Expand_genes(dsstallsrsen1R)
rdsbR = Expand_genes(dsstallsrdbl8R)
rdsdsR = Expand_genes(dsstallsrdsR)

#ds left
dswL = Expand_genes(dsstallswtL)
dssL = Expand_genes(dsstallssen1L)
dsbL = Expand_genes(dsstallsdbl8L)
dsdsL = Expand_genes(dsstallsdsL)

rdswL = Expand_genes(dsstallsrwtL)
rdssL = Expand_genes(dsstallsrsen1L)
rdsbL = Expand_genes(dsstallsrdbl8L)
rdsdsL = Expand_genes(dsstallsrdsL)

#selection 
sw = Expand_genes(selbodywt)
ss = Expand_genes(selbodysen1)
sb = Expand_genes(selbodydbl8)
sds = Expand_genes(selbodyds)

rsw = Expand_genes(selbodyrwt)
rss = Expand_genes(selbodyrsen1)
rsb = Expand_genes(selbodyrdbl8)
rsds = Expand_genes(selbodyrds)

#%%


def concati(rev, forward):
    forwards = np.stack(forward, axis=0 )
    revs = np.stack(rev, axis=0 )
    rev1 = revs[:, ::-1]
    
    new = np.concatenate([forwards,rev1])
    return new

#reverse LEFT, forward right 
HTw = concati(rdswL, dswR)
HTs = concati(rdssL, dssR)
HTb = concati(rdsbL, dsbR)
HTds = concati (rdsdsL,dsdsR)

#reverse right + forward leeft 
HHw = concati(rdswR, dswL)
HHs = concati(rdssR, dssL)
HHb = concati(rdsbR, dsbL)
HHds = concati(rdsdsR, dsdsL)

sa = concati(raw, aw)
sas = concati(ras,ass)
sab = concati(rad,ad)
sads = concati(rads,ads)

dbl8w = concati(rdbw, dbw)
dbl8s = concati(rdbs, dbss)
dbl8b = concati(rdbd, dbd)
dbl8ds = concati(rdbds, dbds)

caw = concati(cwr, cw)
cas = concati(csr, cs)
cab = concati(cbr, cb)
cads = concati(cdsr, cds)

sew = concati(rsw, sw)
ses = concati(rss, ss)
seb = concati(rsb, sb)
seds = concati(rsds, sds)

#%%

import scipy.stats as stats
from scipy.stats import t

x = pd.DataFrame(sa)
m = x[0].mean()

s = x[0].std() 
dof = len(x[0])-1 
confidence = 0.95
t_crit = np.abs(t.ppf((1-confidence)/2,dof))
print(t_crit)
(m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x))) 


def conint(array):
    confidence = 0.95
    trythis = []
    x = pd.DataFrame(array)
    #print(x)
    for column in x:

        m = (x[column]).mean()
        s = x[column].std()

        dof = len(x[column])-1 
        t_crit = np.abs(t.ppf((1-confidence)/2,dof))
        interval = (m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x)))
        
        trythis.append(interval)
        saint = pd.DataFrame(trythis, columns=['lower', 'upper'])
        #oft = saint.T
    return saint

saCON = conint(sa)
sasCON = conint(sas)
sabCON = conint(sab)
sadsCON = conint(sads)

bwCON = conint(dbl8w)
bsCON = conint(dbl8s)
bbCON = conint(dbl8b)
bdsCON = conint(dbl8ds)

HHwCON = conint(HHw)
HHsCON = conint(HHs)
HHbCON = conint(HHb)
HHdsCON = conint(HHds)


HTwCON = conint(HTw)
HTsCON = conint(HTs)
HTbCON = conint(HTb)
HTdsCON = conint(HTds)

cawCON = conint(caw)
casCON = conint(cas)
cabCON = conint(cab)
cadsCON = conint(cads)

sewCON = conint(sew)
sesCON = conint(ses)
sebCON = conint(seb)
sedsCON = conint(seds)


        
    
    
def pileup (thing):
    #thing = -thing
    #print(thing)
    #print(-thing)
    whack = np.stack( thing, axis=0 )
    stuff = whack.mean(axis=0)
    wow = pd.DataFrame(stuff)
    
    return wow
#sen compiled


#hh ht 
HHwline = pileup(HHw)
HHsline = pileup(HHs)
HHbline = pileup(HHb)
HHdsline = pileup(HHds)
HTwline = pileup(HTw)
HTsline = pileup(HTs)
HTbline = pileup(HTb)
HTdsline = pileup(HTds)

#selection
swline = pileup(sew)
ssline = pileup(ses)
sbline = pileup(seb)
sdsline = pileup(seds)


#control gene body
cwline = pileup(caw)
csline = pileup(cas)
cbline = pileup(cab)
cdsline = pileup(cads)

rcwline = pileup(cwr)
rcsline = pileup(csr)
rcbline = pileup(cbr)
rcdsline = pileup(cdsr)


#sen
awline = pileup(sa)
assline = pileup(sas)
adline = pileup(sab)
adsline= pileup(sads)


#dbl8
dbwline = pileup(dbl8w)
dbssline = pileup(dbl8s)
dbdline = pileup(dbl8b)
dbdsline= pileup(dbl8ds)

#%%

figgy, (ax1) =plt.subplots(1, sharey=True)

x = np.arange(0, 1000, 1)
ax1.plot(cwline.index, cwline[0], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(csline.index, csline[0], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(cbline.index, cbline[0], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(cdsline.index, cdsline[0], color = 'tomato', alpha=0.5, linewidth=1)
y11 = np.array(cawCON['lower']).flatten()
y10= np.array(cawCON['upper']).flatten()
ax1.fill_between(x, y10, y11, facecolor='grey', alpha =0.3,)
y11s = np.array(casCON['lower']).flatten()
y10s= np.array(casCON['upper']).flatten()
ax1.fill_between(x, y10s, y11s, facecolor='mediumaquamarine', alpha =0.3)
y11b = np.array(cabCON['lower']).flatten()
y10b = np.array(cabCON['upper']).flatten()
ax1.fill_between(x, y10b, y11b, facecolor='mediumaquamarine', alpha =0.3)
y11ds = np.array(cadsCON['lower']).flatten()
y10ds= np.array(cadsCON['upper']).flatten()
ax1.fill_between(x, y10ds, y11ds, facecolor='mediumaquamarine', alpha =0.3)
ax1.legend(loc='best')



figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(awline.index, awline[0], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(assline.index, assline[0], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(adline.index, adline[0], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(adsline.index, adsline[0], color = 'tomato', alpha=0.8, linewidth=1)
x = np.arange(0, 1000, 1)
y1 = np.array(saCON['lower']).flatten()
y2= np.array(saCON['upper']).flatten()
ax1.fill_between(x, y1, y2, facecolor='grey', alpha =0.2)
y1s =np.array(sasCON['lower']).flatten()
y2s=np.array(sasCON['upper']).flatten()
ax1.fill_between(x, y1s, y2s, facecolor='deepskyblue', alpha =0.2, )
y1b =np.array(sabCON['lower']).flatten()
y2b=np.array(sabCON['upper']).flatten()
#ax1.fill_between(x, y1b, y2b, facecolor='mediumaquamarine', alpha =0.2)
y1ds =np.array(sadsCON['lower']).flatten()
y2ds=np.array(sadsCON['upper']).flatten()
#ax1.fill_between(x, y1ds, y2ds, facecolor='mediumaquamarine', alpha =0.1)
ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
#ax1.set_ylim([-4, 124])
ax1.set_ylabel('normalised reads')
ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.fill_between(x, y1a, y2a, facecolor='mediumaquamarine', alpha =0.2)
#y2 = linerev
#ax2.fill_between(x, y2, facecolor='green', alpha =0.3)
#filkggy, (ax1) =plt.subplots(1,1, sharey=True)

#figgy, (ax1) =plt.subplots(1,1, sharey=True)
ax1.plot(dbwline.index, dbwline[0], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(dbssline.index, dbssline[0], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(dbdline.index, dbdline[0], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(dbdsline.index, dbdsline[0], color = 'tomato', alpha=0.8, linewidth=1)
y3 = np.array(bwCON['lower']).flatten()
y4 = np.array(bwCON['upper']).flatten()
ax1.fill_between(x, y3, y4, facecolor='grey', alpha =0.2)
y3s = np.array(bsCON['lower']).flatten()
y4s = np.array(bsCON['upper']).flatten()
#ax1.fill_between(x, y3s, y4s, facecolor='indianred', alpha =0.1,)
y3b = np.array(bbCON['lower']).flatten()
y4b = np.array(bbCON['upper']).flatten()
ax1.fill_between(x, y3b, y4b, facecolor='gold', alpha =0.2)
y3ds = np.array(bdsCON['lower']).flatten()
y4ds = np.array(bdsCON['upper']).flatten()
#ax1.fill_between(x, y3ds, y4ds, facecolor='indianred', alpha =0.1)
ax1.legend(loc='best')
ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
#ax1.set_ylim([-4, 124])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


#figgy, (ax1) =plt.subplots(1, sharey=True)
#ds left
ax1.plot(HHwline.index, HHwline[0], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(HHsline.index, HHsline[0], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(HHbline.index, HHbline[0], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(HHdsline.index, HHdsline[0], color = 'tomato', alpha=0.8, linewidth=1)
y5 = np.array(HHwCON['lower']).flatten()
y6= np.array(HHwCON['upper']).flatten()
ax1.fill_between(x, y5, y6, facecolor='grey', alpha =0.2)
y5s = np.array(HHsCON['lower']).flatten()
y6s= np.array(HHsCON['upper']).flatten()
#ax1.fill_between(x, y5s, y6s, facecolor='teal', alpha =0.2)
y5b = np.array(HHbCON['lower']).flatten()
y6b= np.array(HHbCON['upper']).flatten()
#ax1.fill_between(x, y5b, y6b, facecolor='teal', alpha =0.2)
y5ds = np.array(HHdsCON['lower']).flatten()
y6ds= np.array(HHdsCON['upper']).flatten()
ax1.fill_between(x, y5ds, y6ds, facecolor='indianred', alpha =0.2)
#ax1.legend(loc='best')


ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')

#figgy, (ax1) =plt.subplots(1, sharey=True)
#ds right
ax1.plot(HTwline.index, HTwline[0], color = 'black', alpha =0.8, linewidth=1)
ax1.plot(HTsline.index, HTsline[0], color = 'deepskyblue', alpha=0.8, linewidth=1)
ax1.plot(HTbline.index, HTbline[0], color = 'orange', alpha=0.8, linewidth=1)
ax1.plot(HTdsline.index, HTdsline[0], color = 'tomato', alpha=0.8, linewidth=1)
y7 = np.array(HTwCON['lower']).flatten()
y8= np.array(HTwCON['upper']).flatten()
ax1.fill_between(x, y7, y8, facecolor='grey', alpha =0.2)
y7s = np.array(HTsCON['lower']).flatten()
y8s= np.array(HTsCON['upper']).flatten()
#ax1.fill_between(x, y7s, y8s, facecolor='salmon', alpha =0.2)
y7b = np.array(HTbCON['lower']).flatten()
y8b= np.array(HTbCON['upper']).flatten()
#ax1.fill_between(x, y7b, y8b, facecolor='salmon', alpha =0.2)
y7ds = np.array(HTdsCON['lower']).flatten()
y8ds= np.array(HTdsCON['upper']).flatten()
ax1.fill_between(x, y7ds, y8ds, facecolor='salmon', alpha =0.2, )
#ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.set_ylim([-4, 124])


#figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(cwline.index, cwline[0], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(csline.index, csline[0], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(cbline.index, cbline[0], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(cdsline.index, cdsline[0], color = 'tomato', alpha=0.5, linewidth=1)
y11 = np.array(cawCON['lower']).flatten()
y10= np.array(cawCON['upper']).flatten()
ax1.fill_between(x, y10, y11, facecolor='grey', alpha =0.3,)
y11s = np.array(casCON['lower']).flatten()
y10s= np.array(casCON['upper']).flatten()
#ax1.fill_between(x, y10s, y11s, facecolor='mediumaquamarine', alpha =0.3)
y11b = np.array(cabCON['lower']).flatten()
y10b = np.array(cabCON['upper']).flatten()
#ax1.fill_between(x, y10b, y11b, facecolor='mediumaquamarine', alpha =0.3)
y11ds = np.array(cadsCON['lower']).flatten()
y10ds= np.array(cadsCON['upper']).flatten()
ax1.fill_between(x, y10ds, y11ds, facecolor='mediumaquamarine', alpha =0.3)
ax1.legend(loc='best')



ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.set_title('forward strand')


figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.plot(swline.index, swline[0], color = 'black', alpha =0.5, linewidth=1)
ax1.plot(ssline.index, ssline[0], color = 'deepskyblue', alpha=0.5, linewidth=1)
ax1.plot(sbline.index, sbline[0], color = 'orange', alpha=0.5, linewidth=1)
ax1.plot(sdsline.index, sdsline[0], color = 'tomato', alpha=0.5, linewidth=1)
y9 = np.array(sewCON['lower']).flatten()
y10= np.array(sewCON['upper']).flatten()
ax1.fill_between(x, y9, y10, facecolor='grey', alpha =0.2)



y9s = np.array(sesCON['lower']).flatten()
y10s= np.array(sesCON['upper']).flatten()
ax1.fill_between(x, y9s, y10s, facecolor='deepskyblue', alpha =0.3)
y9b = np.array(sebCON['lower']).flatten()
y10b= np.array(sebCON['upper']).flatten()
ax1.fill_between(x, y9b, y10b, facecolor='orange', alpha =0.3)
y9ds = np.array(sedsCON['lower']).flatten()
y10ds= np.array(sedsCON['upper']).flatten()
ax1.fill_between(x, y9ds, y10ds, facecolor='tomato', alpha =0.3)
ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


ax1.set_ylim([-4, 124])


#%%
import seaborn as sns
fxdsl,((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fxdsl.tight_layout()
cbar_ax = fxdsl.add_axes([.91, .3, .03, .4])
sns.heatmap(HHw, cmap = 'coolwarm', ax=ax1, vmin= 0, vmax=55, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('S∆D∆ Head-Head collisions')
ax1.set_ylabel('WT')
sns.heatmap(HHs, cmap = 'coolwarm', ax=ax3, vmin= 0, vmax=55, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(HHb, cmap = 'coolwarm', ax=ax5, vmin= 0, vmax=55, cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(HHds, cmap = 'coolwarm', ax=ax7, vmin= 0,vmax=55, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')
ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])

fuckkkkkk, (ax2, ax4, ax6, ax8) = plt.subplots(4,1)
fuckkkkkk.tight_layout()
cbar_ax = fuckkkkkk.add_axes([.91, .3, .03, .4])
sns.heatmap(HTw, cmap = 'coolwarm', ax=ax2, vmin= 0, vmax=55, cbar=True, cbar_ax=cbar_ax)
ax2.set_ylabel('WT')
ax2.axes.yaxis.set_ticks([])
ax2.set_title('S∆D∆ co-directional collisions')
ax2.axes.xaxis.set_ticks([])
sns.heatmap(HTs, cmap = 'coolwarm', ax=ax4, vmin= 0, vmax=55, cbar=True, cbar_ax=cbar_ax)
ax4.axes.yaxis.set_ticks([])
ax4.set_ylabel('Sen1∆')
ax4.axes.xaxis.set_ticks([])
sns.heatmap(HTb, cmap = 'coolwarm', ax=ax6, vmin= 0, vmax=55, cbar=True, cbar_ax=cbar_ax)
ax6.axes.yaxis.set_ticks([])
ax6.set_ylabel('Dbl8∆')  
ax6.axes.xaxis.set_ticks([])
sns.heatmap(HTds, cmap = 'coolwarm', ax=ax8, vmin= 0, vmax=55, cbar=True, cbar_ax=cbar_ax)
ax8.axes.yaxis.set_ticks([])
ax8.set_ylabel('Sen1∆Dbl8∆')

ax8.set_xticks([0,1000])
ax8.set_xticklabels(['TSS','TES'])


#CONTROL HEAT MAP
fix, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fix.tight_layout()
cbar_ax = fix.add_axes([.91, .3, .03, .4])
sns.heatmap(caw, cmap = 'coolwarm', ax=ax1,vmin= 0, vmax=110, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(cas, cmap = 'coolwarm', ax=ax3,vmin= 0, vmax=110,cbar=True, cbar_ax=cbar_ax)

ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(cab, cmap = 'coolwarm', ax=ax5,vmin= 0, vmax=110, cbar=True, cbar_ax= cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
    
sns.heatmap(cads, cmap = 'coolwarm', ax=ax7,vmin= 0, vmax=110,cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])


#DBL8 heatmap !!!!!
fxb, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fxb.tight_layout()
cbar_ax = fxb.add_axes([.91, .3, .03, .4])
sns.heatmap(dbl8w, cmap = 'coolwarm', ax=ax1, vmin= 0, vmax=160,  cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax3, vmin= 0, vmax=160,  cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax5, vmin= 0, vmax=160,  cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=ax7, vmin= 0, vmax=160,  cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])

#SEN1 heatmap bitches 
fx, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fx.tight_layout()
cbar_ax = fx.add_axes([.91, .3, .03, .4])
sns.heatmap(sa, cmap = 'coolwarm', ax=ax1, vmin= 0, vmax=120, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(sas, cmap = 'coolwarm', ax=ax3, vmin= 0, vmax=120, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(sab, cmap = 'coolwarm', ax=ax5, vmin= 0, vmax=120, cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(sads, cmap = 'coolwarm', ax=ax7, vmin= 0, vmax=120, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')


ax7.set_xticks([0,1000])
ax7.set_xticklabels(['TSS','TES'])




