#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 18:21:26 2024

@author: patricfernandez
"""
import pandas as pd
import numpy as np
from os import listdir
#import random
#import math
#from pysam import FastaFile
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from scipy import signal
import seaborn as sns
#import os
#from plotnine import *
import re


def Find(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]


    return genesfor, genesrev, genes

genesfor, genesrev, ggenes = Find("dbl8_stall_sites_direction.txt")


#trythisandcry = ggenes.merge(Genes2, on='ID', how='left')

Genes2 = pd.read_csv('Schizosaccharomyces_pombe_all_chromosomes.gff3',sep='\t', names=('1', '2','3', '4','5', '6','7', '8', '9'),skiprows=1)


for line in Genes2.itertuples():
     index = line[0]
     infor = line[9]
     position = infor.rfind(';')
     if position >0:
         Genes2.at[index, '9'] = infor[3:position]
     else:
         Genes2.at[index, '9'] = infor[3:]

Genes2.drop(labels = ['2', '6', '8'], inplace =True, axis = 1)
Genes2.rename(columns={'1':'chro', '3':'type','9':'ID', '4':'start', '5':'stop', '7':'strand'}, inplace = True)


def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    print(genes)
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]


    return genesfor, genesrev, genes


newcontrolfor, newcontrolrev, newccontrol = Find('new_control.csv')

con1 = newccontrol[(newccontrol['chro'] == 'chr1')]
con2 = newccontrol[(newccontrol['chro'] == 'chr2')]
con3 = newccontrol[(newccontrol['chro'] == 'chr3')]


sen1stall = ggenes[(ggenes['genotype'] == 'sen1D')]
dbl8stall = ggenes[(ggenes['genotype'] == 'dbl8D')]
doublestall = ggenes[(ggenes['genotype'] == 'sen1dbl8DD_unique')]


def Findfeat(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)
    genes.loc[genes['Chromosome'] == "I", 'chro'] = 'chr1'
    genes.loc[genes['Chromosome'] == "II", 'chro'] = 'chr2'
    genes.loc[genes['Chromosome'] == "III", 'chro'] = 'chr3'
    
    featfor = genes[(genes['Strand'] == 'forward')]
    featrev = genes[(genes['Strand'] == 'reverse')]

    
    return featfor, featrev, genes

featfor, featrev, ffeat = Findfeat('protein_coding_gene_list.tsv')
featfory, featrevy, ffeaty = Findfeat('protein_coding_gene_list.tsv')

feat1 = ffeat[(ffeat['chro'] == 'chr1')]
feat2 = ffeat[(ffeat['chro'] == 'chr2')]
feat3 = ffeat[(ffeat['chro'] == 'chr3')]



gene1 = ggenes[(ggenes['chro'] == 'chr1')]
gene2 = ggenes[(ggenes['chro'] == 'chr2')]
gene3 = ggenes[(ggenes['chro'] == 'chr3')]



aaaah = {'start':[3753687,1602264,1070904], 'end': [3789421,1644747,1137003], 'Chromosome': [1,2,3]}
comcentro = pd.DataFrame(data=aaaah)
comcentro['Sbinpos'] = comcentro['start']/50
comcentro['Sbinpos'] = comcentro['Sbinpos'].astype(int)
comcentro['Sbinpos'] = comcentro['Sbinpos']*50 +25
comcentro['Ebinpos'] = comcentro['end']/50
comcentro['Ebinpos'] = comcentro['Ebinpos'].astype(int)
comcentro['Ebinpos'] = comcentro['Ebinpos']*50 +25


d = {'start': [3753687], 'end': [3789421]}
centro1 = pd.DataFrame(data=d)

dd = {'start': [1602264], 'end': [1644747]}
centro2 = pd.DataFrame(data=dd)

ddd = {'start': [1070904], 'end': [1137003]}
centro3 = pd.DataFrame(data=ddd)


te1 = {'start': [1, 5554844], 'end': [29663,5579133]}
telo1 = pd.DataFrame(data=te1)
te2 ={'start': [1, 4500619], 'end': [39186,4539800]}
telo2 = pd.DataFrame(data=te2)
te3 ={'start': [], 'end': []}
telo3 = pd.DataFrame(data=te3)


#%%
#import significance files
def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes = genes.rename(columns={'Unnamed: 0': 'ID'})
    print(genes)
    

    return genes

wtvSEN = Find("wt_vs_sen_wthia.csv")
SENvDBL = Find("wt_vs_dbl_wthia.csv")
wtvDS = Find("wt_vs_ds_wthia.csv")
#DBLvDS = Find("significantDE_DSvsdbl.csv")


wtvSEN = wtvSEN.merge(Genes2, on='ID', how='left')
SENvDBL = SENvDBL.merge(Genes2, on='ID', how='left')
wtvDS = wtvDS.merge(Genes2, on='ID', how='left')
#DBLvDS = DBLvDS.merge(Genes2, on="ID", how='left')


wtvSEN = ggenes.merge(wtvSEN, on='ID', how='left')
SENvDBL = ggenes.merge(SENvDBL, on='ID', how='left')
wtvDS = ggenes.merge(wtvDS, on='ID', how='left')

wtvSEN.dropna(subset=['lfcSE'], inplace=True)
SENvDBL.dropna(subset=['lfcSE'], inplace=True)
wtvDS.dropna(subset=['lfcSE'], inplace=True)






import matplotlib.pyplot as plt
from matplotlib_venn import venn3



wtvSEN_ids = set(wtvSEN['ID'])
SENvDBL_ids = set(SENvDBL['ID'])
wtvDS_ids = set(wtvDS['ID'])

# Create a Venn diagram
venn_labels = {'100': len(wtvSEN_ids - SENvDBL_ids - wtvDS_ids),
               '010': len(SENvDBL_ids - wtvSEN_ids - wtvDS_ids),
               '001': len(wtvDS_ids - wtvSEN_ids - SENvDBL_ids),
               '110': len(wtvSEN_ids & SENvDBL_ids - wtvDS_ids),
               '101': len(wtvSEN_ids & wtvDS_ids - SENvDBL_ids),
               '011': len(SENvDBL_ids & wtvDS_ids - wtvSEN_ids),
               '111': len(wtvSEN_ids & SENvDBL_ids & wtvDS_ids)}

venn_labels_str = {key: str(val) for key, val in venn_labels.items()}

venn3(subsets=(len(wtvSEN_ids), len(SENvDBL_ids), len(wtvSEN_ids & SENvDBL_ids),
               len(wtvDS_ids), len(wtvSEN_ids & wtvDS_ids), len(SENvDBL_ids & wtvDS_ids),
               len(wtvSEN_ids & SENvDBL_ids & wtvDS_ids)),
      set_labels=('wtvSEN', 'SENvDBL', 'wtvDS'))

plt.title("Venn Diagram of IDs")
plt.show()



set_colors = ('blue', 'green', 'red')

# Create a Venn diagram with custom colors
venn3(subsets=(len(wtvSEN_ids), len(SENvDBL_ids), len(wtvSEN_ids & SENvDBL_ids),
               len(wtvDS_ids), len(wtvSEN_ids & wtvDS_ids), len(SENvDBL_ids & wtvDS_ids),
               len(wtvSEN_ids & SENvDBL_ids & wtvDS_ids)),
      set_labels=('wtvSEN', 'WTvDBL', 'wtvDS'),
      set_colors=set_colors)

plt.title("Customized Venn Diagram of IDs")
plt.show()


wtvSEN_ids = set(wtvSEN['ID'])
SENvDBL_ids = set(SENvDBL['ID'])
wtvDS_ids = set(wtvDS['ID'])

# Calculate the sizes of each set
only_wtvSEN = len(wtvSEN_ids - SENvDBL_ids - wtvDS_ids)
only_SENvDBL = len(SENvDBL_ids - wtvSEN_ids - wtvDS_ids)
only_wtvDS = len(wtvDS_ids - wtvSEN_ids - SENvDBL_ids)
intersection_wtvSEN_SENvDBL = len(wtvSEN_ids & SENvDBL_ids - wtvDS_ids)
intersection_wtvSEN_wtvDS = len(wtvSEN_ids & wtvDS_ids - SENvDBL_ids)
intersection_SENvDBL_wtvDS = len(SENvDBL_ids & wtvDS_ids - wtvSEN_ids)
intersection_all = len(wtvSEN_ids & SENvDBL_ids & wtvDS_ids)

# Custom colors for each set
set_colors = ('blue', 'green', 'red')

# Create a Venn diagram with custom colors and labels
venn3(subsets=(only_wtvSEN, only_SENvDBL, intersection_wtvSEN_SENvDBL,
               only_wtvDS, intersection_wtvSEN_wtvDS, intersection_SENvDBL_wtvDS,
               intersection_all),
      set_labels=('WT vs. Sen1D\n' + str(len(wtvSEN_ids)),
                  'WT vs. Dbl8D\n' + str(len(SENvDBL_ids)),
                  'WT vs. Sen1DDbl8D\n' + str(len(wtvDS_ids))),
      set_colors=set_colors)

plt.title("Customized Venn Diagram of IDs")
plt.show()







wtvSEN1 = wtvSEN[(wtvSEN['chro'] == 'I')]
wtvSEN2 = wtvSEN[(wtvSEN['chro'] == 'II')]
wtvSEN3 = wtvSEN[(wtvSEN['chro'] == 'III')]

wtvDS1 = wtvDS[(wtvDS['chro'] == 'I')]
wtvDS2 = wtvDS[(wtvDS['chro'] == 'II')]
wtvDS3 = wtvDS[(wtvDS['chro'] == 'III')]

SENvDBL1 = SENvDBL[(SENvDBL['chro'] == 'I')]
SENvDBL2 = SENvDBL[(SENvDBL['chro'] == 'II')]
SENvDBL3 = SENvDBL[(SENvDBL['chro'] == 'III')]

DBLvDS1 = DBLvDS[(DBLvDS['chro'] == 'I')]
DBLvDS2 = DBLvDS[(DBLvDS['chro'] == 'II')]
DBLvDS3 = DBLvDS[(DBLvDS['chro'] == 'III')]




dtcp1yay_aa_list = wtvSEN1['ID'].to_list()
dtcp2yay_aa_list = wtvSEN2['ID'].to_list()
dtcp3yay_aa_list = wtvSEN3['ID'].to_list()

gene1_list = wtvDS1['ID'].to_list()
gene2_list = wtvDS2['ID'].to_list()
gene3_list = wtvDS3['ID'].to_list()



def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

d1vg1 = intersection (dtcp1yay_aa_list, gene1_list)
d2vg2 = intersection (dtcp2yay_aa_list, gene2_list)
d3vg3 = intersection (dtcp3yay_aa_list, gene3_list)


len(gene1_list) + len(gene2_list) +len(gene3_list)
len(dtcp1yay_aa_list) + len(dtcp2yay_aa_list) +len(dtcp3yay_aa_list)
len(d1vg1) + len(d2vg2) +len(d3vg3)


wtvDS1aa = wtvDS1[~wtvDS1['ID'].isin(d1vg1)]
wtvSEN1 = wtvSEN1[~wtvSEN1['ID'].isin(d1vg1)]

def fucker(lst1, lst2):
    lst3 = [value for value in lst1 if not value in lst2]
    return lst3

d1vg1_NOT = intersection (dtcp1yay_aa_list, gene1_list)
d2vg2_NOT = intersection (dtcp2yay_aa_list, gene2_list)
d3vg3_NOT = intersection (dtcp3yay_aa_list, gene3_list)


len(d1vg1_NOT) + len(d2vg2_NOT) +len(d3vg3_NOT)


#%%
    
#let's keep the genes that are in common between the two 
wtvDS1 = wtvDS1[~wtvDS1['ID'].isin(d1vg1)]
wtvSEN1 = wtvSEN1[~wtvSEN1['ID'].isin(d1vg1)]


#what i would like to do now is to check rpelicon of the DE genes
def Create_df(df,dr,ef,er,w):
    df = pd.read_csv(df) # Reads a csv file

    df['chro'] = df['chro'].replace({"I": "chr1", "II": "chr2", "III": "chr3", 'mitochondrial': 'mito'})
    df['count'].replace(0,1, inplace = True) # This replaces zeros with ones. inplace = true saves the value to the dataframe (not useing a copy)
    df.rename(columns = {"count" : "df_count"}, inplace = True) # renames the colunm to appropriate counts
    
    dr = pd.read_csv(dr, usecols=[2]) # Read only the counts column from the next file

    dr['count'].replace(0,1, inplace = True)
    dr.rename(columns = {'count' : 'dr_count'}, inplace = True) # renames the colunm to appropriate counts 
    
    ef = pd.read_csv(ef, usecols=[2])
    ef['count'].replace(0,1, inplace = True)   
    ef.rename(columns = {'count' : 'ef_count'}, inplace = True)
    
    er = pd.read_csv(er, usecols=[2])
    er['count'].replace(0,1, inplace = True)
    er.rename(columns = {'count' : 'er_count'}, inplace = True)
    
    all_data = pd.concat([df, dr, ef, er], axis=1, join='outer') # Create a single dataframe by merging the 4 dataframes
    
    

    

  
# Here we add new colums to the dataframe
# First we normalise the data. ie each line in the column is simply the value divided by the sum of the column
# Note: pandas automatically itterate through the rows.    
    all_data['norm_df'] = all_data['df_count']/all_data['df_count'].sum()
    all_data['norm_dr'] = all_data['dr_count']/all_data['dr_count'].sum()
    all_data['norm_ef'] = all_data['ef_count']/all_data['ef_count'].sum()
    all_data['norm_er'] = all_data['er_count']/all_data['er_count'].sum()

# Next we calculate the ratios for each strand and assign to a new colunm
    all_data['ratio_delta_f'] = all_data['norm_df']/(all_data['norm_df'] + all_data['norm_ef'])
    all_data['ratio_delta_r'] = all_data['norm_dr']/(all_data['norm_dr'] + all_data['norm_er'])
    all_data['ratio_epsilon_f'] = all_data['norm_ef']/(all_data['norm_ef'] + all_data['norm_df'])
    all_data['ratio_epsilon_r'] = all_data['norm_er']/(all_data['norm_er'] + all_data['norm_dr'])


# Now we have column for pol delta useage for the duplex
    all_data['d_usage'] = (all_data['ratio_delta_f'] + all_data['ratio_delta_r']) / 2

# now we a column for the percentage of right-moving forks
   # all_data['right_forks']  = all_data['ratio_epsilon_f']*2 - 1
    all_data['right_forks'] = (all_data['norm_ef'] - all_data['norm_er'])/ (all_data['norm_ef'] + all_data['norm_er'])
    
    
# now we will produce a new colum for each sliding window average for each of the calculated columns
# Note: centre = True, means that the data is summed from both left and right. False means its the sum of the last of the number of values.
    
    all_data['smoo_ratio_d_f'] = all_data['ratio_delta_f'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_d_r'] = all_data['ratio_delta_r'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_e_f'] = all_data['ratio_epsilon_f'].rolling(window = w, center=True).mean()
    all_data['smoo_ratio_e_r'] = all_data['ratio_epsilon_r'].rolling(window = w, center=True).mean()
    all_data['smoo_d_usage'] = all_data['d_usage'].rolling(window = w, center=True).mean()
    all_data['smoo_right_forks'] = all_data['right_forks'].rolling(window = w, center=True).mean()
    

# now create a differential collunm and a discreate orgin (centre = 1 bin) column.
# Note the script removes differentials not present in both strands and which are singular 
    
    # take unsmoothed pol usage data into two pandas arrays
    ser_etop = all_data['ratio_epsilon_f']
    ser_ebot = all_data['ratio_epsilon_r']

    # Calculate the differentials 
    ser_topd = ser_etop.diff()
    ser_botd = ser_ebot.diff()
    # Reverse the sign on the bottom strand differentials to match top strand
    ser_botd = ser_botd*-1

    # curently dont use a roling average, but here they are if needed
    #ser_topd = ser_topd.rolling(window = 3, center=True).mean()
    #ser_botd = ser_botd.rolling(window = 3, center=True).mean()

    # Removing all the negative valuse to zero
    ser_topd.clip(0, 1, inplace=True)
    ser_botd.clip(0, 1, inplace=True)
    
    # If the value in one is zero, then seting it in zero in both datasets (this removes a lot of noise)
    for i in range(len(ser_topd)):
        if ser_topd.iloc[i] == 0:
            ser_botd.iloc[i] = 0
    for i in range(len(ser_botd)):
        if ser_botd.iloc[i] == 0:
            ser_topd.iloc[i] = 0
    
    # Now we add the two things together and divide by two - i.e we have an average origin activity
    ser_cumd = (ser_topd + ser_botd)/2     

    # Now we want to calculate the quantile of all the positive data.
    templist = np.array([])
    for i in range(1,len(ser_cumd)):
        if ser_cumd.iloc[i] != 0: 
            templist = np.append(templist, ser_cumd.iloc[i])
    # set a cutoff threshold at of the top 10% of all positive values
    cutoff = np.quantile(templist, 0.9)

    # now if a single +ve value (i.e at least one zero either side) is below this cutoff, then set to zero. This again removes noise.
    for i in range(1, len(ser_cumd)-1):
        if ser_cumd.iloc[i] != 0:
            if ser_cumd.iloc[i-1] == 0 and ser_cumd.iloc[i+1] == 0 and ser_cumd.iloc[i]< cutoff:
                ser_cumd.iloc[i] = 0

    # Now we save the data back to the daframe. This is a list of differentials (origin activity) per bin.
    all_data['differentials'] = ser_cumd

    # Next we want to position "zonal" origins to a single point (bin) and give them a cumulative value (efficiency of the zone)
    start=0
    stop=0
    origins = ser_cumd.clip(0, 0) # creates a new pd.series of zero values same length of ser_cumd
    for i in range(len(ser_cumd)):
        if i <= stop: # simply prevents itterating over the non-zero block
            continue # continue goes back to the for loop
        else:
            if ser_cumd.iloc[i] != 0:
                start = i        
                for z in range(start+1, len(ser_cumd)): # find the ned of the non-zero block
                    if ser_cumd.iloc[z] == 0:
                        stop = z
                        tot = 0
                        tem = 0
                        tem_loc = 0
                        for v in range(i,z): # adds the non-zero block numbers together and identifies the location of the highest value
                            tot = tot + ser_cumd.iloc[v]
                            if ser_cumd.iloc[v] > tem:
                                tem = ser_cumd.iloc[v]
                                tem_loc = v
                        origins.iloc[tem_loc] = tot # adds the total of the non-zero block to the bin with the highest individula differential
                        break

    all_data['origin'] = origins

# create three separate dataframes for the three chromosomes and return these too
 
    chrI = all_data.loc[all_data['chro']=='chr1']
    chrII = all_data.loc[all_data['chro']=='chr2']
    chrIII = all_data.loc[all_data['chro']=='chr3']
    return chrI, chrII, chrIII, all_data


owtchrI, owtchrII, owtchrIII, owt = Create_df('RZ259-RhArepr-d.e1.f-w300.count.csv', 'RZ259-RhArepr-d.e1.r-w300.count.csv', 'RZ267-RhArepr-e.e1.f-w300.count.csv', 'RZ267-RhArepr-e.e1.r-w300.count.csv', 5)
#sen1chrI, sen1chrII, sen1chrIII 


owt1 = owt[(owt['chro'] == 'chr1')]
owt2 = owt[(owt['chro'] == 'chr2')]
owt3 = owt[(owt['chro'] == 'chr3')]

owt1n = owt1[owt1['origin'] > 0.229]
owt2n = owt2[owt2['origin'] > 0.229]
owt3n = owt3[owt3['origin'] > 0.229]



def create_replicon_dataframe(df):
    # Sort the DataFrame by the 'pos' column
    df = df.sort_values(by='pos')

    # Initialize lists to store replicon information
    start_positions = []
    stop_positions = []
    replicon_ids = []

    # Iterate through the sorted DataFrame to find replicons
    for i in range(len(df) - 1):
        start_positions.append(df.iloc[i]['pos'])
        stop_positions.append(df.iloc[i + 1]['pos'])
        replicon_ids.append(i + 1)  # ID starts from 1

    # Create a new DataFrame from the replicon information
    replicon_df = pd.DataFrame({
        'start': start_positions,
        'stop': stop_positions,
        'id': replicon_ids
    })

    return replicon_df

repliconchr1 = create_replicon_dataframe(owt1n)
repliconchr2 = create_replicon_dataframe(owt2n)
repliconchr3 = create_replicon_dataframe(owt3n)



def assign_id_to_stalls(replicon_df, gene_df, b):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if (gene_row['start'] >= replicon_row['start'] and
                gene_row[b] <= replicon_row['stop']):
                matching_replicon_ids.append(replicon_row['id'])
                gene_list.append(gene_row['ID'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'id': matching_replicon_ids, 'ID': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='ID', how='left')
    print(id_and_ID)
    
    

   # gene_df = gene_df.append({'id':matching_replicon_ids}, ignore_index=True)
  #  matching_replicon_ids = pd.DataFrame(matching_replicon_ids)
   # print(matching_replicon_ids)
    #gene_df = gene_df.concat({'id':matching_replicon_ids}, ignore_index=True)
    #filtered_replicon_df = replicon_df[replicon_df['id'].isin(matching_replicon_ids)]

    return gene_df

gene1 = assign_id_to_stalls(repliconchr1, gene1, 'end')
gene2 = assign_id_to_stalls(repliconchr2, gene2, 'end')
gene3 = assign_id_to_stalls(repliconchr3, gene3, 'end')


con1 = assign_id_to_stalls(repliconchr1, con1, 'end')
con2 = assign_id_to_stalls(repliconchr2, con2, 'end')
con3 = assign_id_to_stalls(repliconchr3, con3, 'end')


wtvDS1 = assign_id_to_stalls(repliconchr1, wtvDS1, 'stop')
wtvDS2 = assign_id_to_stalls(repliconchr2, wtvDS2, 'stop')
wtvDS3 = assign_id_to_stalls(repliconchr3, wtvDS3, 'stop')

wtvSEN1 = assign_id_to_stalls(repliconchr1, wtvSEN1, 'stop')
wtvSEN2 = assign_id_to_stalls(repliconchr2, wtvSEN2, 'stop')
wtvSEN3 = assign_id_to_stalls(repliconchr3, wtvSEN3, 'stop')

SENvDBL1
SENvDBL1 = assign_id_to_stalls(repliconchr1, SENvDBL1, 'stop')
SENvDBL2 = assign_id_to_stalls(repliconchr2, SENvDBL2, 'stop')
SENvDBL3 = assign_id_to_stalls(repliconchr3, SENvDBL3, 'stop')


DBLvDS1
DBLvDS1 = assign_id_to_stalls(repliconchr1, DBLvDS1, 'stop')
DBLvDS2 = assign_id_to_stalls(repliconchr2, DBLvDS2, 'stop')
DBLvDS3 = assign_id_to_stalls(repliconchr3, DBLvDS3, 'stop')



wtvDS1 = wtvDS1.rename(columns={'stop': 'end'})
wtvDS2 = wtvDS2.rename(columns={'stop': 'end'})
wtvDS3 = wtvDS3.rename(columns={'stop': 'end'})


wtvSEN1 = wtvSEN1.rename(columns={'stop': 'end'})
wtvSEN2 = wtvSEN2.rename(columns={'stop': 'end'})
wtvSEN3 = wtvSEN3.rename(columns={'stop': 'end'})

SENvDBL1 = SENvDBL1.rename(columns={'stop': 'end'})
SENvDBL2 = SENvDBL2.rename(columns={'stop': 'end'})
SENvDBL3 = SENvDBL3.rename(columns={'stop': 'end'})

DBLvDS1 = DBLvDS1.rename(columns={'stop': 'end'})
DBLvDS2 = DBLvDS2.rename(columns={'stop': 'end'})
DBLvDS3 = DBLvDS3.rename(columns={'stop': 'end'})




wtvSEN1 = wtvSEN[(wtvSEN['chro'] == 'I')]
wtvSEN2 = wtvSEN[(wtvSEN['chro'] == 'II')]
wtvSEN3 = wtvSEN[(wtvSEN['chro'] == 'III')]

wtvDS1 = wtvDS[(wtvDS['chro'] == 'I')]
wtvDS2 = wtvDS[(wtvDS['chro'] == 'II')]
wtvDS3 = wtvDS[(wtvDS['chro'] == 'III')]

def Chromosome_plot (data1, data2, data3, data4, featurex, genee, thing, c, con, odata1):
    ff, (ax1, ax2, ax3, ax5) = plt.subplots(4,1, sharex=True, sharey=True)
    #ax1.set_title('WT')
    ax1.plot(data1['pos'], data1['norm'], color ='black', alpha=0.8)
    ax1.plot(data1['pos'], data1['norm2'], color ='dimgrey', alpha=0.8)
    ax2.plot(data2['pos'], data2['norm'], color ='steelblue', alpha=0.8)
    ax2.plot(data2['pos'], data2['norm2'], color ='deepskyblue', alpha=0.6)
  #  axb.plot(data3['pos'], data3['norm'], color ='orange', alpha=0.8)
  #  axb.plot(data3['pos'], data3['norm2'], color ='gold', alpha=0.8)
    ax3.plot(data4['pos'], data4['norm'], color ='tomato', alpha=0.8)
    ax3.plot(data4['pos'], data4['norm2'], color ='coral', alpha=0.8)
   # ax1.set_ylim(0,200)
    ax1.set_ylabel('NRD')
    
    for be in odata1.itertuples(index=False, name = None):
        if be[23] > 0.228:
            ax1.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
            ax2.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
            ax3.axvspan(be[1], be[1], 0.0, 1, color = 'black', alpha =0.7)
            

    ax5.set_xlabel('Chromosome position')

                  
    for fe in featurex.itertuples(index=False, name=None):
      #  ax5.annotate(fe[0], xy = [fe[2],0.45])  
        if fe[5] == 'reverse':
            if fe[7] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                if fe[6] == 'leftward':
                    ax5.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                
            #    ax5.annotate(fe[0], xy = [fe[2],-0.15])    
            
        elif fe[5] == 'forward':
            if fe[7] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                if fe[6] == 'rightward':
                
                    ax5.axvspan(fe[2],fe[3],0.5,0.8,color="tomato",alpha=0.5)
                    ax5.set_ylabel('Gene annotations')


    for gee in thing.itertuples(index=False, name=None):
        if gee[7] == 'gene':
            if gee[10] == '-':
                if gee[2] > 0:
                    ax5.axvspan(gee[8],gee[9],0.3,0.5,color="black",alpha=0.99)
                if gee[2] < 0:
                        ax5.axvspan(gee[8],gee[9],0.3,0.5,color="black",alpha=0.3)
            elif gee[10] == '+':
                if gee[2] > 0:
                    ax5.axvspan(gee[8],gee[9],0.5,0.8,color="black",alpha=0.99)
                if gee[2] < 0:
                    ax5.axvspan(gee[8],gee[9],0.5,0.8,color="black",alpha=0.3)
                
            
    for ge in genee.itertuples(index=False, name=None):
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax5.axvspan(ge[4],ge[5],0.2,0.3,color="rosybrown",alpha=0.3)
            elif ge[7] == 'forward':
                ax5.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)

    for co in con.itertuples(index=False, name=None):
        
        if co[4] == 'forward':
            ax5.axvspan(co[2],co[3],0.5,0.8,color="seagreen",alpha=0.3)
        elif co[4] == 'reverse':
            ax5.axvspan(co[2],co[3],0.3,0.5,color="seagreen",alpha=0.3)
    return ff

chromosome1 = Chromosome_plot(wtchr1, sen1chr1, dbl8chr1, dschr1, gene1, feat1, wtvSEN1, 'I', con1, owt1n)



def building_layers(replicon_df, gene_df, delta):
    
    matching_replicon_ids = []
    mathcing_rep_and_d_id = []
    topdeltapeaks = []
  #  matching_gene_ids = []  # New list to store gene IDs
    
    # Iterate over rows in gene_df
    #this line of code asks if a stall gene occurs within a given replicon and iterates, 
    #creates boolean mask
    for _, gene_row in gene_df.iterrows():
        # Check the condition using vectorized operations
        mask = (replicon_df['start'] <= gene_row['start']) & (gene_row['end'] <= replicon_df['stop'])
   #     print(mask)
        
        # Extract matching replicon IDs and add to the list
        #this takes which replicon_id the stall occurs in 
        
        matching_replicon_ids.extend(replicon_df.loc[mask, 'id'].tolist())
        

    # Create a DataFrame with matching IDs and 'YES' for 'stall_containing'
    temp = pd.DataFrame({'id': matching_replicon_ids, 'stall_containing': 'YES'})
    
    
    # Merge the temporary DataFrame with replicon_df
    replicon_df = replicon_df.merge(temp, on='id', how='left')

    # Fill NaN values in 'stall_containing' with 'NO'
    replicon_df['stall_containing'] = replicon_df['stall_containing'].fillna('NO')
    
       

    for _, delta_row in delta.iterrows():
        

        mask2 = (replicon_df['start'] <= delta_row['start']) & (delta_row['stop'] <= replicon_df['stop'])

        mathcing_rep_and_d_id.extend(replicon_df.loc[mask2, 'id'].tolist())

    temp2 = pd.DataFrame({'id': mathcing_rep_and_d_id, 'DE_containing': 'YES'})
    replicon_df = replicon_df.merge(temp2, on='id', how='left')
    replicon_df['DE_containing'] = replicon_df['DE_containing'].fillna('NO')
    
    
   
    #what i need to do now is add another later, whereby if mask = true, i also add peak_no to list
    
    df_no_duplicates = replicon_df.drop_duplicates(subset='id', keep='first')

        


    return df_no_duplicates

# Example usage
repliconchr1 = building_layers(repliconchr1, gene1, wtvDS1)
repliconchr2 = building_layers(repliconchr2, gene2, wtvDS2)
repliconchr3 = building_layers(repliconchr3, gene3, wtvDS3)

repliconchr1 = building_layers(repliconchr1, gene1, wtvSEN1)
repliconchr2 = building_layers(repliconchr2, gene2, wtvSEN2)
repliconchr3 = building_layers(repliconchr3, gene3, wtvSEN3)

repliconchr1 = building_layers(repliconchr1, gene1, SENvDBL1)
repliconchr2 = building_layers(repliconchr2, gene2, SENvDBL2)
repliconchr3 = building_layers(repliconchr3, gene3, SENvDBL3)

repliconchr1 = building_layers(repliconchr1, gene1, DBLvDS1)
repliconchr2 = building_layers(repliconchr2, gene2, DBLvDS2)
repliconchr3 = building_layers(repliconchr3, gene3, DBLvDS3)


repliconchr1 = create_replicon_dataframe(owt1n)
repliconchr2 = create_replicon_dataframe(owt2n)
repliconchr3 = create_replicon_dataframe(owt3n)



repliconchr1 = building_layers(repliconchr1, con1, wtvDS1)
repliconchr2 = building_layers(repliconchr2, con2, wtvDS2)
repliconchr3 = building_layers(repliconchr3, con3, wtvDS3)


repliconchr1 = building_layers(repliconchr1, con1, wtvSEN1)
repliconchr2 = building_layers(repliconchr2, con2, wtvSEN2)
repliconchr3 = building_layers(repliconchr3, con3, wtvSEN3)


repliconchr1 = building_layers(repliconchr1, con1, SENvDBL1)
repliconchr2 = building_layers(repliconchr2, con2, SENvDBL2)
repliconchr3 = building_layers(repliconchr3, con3, SENvDBL3)

#here
repliconchr1 = building_layers(repliconchr1, con1, DBLvDS1)
repliconchr2 = building_layers(repliconchr2, con2, DBLvDS2)
repliconchr3 = building_layers(repliconchr3, con3, DBLvDS3)





repliconchr1 = create_replicon_dataframe(owt1n)
repliconchr2 = create_replicon_dataframe(owt2n)
repliconchr3 = create_replicon_dataframe(owt3n)




count = 0

for x in repliconchr3.itertuples():
        if (x[4] == 'YES'):
            count += 1
            print(count)



count = 0

for x in repliconchr3.itertuples():
    if (x[5] == 'NO'):
        if (x[4] == 'YES'):
            count += 1
            print(count)

        
#####for stalls 
#always yes for stall
  #wtvds      
#chr1 = 33, no = 19
#chr2= 35, n0 = 11
#chr3 = 17, no =9

#wtvsen1
#chr1 = 24, no = 28
#chr2= 18, n0 = 28
#chr3 = 11, no =15

#sen vs dbl8
#chr1 = 18, no = 34
#chr2= 11, n0 = 35
#chr3 = 9, no =17

#dbl8 vs DS
#chr1 = 21, no = 31
#chr2= 19, n0 = 27
#chr3 = 13, no =13

###for controls
  #wtvds      
#chr1 = 46, no = 15
#chr2= 34, no =9
#chr3 = 16, no= 5

  #wtvsen      
#chr1 = 24, no = 37
#chr2= 26, no = 17
#chr3 = 9, no = 12

  #sen vs dbl      
#chr1 = 17, no = 44
#chr2= 12, no =31
#chr3 = 6, no =15

  #dbl8 vs ds      
#chr1 = 28
#chr2= 25
#chr3 = 10




count = 0

for x in repliconchr1.itertuples():
    if (x[4] == 'YES'):
            count += 1
            print(count)

#genes in r1: 61
#genes in r2 : 43
#genes in r3 = 21
 
#controls       
groups = ['WT vs. DS', 'WT vs. SEN1', 'WT vs. DBL8']
values1 = [96/125*100, 59/125*100, 35/125*100]
values2 = [29/125*100, 66/125*100, 90/125*100]

#stalls
groups = ['WT vs. DS', 'WT vs. SEN1', 'SEN1 vs. DBL8']
values1 = [85/124*100, 53/124*100, 38/124*100]
values2 = [39/124*100, 71/124*100, 86/124*100]


fig, ax = plt.subplots()

# Stacked bar chart
ax.bar(groups, values1, color = "seagreen", alpha = 0.8, label = "ddifferentiall expressed")
ax.bar(groups, values2, bottom = values1, color = "black", alpha = 0.8,  label = "equally expressed")

for bar in ax.patches:
  ax.text(bar.get_x() + bar.get_width() / 2,
          bar.get_height() / 2 + bar.get_y(),
          round(bar.get_height()), ha = 'center',
          color = 'w', weight = 'bold', size = 10)

total_values = np.add(values1, values2)

# Total values labels
for i, total in enumerate(total_values):
  ax.text(i, total + 0.5, round(total),
          ha = 'center', weight = 'bold', color = 'black')


#ax.legend()
plt.show()




import numpy as np
from scipy.stats import chi2_contingency

# Observed frequencies
observed = np.array([[125 * 0.43, 124 * 0.47],
                     [125 * 0.31, 124 * 0.28],
                     [125 * 0.69, 124 * 0.77]])

# Perform chi-squared test
chi2, p, dof, expected = chi2_contingency(observed)

print("Chi-squared statistic:", chi2)
print("p-value:", p)
print("Degrees of freedom:", dof)
print("Expected frequencies:")
print(expected)

#%%

wtvSEN = wtvSEN.merge(Genes2, on='ID', how='left')
SENvDBL = SENvDBL.merge(Genes2, on='ID', how='left')
wtvDS = wtvDS.merge(Genes2, on='ID', how='left')
#DBLvDS = DBLvDS.merge(Genes2, on="ID", how='left')

Genes2 = pd.merge(Genes2, wtvSEN[['ID', 'log2FoldChange']], on='ID', how='left')
Genes2 = pd.merge(Genes2, SENvDBL[['ID', 'log2FoldChange']], on='ID', how='left')
Genes2 = pd.merge(Genes2, wtvDS[['ID', 'log2FoldChange']], on='ID', how='left')


Genes2.rename(columns={'log2FoldChange_x':'log2FoldChange_sen', 'log2FoldChange_y':'log2FoldChange_dbl','log2FoldChange':'log2FoldChange_ds'}, inplace = True)

Genes2['log2FoldChange_sen'].fillna(0, inplace=True)
Genes2['log2FoldChange_dbl'].fillna(0, inplace=True)
Genes2['log2FoldChange_ds'].fillna(0, inplace=True)

Genes2.loc[Genes2['chro'] == 'I', 'chro'] = 'chr1'
Genes2.loc[Genes2['chro'] == 'II', 'chro'] = 'chr2'
Genes2.loc[Genes2['chro'] == 'III', 'chro'] = 'chr3'


Genes2y = Genes2[Genes2['type'].isin(['gene'])]

firstPA['relativepos'] = (firstPA['pos'] - firstPA['start_y'])/(firstPA['stop'] - firstPA['start_y'])

def score(file,stalls):
    
    new_df = pd.DataFrame()
    for b in np.arange(0.0, 1, 0.1):
        new_df[f'{b:.1f}'] = b 

    for i in range(len(stalls)):
        tempstart = stalls.iloc[i]["start"]
        tempend = stalls.iloc[i]["end"]
        tempchro = stalls.iloc[i]["chro"]
        tempendplus = tempend+10000
        

        tempsubset = file.loc[(file['start'] >= tempend) & (file['stop'] <= tempendplus) & (file['chro'] == tempchro)]
       # print(tempsubset)
       # for j in range(len(tempsubset)):
        tempsubset['realtivepos'] = ((tempsubset['stop'] - tempend)/(tempendplus - tempend)).round(1)
       # print(tempsubset)
        #for _, row in tempsubset.iterrows():
         #       relative_pos = row['realtivepos']
          #      log2FoldChange_ds_value = row['log2FoldChange_ds']
        for c in range(len(tempsubset)):
            deadrelly = tempsubset.iloc[c]['realtivepos']
            log2FoldChange_ds_value = tempsubset.iloc[c]['log2FoldChange_ds']
            
            if deadrelly in new_df.columns:
                print('a')
                new_df.loc[c, deadrelly] = log2FoldChange_ds_value

                # Update new_df based on relative_pos
               # print(new_df.columns)
  #              if relative_pos in new_df.columns:
              #      print('yes')
                    

        #these next two lines are just here for me to check the output
         #   if i <= 10:
         #       print(tempsubset)
        
    return new_df

#trytrytryagain = score(Genes2y, ggenes)

newccontrol['genotype'] = 'control'

compiled = pd.concat([newccontrol,ggenes])


def score(file,stalls):
    
    new_df = pd.DataFrame()

    for i in range(len(stalls)):
        tempID = stalls.iloc[i]["ID"]
        tempgeno = stalls.iloc[i]["genotype"]
        tempstart = stalls.iloc[i]["start"]
        tempend = stalls.iloc[i]["end"]
        tempchro = stalls.iloc[i]["chro"]
        tempendplus = tempend+10000
        

        tempsubset = file.copy().loc[(file['start'] >= tempend) & (file['stop'] <= tempendplus) & (file['chro'] == tempchro)]
        tempsubset['realtivepos'] = ((tempsubset['stop'] - tempend)/(tempendplus - tempend)).round(1)
        tempsubset['locus'] = tempID
        tempsubset['genotype'] = tempgeno
        
      #  if i <= 10:
       #         print(tempsubset)
                
        new_df = pd.concat([new_df,tempsubset])

        
    return new_df

trytrytryagain = score(Genes2y, compiled)
trytrytryagain["genotype"] = pd.Categorical(trytrytryagain["genotype"],
                                                 ordered = True,
                                                 categories = ["sen1D", "dbl8D", "sen1dbl8DD_unique", "control"])


trytrytryagainagainds = trytrytryagain.groupby(['locus', 'realtivepos', 
                                                'genotype']).agg({'log2FoldChange_ds': 'mean'})
#result_meanup = trytrytryagainagain.groupby(['realtivepos', 'genotype'])['log2FoldChange_ds'].mean()
#result_mean_df2 = result_meanup.reset_index()




custom_palette = ["steelblue", "orange", 'tomato', 'seagreen']  

sns.set_palette(custom_palette)




sns.lineplot(data=trytrytryagainagainds, x="realtivepos", y="log2FoldChange_ds", hue = 'genotype', legend=False)
plt.ylim(-0.7,0.7)
plt.show()


def scoredownstream(file,stalls):
    
    new_df = pd.DataFrame()

    for i in range(len(stalls)):
        tempID = stalls.iloc[i]["ID"]
        tempgeno = stalls.iloc[i]["genotype"]
        tempstart = stalls.iloc[i]["start"]
        tempend = stalls.iloc[i]["end"]
        tempchro = stalls.iloc[i]["chro"]
        tempendplus = tempstart-10000
        

        tempsubset = file.copy().loc[(file['stop'] <= tempstart) & (tempendplus <= file['start']) & (file['chro'] == tempchro)]
        tempsubset['realtivepos'] = ((tempsubset['stop'] - tempendplus)/(tempstart - tempendplus)).round(1)
        tempsubset['locus'] = tempID
        tempsubset['genotype'] = tempgeno
        
      #  if i <= 10:
       #         print(tempsubset)
                
        new_df = pd.concat([new_df,tempsubset])

        
    return new_df

trytrytryagaindownstream = scoredownstream(Genes2y, compiled)


trytrytryagaindownstream["genotype"] = pd.Categorical(trytrytryagaindownstream["genotype"],
                                                 ordered = True,
                                                 categories = ["sen1D", "dbl8D", "sen1dbl8DD_unique", "control"])




trytrytryagaindownstreamds = trytrytryagaindownstream.groupby(['locus', 'realtivepos', 'genotype']).agg({'log2FoldChange_sen': 'mean'})
#result_mean = trytrytryagaindownstream.groupby(['realtivepos', 'genotype'])['log2FoldChange_ds'].mean()
#result_mean_df = result_mean.reset_index()




sns.lineplot(data=trytrytryagaindownstreamds, x="realtivepos", y="log2FoldChange_sen", hue = 'genotype', legend=False)
plt.ylim(-0.7,0.7)
plt.show()


#%%



def Chromosome_plot (featurex, centro, telo, genee, c, con, geneeb):
    ff, (ax5) = plt.subplots(1,1, sharex=True)

                  
    for fe in featurex.itertuples(index=False, name=None):
      #  ax5.annotate(fe[0], xy = [fe[2],0.45])  
        if fe[5] == 'reverse':
            if fe[7] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                ax5.axvspan(fe[2],fe[3],0.3,0.5,color="tomato",alpha=0.5)
                
            #    ax5.annotate(fe[0], xy = [fe[2],-0.15])    
            
        elif fe[5] == 'forward':
            if fe[7] == 'sen1D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="steelblue",alpha=0.5)
            if fe[7] == 'dbl8D':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="orange",alpha=0.5)
            if fe[7] == 'sen1dbl8DD_unique':
                ax5.axvspan(fe[2],fe[3],0.5,0.8,color="tomato",alpha=0.5)
                ax5.set_ylabel('Gene annotations')
                
    for c in centro.itertuples(index=False, name=None):
            ax5.axvspan(c[0],c[1],0.1,0.15,color="black",alpha=0.8)
            
    for t in telo.itertuples(index=False, name=None):
            ax5.axvspan(t[0],t[1],0.1,0.15,color="grey",alpha=0.8)
            
    for ge in genee.itertuples(index=False, name=None):
        if ge[3] == 'protein coding':
            if ge[7] == 'reverse':
                ax5.axvspan(ge[4],ge[5],0.2,0.3,color="rosybrown",alpha=0.3)
            elif ge[7] == 'forward':
                ax5.axvspan(ge[4],ge[5],0.8,0.9,color="rosybrown",alpha=0.3)

    for geb in geneeb.itertuples(index=False, name=None):
        if geb[2] < 0:
            if geb[10] == '-':
                ax5.axvspan(geb[8],geb[8],0.3,0.5,color="black",alpha=0.3)
            elif geb[10] == '+':
                ax5.axvspan(geb[8],geb[8],0.5,0.8,color="black",alpha=0.3)

    for geb in geneeb.itertuples(index=False, name=None):
        if geb[2] > 0:
            if geb[10] == '-':
                ax5.axvspan(geb[8],geb[8],0.3,0.5,color="green",alpha=0.3)
            elif geb[10] == '+':
                ax5.axvspan(geb[8],geb[8],0.5,0.8,color="green",alpha=0.3)
                                          



    for co in con.itertuples(index=False, name=None):
        
        if co[4] == 'forward':
            ax5.axvspan(co[2],co[3],0.5,0.8,color="blue",alpha=0.3)
        elif co[4] == 'reverse':
            ax5.axvspan(co[2],co[3],0.3,0.5,color="blue",alpha=0.3)
    return ff

chromosome1 = Chromosome_plot(gene1, centro1, telo1, feat1, 'I', con1, wtvDS1)
chromosome2 = Chromosome_plot(gene2, centro2, telo2, feat2, 'II', con2, wtvSEN2)
chromosome2 = Chromosome_plot(gene3, centro3, telo3, feat3, 'III', con3, wtvSEN3)


#%%



def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes = genes.rename(columns={'Unnamed: 0': 'ID'})
    genes['ID'] = genes['ID'].str.replace('a', '.')

    print(genes)
    

    return genes

wtwt = Find("wtwtcounts.csv")
sensen = Find("sensencounts.csv")
dbldbl = Find("dbldblcounts.csv")
dsds = Find("doubleboudlecounts.csv")

Genes2['length_kb'] = (Genes2['stop'] - Genes2['start'])/1000



wtwt = wtwt.merge(Genes2, on='ID', how='left')
sensen = sensen.merge(Genes2, on='ID', how='left')
dbldbl = dbldbl.merge(Genes2, on='ID', how='left')
dsds = dsds.merge(Genes2, on="ID", how='left')

wtwttrytrytry =ggenes.merge(wtwt, on='ID', how='left')
#so now i need to form a function to sum total counts from raw HT seq output
#normalise to millions 
#so take path as the all one 


def Find(file):
    count_df = pd.read_csv(file, delimiter="\t", header=None)
    count_df_mapped = count_df[~count_df[0].str.startswith('__')]
    total_mapped_reads_per_sample = count_df_mapped[1].sum()
    total_mapped_reads_in_millions = total_mapped_reads_per_sample / 1e6
    print(total_mapped_reads_per_sample)
    return total_mapped_reads_in_millions

WT_HT = Find("RZ316_w_thia_wt_treated.txt")
WT_HTa = Find("RZ316a_w_thia_wt_treated.txt")
sen_HT = Find("RZ317_w_thia_sen_treated.txt")
sen_HTa = Find("RZ317a_w_thia_sen_treated.txt")
dbl_HT = Find("RZ318_w_thia_dbl_treated.txt")
dbl_HTa = Find("RZ318a_w_thia_dbl_treated.txt")
ds_HT = Find("RZ319_w_thia_ds_treated.txt")
ds_HTa = Find("RZ319a_w_thia_ds_treated.txt")


WT_HT = Find("RZ316_wo_thia_wt_treated_modified.txt")
WT_HTa = Find("RZ316a_wo_thia_wt_treated_modified.txt")
sen_HT = Find("RZ317_wo_thia_sen_treated_modified.txt")
sen_HTa = Find("RZ317a_wo_thia_sen_treated_modified.txt")
dbl_HT = Find("RZ318_wo_thia_dbl_treated_modified.txt")
dbl_HTa = Find("RZ318a_wo_thia_dbl_treated_modified.txt")
ds_HT = Find("RZ319_wo_thia_ds_treated_modified.txt")
ds_HTa = Find("RZ319a_wo_thia_ds_treated_modified.txt")

def calculate_FPKM(df, total_mapped_reads, total_mapped_readsa, col1, col2):
    
    
    df['FPKM_1'] = (df[col1] * 1e9) / (df['length_kb'] * total_mapped_reads)
    df['FPKM_2'] = (df[col2] * 1e9) / (df['length_kb'] * total_mapped_readsa)
    df['FPKM_av'] = (df['FPKM_1'] + df['FPKM_2']) / 2
    df.loc[df['chro'] == "I", 'chro'] = 'chr1'
    df.loc[df['chro'] == "II", 'chro'] = 'chr2'
    df.loc[df['chro'] == "III", 'chro'] = 'chr3'
    df.loc[df['strand'] == "+", 'coding_strand'] = 'forward'
    df.loc[df['strand'] == "-", 'coding_strand'] = 'reverse'


    return df

wtwt= calculate_FPKM(wtwt,WT_HT, WT_HTa, "WT", "WT.")
sensen =calculate_FPKM(sensen, sen_HT, sen_HTa, "ds", "ds.")
dbldbl =calculate_FPKM(dbldbl, dbl_HT, dbl_HTa, "ds", "ds.")
dsds= calculate_FPKM(dsds, ds_HT, ds_HTa,"ds", "ds.")


allthegoodg = [gene1, gene2, gene3]
allthestalls = pd.concat(allthegoodg)

allthegoodgc = [con1, con2, con3]
allthecontrolss = pd.concat(allthegoodgc)

#merged_df = pd.merge(df1, df2[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_custom_suffix'))
allthestallswt = pd.merge(allthestalls, wtwt[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_wt'))
allthestallsen = pd.merge(allthestalls, sensen[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_sen'))
allthestallsdbl = pd.merge(allthestalls, dbldbl[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_dbl'))
allthestallsds = pd.merge(allthestalls, dsds[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_ds'))

newccontrol = pd.merge(newccontrol, wtwt[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_wt'))

allthestalls = pd.merge(allthestalls, wtwt[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_wt'))
allthestalls = pd.merge(allthestalls, sensen[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_sen'))
allthestalls = pd.merge(allthestalls, dbldbl[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_dbl'))
allthestalls = pd.merge(allthestalls, dsds[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_ds'))






melted_df = allthestalls.melt(id_vars=['ID', 'gene', 'start', 'end', 'chro', 'coding_strand', 'stalled_fork',
       'genotype', 'FPKM_av', 'FPKM_av_sen', 'FPKM_av_dbl', 'FPKM_av_ds'],
                               var_name='source',
                               value_vars= ['FPKM_av', 'FPKM_av_sen', 'FPKM_av_dbl', 'FPKM_av_ds'],
                               value_name='FPKM_melt')



melted_df = allthestalls.melt(id_vars=['ID', 'gene', 'start', 'end', 'chro', 'coding_strand', 'stalled_fork', 'genotype'],
                               var_name='source_e',
                               value_vars= ['FPKM_av', 'FPKM_av_sen', 'FPKM_av_dbl', 'FPKM_av_ds'],
                               value_name='log2_epb')


melted_df['source_e'] = melted_df['source_e'].str.split('_').str[-1]





custom_palette = ["gray","steelblue", "orange", 'tomato' ]  

sns.set_palette(custom_palette)

sns.boxplot(
    data=melted_df, x="source_e", y="log2_epb",
    hue="source_e",showfliers=False
)


sns.violinplot(
    data=melted_df, x="genotype", y="log2_epb",
    hue="source_e",fill=False, legend=False
)

allthecontrolss = pd.merge(allthecontrolss, wtwt[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_wt'))
allthecontrolss = pd.merge(allthecontrolss, sensen[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_sen'))
allthecontrolss = pd.merge(allthecontrolss, dbldbl[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_dbl'))
allthecontrolss = pd.merge(allthecontrolss, dsds[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_ds'))

allthecontrolss['genotype'] = 'control'



df_all = pd.concat([allthecontrolss, allthestalls], axis=0)


melted_df_controls = allthecontrolss.melt(id_vars=['ID', 'chro', 'start', 'end', 'coding_strand', 'cpc', 'type'],
                               var_name='source_e',
                               value_vars= ['FPKM_av', 'FPKM_av_sen', 'FPKM_av_dbl', 'FPKM_av_ds'],
                               value_name='log2_epb')


sns.violinplot(
    data=melted_df_controls, x="source_e", y="log2_epb",
    hue="source_e",fill=False, legend=False
)


sns.boxplot(
    data=melted_df_controls, x="source_e", y="log2_epb",
    hue="source_e",showfliers=False
)



Genes2

#fuckfuck = pd.merge(allthestalls, Genes2[['ID', 'length_kb']], on='ID', how='inner', suffixes=('', '_ds'))


ffeaty.rename(columns={'Systematic ID':'ID'}, inplace = True)

ffeaty = pd.merge(ffeaty, wtwt[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_wt'))
ffeaty = pd.merge(ffeaty, sensen[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_sen'))
ffeaty = pd.merge(ffeaty, dbldbl[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_dbl'))
ffeaty = pd.merge(ffeaty, dsds[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_ds'))

melted_df_feat = ffeaty.melt(id_vars=['ID', 'Gene name', 'Product description', 'Feature type',
       'Start position', 'End position', 'Chromosome', 'Strand', 'chro'],
                               var_name='source_e',
                               value_vars= ['FPKM_av', 'FPKM_av_sen', 'FPKM_av_dbl', 'FPKM_av_ds'],
                               value_name='log2_epb')


melted_df_feat['log2_FPKM'] = np.log2(melted_df_feat['log2_epb'])

sns.swarmplot(
    data=melted_df_feat, x="source_e", y="log2_epb",
    hue="source_e",size=0.5
)




df_all = pd.concat([allthecontrolss, allthestalls], axis=0)

melted_df_feat = df_all.melt(id_vars=['ID', 'chro', 'start', 'end', 'coding_strand', 'cpc', 'type',
       'genotype',
       'gene', 'stalled_fork'],
                               var_name='source_e',
                               value_vars= ['FPKM_av', 'FPKM_av_sen', 'FPKM_av_dbl', 'FPKM_av_ds'],
                               value_name='FPKM')

dfbigbigbigbgi = melted_df_feat[~melted_df_feat.isin([-np.inf]).any(axis=1)]


sns.boxplot(
    data=melted_df_feat, x="genotype", y="FPKM",
    hue="source_e",showfliers=False
)

melted_df_feat['FPKM_log2'] = np.log2(melted_df_feat['FPKM_log2'])
                                      
sns.boxplot(
    data=melted_df_feat, x="genotype", y="FPKM",
    hue="source_e", showfliers=False, fill=False, gap=.3, legend=False)

# Overlay individual data points
sns.stripplot(
    data=melted_df_feat, x="genotype", y="FPKM",
    hue="source_e", dodge=True, marker="o", alpha=0.3, legend=False)

# Show the plot

plt.title('Box Plot with Individual Data Points')
plt.show()


#%%

#now now sweet chicken puff, you should check what the expression level of the differentiall
#espressed genes are. 

Genes2y = pd.merge(Genes2y, wtwt[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_wt'))
Genes2y = pd.merge(Genes2y, sensen[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_sen'))
Genes2y = pd.merge(Genes2y, dbldbl[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_dbl'))
Genes2y = pd.merge(Genes2y, dsds[['ID', 'FPKM_av']], on='ID', how='inner', suffixes=('', '_ds'))




melted_dfyy = Genes2y.melt(id_vars=['chro', 'type', 'start', 'stop', 'strand', 'ID', 'log2FoldChange_sen',
       'log2FoldChange_dbl', 'log2FoldChange_ds',],
                               var_name='source_e',
                               value_vars= ['FPKM_av', 'FPKM_av_sen', 'FPKM_av_dbl', 'FPKM_av_ds'],
                               value_name='FPKMs')


df_ds = melted_dfyy[melted_dfyy['log2FoldChange_ds'] != 0]
df_ds['FPKM_log2'] = np.log2(df_ds['FPKMs'])
                                      
sns.boxplot(
    data=df_ds, x="source_e", y="FPKM_log2",
    hue="source_e", showfliers=False, fill=False, legend=False)

# Overlay individual data points
sns.stripplot(
    data=df_ds, x="source_e", y="FPKM_log2",
    hue="source_e", dodge=True, marker="o", alpha=0.3, legend=False)

# Show the plot

plt.title('Box Plot with Individual Data Points')
plt.ylim(25,45)
plt.show()


df_dbl = melted_dfyy[melted_dfyy['log2FoldChange_dbl'] != 0]
df_dbl['FPKM_log2'] = np.log2(df_dbl['FPKMs'])
sns.boxplot(
    data=df_dbl, x="source_e", y="FPKM_log2",
    hue="source_e", showfliers=False, fill=False, legend=False)

# Overlay individual data points
sns.stripplot(
    data=df_dbl, x="source_e", y="FPKM_log2",
    hue="source_e", dodge=True, marker="o", alpha=0.3, legend=False)

# Show the plot

plt.title('Box Plot with Individual Data Points')
plt.show()


df_sen = melted_dfyy[melted_dfyy['log2FoldChange_sen'] != 0]
df_sen['FPKM_log2'] = np.log2(df_sen['FPKMs'])
sns.boxplot(
    data=df_sen, x="source_e", y="FPKM_log2",
    hue="source_e", showfliers=False, fill=False, legend=False)

# Overlay individual data points
sns.stripplot(
    data=df_sen, x="source_e", y="FPKM_log2",
    hue="source_e", dodge=True, marker="o", alpha=0.3, legend=False)

# Show the plot

plt.title('Box Plot with Individual Data Points')
plt.show()


#quickly try this
df_ds['spression'] = 'ds'
df_dbl['spression'] = 'dbl'
df_sen['spression'] ='sen'

whoopsie = pd.concat([df_ds, df_dbl, df_sen], axis=0)

sns.boxplot(
    data=whoopsie, x="spression", y="FPKM_log2",
    hue="source_e", showfliers=False, fill=False, legend=False)

# Overlay individual data points
sns.stripplot(
    data=whoopsie, x="spression", y="FPKM_log2",
    hue="source_e", dodge=True, marker="o", alpha=0.3, legend=False)

# Show the plot

plt.title('Box Plot with Individual Data Points')
plt.show()





dfbigbigbigbgi['fuck'] = 'stall&control'



whoopsie['genotype'] = 'dif_exp'





bigbird = pd.concat([whoopsie, dfbigbigbigbgi], axis=0)

bigbird.loc[bigbird['genotype'] == "sen1D", 'genotype'] = 'stall'
bigbird.loc[bigbird['genotype'] == "dbl8D", 'genotype'] = 'stall'
bigbird.loc[bigbird['genotype'] == "sen1dbl8DD_unique", 'genotype'] = 'stall'

sns.violinplot(data=bigbird, y="genotype", x="FPKM_log2", hue="genotype")

sns.boxplot(
    data=bigbird, x="genotype", y="FPKM_log2",
    hue="source_e", showfliers=False, fill=False, legend=False)

# Overlay individual data points
sns.stripplot(
    data=bigbird, x="genotype", y="FPKM_log2",
    hue="source_e", dodge=True, marker="o", alpha=0.3, legend=False)

# Show the plot

plt.title('Box Plot with Individual Data Points')
plt.show()


custom_palette = ["grey", "seagreen", 'tomato']  

sns.set_palette(custom_palette)
sns.boxplot(
    data=bigbird, y="genotype", x="FPKM_log2",
    hue="genotype", showfliers=False, fill=False, legend=False)
sns.stripplot(
    data=bigbird, y="genotype", x="FPKM_log2",
    hue="genotype", marker="o", alpha=0.15, legend=False)





#####
df_dbl['FPKM_log2'] = np.log2(df_dbl['FPKMs'])


prefixes = ["SPNCRNA", "SPSNORNA", "SPCC", "SPBP", "SPBC", "SPCP", "SPAC", "SPAP"]

# Initialize a dictionary to store the count for each prefix
prefix_count = {prefix: 0 for prefix in prefixes}

# Iterate over the 'ID' column and update the count for each prefix
for gene_name in df_ds['ID']:
    for prefix in prefixes:
        if gene_name.startswith(prefix):
            prefix_count[prefix] += 1
            break  # Stop checking prefixes once a match is found

# Print the count for each prefix
for prefix, count in prefix_count.items():
    print(f"Number of genes starting with {prefix}: {count}")
    

    
    
    
sample_data = {
    'df1': {'ncRAN': 355, 'protein coding': 257},
    'df2': {'ncRAN': 95, 'protein coding': 16},
    'df3': {'ncRAN': 176, 'protein coding': 113}
}

groups = ['WT vs. DS', 'WT vs. DBL', 'WT vs. SEN']
values1 = [355/612*100, 95/111*100, 176/289*100]
values2 = [257/612*100, 16/111*100, 113/289*100]



fig, ax = plt.subplots()

# Stacked bar chart
ax.bar(groups, values1, color = "#1D5B79", alpha = 0.8, label = "ddifferentiall expressed")
ax.bar(groups, values2, bottom = values1, color = "#A6CF98", alpha = 0.8,  label = "equally expressed")

for bar in ax.patches:
  ax.text(bar.get_x() + bar.get_width() / 2,
          bar.get_height() / 2 + bar.get_y(),
          round(bar.get_height()), ha = 'center',
          color = 'w', weight = 'bold', size = 10)

total_values = np.add(values1, values2)

# Total values labels
for i, total in enumerate(total_values):
  ax.text(i, total + 0.5, round(total),
          ha = 'center', weight = 'bold', color = 'black')


#ax.legend()
plt.show()


["#1D5B79", "#EF6262", "#A6CF98", "#557C55" ]
prefixes = ['ncRAN', 'snoRNA', 'protein coding']


#%%

#so we have a huge problem in that not all protein coding genes are aliigning oops FUCK 
#so let's write something whereby if "start" and "stop" of wtwt 
#overlaps with a gene body (where Genes2['type'] == 'gene')
#we tell python to asign those reads to that gene 
#i think we can do this simply in allthestalls


#fuckfuck = pd.merge(allthestalls, Genes2[['ID', 'length_kb']], on='ID', how='inner', suffixes=('', '_ds'))
#that verifies that the stall genes are in Genes2


def assign_id_to_stalls(replicon_df, gene_df):
    matching_replicon_ids = []
    gene_list = []

    for _, gene_row in gene_df.iterrows():
        for _, replicon_row in replicon_df.iterrows():
            if replicon_row['type'] == 'gene':
                if (gene_row['start'] >= replicon_row['start'] and
                    gene_row['stop'] <= replicon_row['stop']):
                    matching_replicon_ids.append(gene_row['FPKM_av'])
                    gene_list.append(replicon_df['ID'])
    print(len(matching_replicon_ids))
    print(len(gene_list))
    
    id_and_ID = pd.DataFrame({'fpkm': matching_replicon_ids, 'ID': gene_list})
    gene_df = gene_df.merge(id_and_ID, on='ID', how='left')
    print(id_and_ID)
    

    return gene_df

wtwtfuck = assign_id_to_stalls(Genes2, wtwt)
#later filter for nans in the column i'm going to make















