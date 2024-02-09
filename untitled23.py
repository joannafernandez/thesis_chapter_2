#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 21:46:09 2024

@author: patricfernandez
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

#%%
#Derine a function to make dataframes and normalise

def Collect_data(file1,file2,fact):
    forward = pd.read_csv(file1)
    reverse = pd.read_csv(file2)
    forward['score'] = forward['count'] + reverse['count']
 #   print(forward)
    print (forward['score'].sum())
    forward['norm']= (forward['score'] / forward['score'].sum())*fact
   # print(forward)

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



wtchr1, wtchr2, wtchr3, wt = Collect_data("1-RNApol2-ChIP-RZ316.cen.f-w50.count.csv",
                                      '1-RNApol2-ChIP-RZ316.cen.r-w50.count.csv', 1000000)

sen1chr1,sen1chr2,sen1chr3, sen1 = Collect_data("3-RNApol2-ChIP-RZ317.cen.f-w50.count.csv",
                                          '3-RNApol2-ChIP-RZ317.cen.r-w50.count.csv', 1000000)

dbl8chr1,dbl8chr2,dbl8chr3, dbl8 = Collect_data("5-RNApol2-ChIP-RZ318.cen.f-w50.count.csv",
                                          '5-RNApol2-ChIP-RZ318.cen.r-w50.count.csv', 1000000)

dschr1,dschr2,dschr3, ds = Collect_data("7-RNApol2-ChIP-RZ319.cen.f-w50.count.csv",
                                    '7-RNApol2-ChIP-RZ319.cen.r-w50.count.csv', 1000000)


#with rna seq data insteads
wtchr1, wtchr2, wtchr3, wt = Collect_data("10-RNAseq-RZ316-w-thia.cen.f-w50.count.csv",
                                      '10-RNAseq-RZ316-w-thia.cen.r-w50.count.csv', 1000000)

sen1chr1,sen1chr2,sen1chr3, sen1 = Collect_data("12-RNAseq-RZ317-w-thia.cen.f-w50.count.csv",
                                          '12-RNAseq-RZ317-w-thia.cen.r-w50.count.csv', 1000000)

dbl8chr1,dbl8chr2,dbl8chr3, dbl8 = Collect_data("14-RNAseq-RZ318-w-thia.cen.f-w50.count.csv",
                                          '14-RNAseq-RZ318-w-thia.cen.r-w50.count.csv', 1000000)

dschr1,dschr2,dschr3, ds = Collect_data("16-RNAseq-RZ319-w-thia.cen.f-w50.count.csv",
                                    '16-RNAseq-RZ319-w-thia.cen.r-w50.count.csv', 1000000)

###getting read infor from the repeats 

wtchr1, wtchr2, wtchr3, wt = Collect_data("316_RNAP2_chip_2.cen.f-w100.count.csv",
                                      '316_RNAP2_chip_2.cen.r-w100.count.csv', 1000000)

sen1chr1,sen1chr2,sen1chr3, sen1 = Collect_data("317_RNAP2_chip_2.cen.f-w100.count.csv",
                                          '317_RNAP2_chip_2.cen.r-w100.count.csv', 1000000)

dbl8chr1,dbl8chr2,dbl8chr3, dbl8 = Collect_data("318_RNAP2_chip_2.cen.f-w100.count.csv",
                                                "318_RNAP2_chip_2.cen.r-w100.count.csv", 1000000)

dschr1,dschr2,dschr3, ds = Collect_data("319_RNAP2_chip_2.cen.f-w100.count.csv",
                                    '319_RNAP2_chip_2.cen.r-w100.count.csv', 1000000)


#%%%

def Collect_data(file1,file2,fact, file1REP, file2REP):
    
    forward = pd.read_csv(file1)
    reverse = pd.read_csv(file2)
    forward['score'] = forward['count'] + reverse['count']
    forward['norm']= (forward['score'] / forward['score'].sum())*fact
    
    
    forwardREP = pd.read_csv(file1REP)
    reverseREP = pd.read_csv(file2REP)
    forward['score2'] = forwardREP['count'] + reverseREP['count']
    forward['norm2']= (forward['score2'] / forward['score2'].sum())*fact
    
    forward['normalised'] = (forward['norm'] + forward['norm2'])/2
    
    print(forward)   


    chr1= forward[(forward["chro"] == "chr1")]
    chr2= forward[(forward["chro"] == "chr2")]
    chr3= forward[(forward["chro"] == "chr3")]
    return chr1,chr2,chr3, forward


wtchr1, wtchr2, wtchr3, wt = Collect_data("10-RNAseq-RZ316-w-thia.cen.f-w50.count.csv",
                                      '10-RNAseq-RZ316-w-thia.cen.r-w50.count.csv', 1000000,
                                      "316T_RNAseq_bowtie.cen_replicte2.f-w50.count.csv", 
                                      '316T_RNAseq_bowtie.cen_replicte2.r-w50.count.csv')


sen1chr1,sen1chr2,sen1chr3, sen1 = Collect_data("12-RNAseq-RZ317-w-thia.cen.f-w50.count.csv",
                                          '12-RNAseq-RZ317-w-thia.cen.r-w50.count.csv', 1000000,
                                          "317T_RNAseq_bowtie.cen_replicte2.f-w50.count.csv",
                                          '317T_RNAseq_bowtie.cen_replicte2.r-w50.count.csv',)

dbl8chr1,dbl8chr2,dbl8chr3, dbl8 = Collect_data("14-RNAseq-RZ318-w-thia.cen.f-w50.count.csv",
                                          '14-RNAseq-RZ318-w-thia.cen.r-w50.count.csv', 1000000,
                                          "318T_RNAseq_bowtie.cen_replicte2.f-w50.count.csv",
                                          "318T_RNAseq_bowtie.cen_replicte2.r-w50.count.csv")

dschr1,dschr2,dschr3, ds = Collect_data("16-RNAseq-RZ319-w-thia.cen.f-w50.count.csv",
                                    '16-RNAseq-RZ319-w-thia.cen.r-w50.count.csv', 1000000,
                                    "319T_RNAseq_bowtie.cen_replicte2.f-w50.count.csv",
                                    '319T_RNAseq_bowtie.cen_replicte2.r-w50.count.csv')

df_all = pd.concat([df1_subset, df2_subset, df3_subset, df4_subset], axis=0)

# Calculate the correlation matrix
correlation_matrix = new_df.corr()

# Create a heatmap
s#ns.set(style="white")  # Set the style, optional
plt.figure(figsize=(10, 8))  # Set the figure size, optional
sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", fmt=".2f", linewidths=.5, vmin= 0.50, vmax = 1)

# Set the title
plt.title('Correlation Heatmap of Normalized Values')

#%%

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


Genes2_3 = Genes2[(Genes2['type'] == 'three_prime_UTR')]
Genes2_3['ID'] = Genes2_3['ID'].str.split(':').str[0]
Genes2_3['ID'] = Genes2_3['ID'].apply(lambda x: '.'.join(x.split('.')[:2]))

Genes2_3['length'] = Genes2_3['stop'] - Genes2_3['start']
mean_lenght_of_utr3 = Genes2_3['length'].mean()

sns.boxplot(data=Genes2_3, x='length')


def Find(file):
    genes = pd.read_csv(file, delimiter="\t")
    print(genes)

    return  genes

ggenes = Find("dbl8_stall_sites_direction.txt")
Painfo = Find("GSM2496939_S_Pombe.WT.MM.pA.read_counts.txt")


def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    print(genes)

    return  genes

controlgenesmp = Find('maxpas.csv')
controlfirstpas = Find('control_firstPA.csv')
control_ayo = Find('control_ayo.csv')

def Find(file):
    genes = pd.read_csv(file, delimiter="\t")
    chr1 = genes[(genes['AB325691'] == 'I')]
    chr2 = genes[(genes['AB325691'] == 'II')]
    chr3 = genes[(genes['AB325691'] == 'III')]
    
    
    theone = pd.concat([chr1, chr2, chr3])
    theone.rename(columns={'AB325691':'chro', '-':'coding_strand','8791':'pos', '1831':'count'}, inplace = True)


    return  theone

Painfo = Find("GSM2496939_S_Pombe.WT.MM.pA.read_counts.txt")



def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    print(genes)
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]


    return genesfor, genesrev, genes


controlfor, controlrev, newccontrol = Find('new_control.csv')

ggenes = ggenes.merge(Genes2_3, on="ID", how='left')
newccontrol = newccontrol.merge(Genes2_3, on="ID", how='left')
ggenes['length_UTR'] = ggenes['stop'] - ggenes['start_y']
mean_ggeensnsnenske = ggenes['length_UTR'].mean()

sns.swarmplot(ggenes, x = 'length_UTR', y= 'genotype', hue= 'genotype')

#genesfor = ggenes[(ggenes['coding_strand'] == 'forward')]
#genesrev = ggenes[(ggenes['coding_strand'] == 'reverse')]

ggenes.dropna(subset=['start_y'],inplace=True)
newccontrol.dropna(subset=['start_y'],inplace=True)

ggenes.loc[ggenes['coding_strand'] == 'forward', 'coding_strand'] = '+'
ggenes.loc[ggenes['coding_strand'] == 'reverse', 'coding_strand'] = '-'

def editforfunc(genes):
    genesfor = genes[(genes['coding_strand'] == 'forward')]
    genesrev = genes[(genes['coding_strand'] == 'reverse')]


    genesfor['Sbinpos'] = genesfor['start_y']/50
    genesfor['Sbinpos'] = genesfor['Sbinpos'].astype(int)
    genesfor['Sbinpos'] = genesfor['Sbinpos']*50 +25
    genesfor['Ebinpos'] = genesfor['stop']/50
    genesfor['Ebinpos'] = genesfor['Ebinpos'].astype(int)
    genesfor['Ebinpos'] = genesfor['Ebinpos']*50 +25


    genesrev['Sbinposr'] = genesrev['stop']/50
    genesrev['Sbinposr'] = genesrev['Sbinposr'].astype(int)
    genesrev['Sbinposr'] = genesrev['Sbinposr']*50 +25
    genesrev['Ebinposr'] = genesrev['start_y']/50
    genesrev['Ebinposr'] = genesrev['Ebinposr'].astype(int)
    genesrev['Ebinposr'] = genesrev['Ebinposr']*50 +25

    return genesfor, genesrev, genes

genesfor, genesrev, ggenes = editforfunc(ggenes)
controlfor, controlrev, controlcontrol = editforfunc(newccontrol)




gene1 = ggenes[(ggenes['chro_x'] == 'chr1')]
gene2 = ggenes[(ggenes['chro_x'] == 'chr2')]
gene3 = ggenes[(ggenes['chro_x'] == 'chr3')]


con1 = newccontrol[(newccontrol['chro_x'] == 'chr1')]
con2 = newccontrol[(newccontrol['chro_x'] == 'chr2')]
con3 = newccontrol[(newccontrol['chro_x'] == 'chr3')]


PA1 = Painfo[(Painfo['chro'] == 'I')]
PA2 = Painfo[(Painfo['chro'] == 'II')]
PA3 = Painfo[(Painfo['chro'] == 'III')]

#%%

PA1['id'] = range(1, len(PA1) + 1)
PA2['id'] = range(1, len(PA2) + 1)
PA3['id'] = range(1, len(PA3) + 1)

def assign_gene_id(row, gene_df):
    pos = row['pos']
    for _, gene_row in gene_df.iterrows():
        start_y = gene_row['start_y']
        stop = gene_row['stop']
        if start_y <= pos <= stop:
            return gene_row['ID']
    return None

#version two with coding strand 
def assign_gene_id(row, gene_df):
    pos = row['pos']
    coding_strand = row['coding_strand']
    
    for _, gene_row in gene_df.iterrows():
        start_y = gene_row['start_y']
        stop = gene_row['stop']
        if start_y <= pos <= stop and coding_strand == gene_row['coding_strand']:
            return gene_row['ID']
    
    return None
# Assuming your dataframes are named PA1 and gene1
PA1['ID'] = PA1.apply(lambda row: assign_gene_id(row, gene1), axis=1)
PA2['ID'] = PA2.apply(lambda row: assign_gene_id(row, gene2), axis=1)
PA3['ID'] = PA3.apply(lambda row: assign_gene_id(row, gene3), axis=1)

PA1.dropna(subset=['ID'], inplace=True)
PA2.dropna(subset=['ID'], inplace=True)
PA3.dropna(subset=['ID'], inplace=True)

################################
id_counts = PA1.groupby('ID').size().reset_index(name='count_pas')
id_counts2 = PA2.groupby('ID').size().reset_index(name='count_pas')
id_counts3 = PA3.groupby('ID').size().reset_index(name='count_pas')

gene1y = gene1.merge(id_counts, on='ID', how='left')
gene2 = gene2.merge(id_counts2, on='ID', how='left')
gene3 = gene3.merge(id_counts3, on='ID', how='left')

ayo = pd.concat([gene1y, gene2, gene3], ignore_index=True)
# Display the result
print(id_counts)


##########this finds the first PA site listed for the gene 
min_pos_indices = PA1.groupby('ID')['pos'].idxmin()
min_pos_indices2 = PA2.groupby('ID')['pos'].idxmin()
min_pos_indices3 = PA3.groupby('ID')['pos'].idxmin()

# Create a new DataFrame with only the rows corresponding to the minimum "pos" in each group
min_pos_df_1 = PA1.loc[min_pos_indices]
min_pos_df_2 = PA2.loc[min_pos_indices2]
min_pos_df_3 = PA3.loc[min_pos_indices3]

firstPA = pd.concat([min_pos_df_1, min_pos_df_2, min_pos_df_3], ignore_index=True)
firstPA = pd.merge(firstPA, ggenes[['ID', 'genotype', 'start_y', 'stop', 'strand', 'stalled_fork']], on='ID', how='left')
firstPA['relativepos'] = (firstPA['pos'] - firstPA['start_y'])/(firstPA['stop'] - firstPA['start_y'])
firstPA['length'] = (firstPA['stop'] - firstPA['start_y'])
sns.swarmplot(data=firstPA,x='genotype', y='length', hue='genotype')



custom_palette = ["steelblue", "orange", 'tomato', 'teal']  

sns.set_palette(custom_palette)



controlfirstpas['genotype'] = 'control'
firstPAy = pd.concat([firstPA, controlfirstpas], ignore_index=True)


firstPAy["genotype"] = pd.Categorical(firstPAy["genotype"],
                                                 ordered = True,
                                                 categories = ["sen1D", "dbl8D", "sen1dbl8DD_unique", 'control'])


sns.swarmplot(data=firstPAy,x='genotype', y='relativepos', hue='genotype')
plt.yticks([0,0.2,0.4,0.6,0.8,1], ['UTR start','0.2','0.4','0.6','0.8','UTR end'])

# Show the plot
plt.show()
ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')

##############this find the max count associated with a PA in the listed genes
max_count_indices = PA1.groupby('ID')['count'].idxmax()
max_count_indices2 = PA2.groupby('ID')['count'].idxmax()
max_count_indices3 = PA3.groupby('ID')['count'].idxmax()

result_PA1 = PA1.loc[max_count_indices]
result_PA2 = PA2.loc[max_count_indices2]
result_PA3 = PA3.loc[max_count_indices3]

maxcountPA = pd.concat([result_PA1, result_PA2, result_PA3], ignore_index=True)
maxcountPA = pd.merge(maxcountPA, ggenes[['ID', 'genotype', 'start_y', 'stop', 'strand', 'stalled_fork']], on='ID', how='left')
maxcountPA['relativepos'] = (maxcountPA['pos'] - maxcountPA['start_y'])/(maxcountPA['stop'] - maxcountPA['start_y'])
maxcountPA['length'] = (maxcountPA['stop'] - maxcountPA['start_y'])
sns.swarmplot(data=maxcountPA,x='genotype', y='relativepos', hue='genotype')
sns.kdeplot(data=maxcountPA, y='relativepos', hue='genotype')



controlgenesmp['genotype'] = 'control'
maxcountPAy = pd.concat([maxcountPA, controlgenesmp], ignore_index=True)


maxcountPAy["genotype"] = pd.Categorical(maxcountPAy["genotype"],
                                                 ordered = True,
                                                 categories = ["sen1D", "dbl8D", "sen1dbl8DD_unique", 'control'])


sns.swarmplot(data=maxcountPAy,x='genotype', y='relativepos', hue='genotype')
plt.yticks([0,0.2,0.4,0.6,0.8,1], ['UTR start','0.2','0.4','0.6','0.8','UTR end'])

ax1.tick_params(axis='both', labelsize=24)
#plotting the number of PAS in each gene subset 


control_ayo
control_ayo['genotype'] = 'control'
ayoy = pd.concat([ayo, control_ayo], ignore_index=True)


ayoy["genotype"] = pd.Categorical(ayoy["genotype"],
                                                 ordered = True,
                                                 categories = ["sen1D", "dbl8D", "sen1dbl8DD_unique", 'control'])


sns.boxplot(data=ayoy, x="count_pas", hue="genotype", fill=False, gap=.1)





maxcountPA_subset = maxcountPA[['pos', 'start_y', 'stop', 'genotype']].add_suffix('_max')

# Select columns from df2 with the desired suffix
firstPA_subset = firstPA[['pos']].add_suffix('_first')

# Concatenate the subsets into a new DataFrame
new_df = pd.concat([maxcountPA_subset, firstPA_subset], axis=1)

# Display the result
print(new_df)

new_df['relativeposmax'] = (new_df['pos_max'] - new_df['start_y_max'])/(new_df['stop_max'] - new_df['start_y_max'])

new_df['relativeposfirst'] = (new_df['pos_first'] - new_df['start_y_max'])/(new_df['stop_max'] - new_df['start_y_max'])





new_dfmelty = new_df.melt(id_vars=['pos_max', 'start_y_max', 'stop_max', 'pos_first', 'genotype_max'],
                               var_name= 'source',
                               value_vars= ['relativeposmax',
                                      'relativeposfirst'],
                               value_name='epb_')

sns.pointplot(data=new_dfmelty, x="source", y="epb_", hue='genotype_max')

sns.swarmplot(data=new_dfmelty, x="source", y="epb_")
sns.kdeplot(data=new_dfmelty, x="epb_", hue='source')



#%%

######this is from chatgtp 
ffeat = ffeat.rename(columns={'Systematic ID': 'ID'})
grouped_ffeat = ffeat.groupby('chro').apply(lambda x: x.sort_values(by='Start position')).reset_index(drop=True)

# Function to check status based on conditions
def check_status(row):
    if row['Start position'] > row['pos'] + 1000:
        return 'fail'
    else:
        return 'pass'

# Apply the function to each row in maxcountPA
maxcountPA['status'] = maxcountPA.apply(lambda row: check_status(row), axis=1)

# Display the result
print(maxcountPA)


def check_status(row, ffeat_row):
    if ffeat_row['Start position'] > row['pos'] + 1000:
        return 'fail'
    else:
        return 'pass'

# Apply the function to each row in maxcountPA
maxcountPA['status'] = maxcountPA.apply(lambda row: check_status(row, grouped_ffeat[
    (grouped_ffeat['chro'] == row['chro']) & (grouped_ffeat['ID'] == row['ID'])
].iloc[0]), axis=1)

# Display the result
print(maxcountPA)


def check_status(row, ffeat_df):
    current_pos = row['pos']
    
    # Get the rows from ffeat_df that match the condition
    condition = (ffeat_df['chro'] == row['chro']) & (ffeat_df['ID'] == row['ID'])
    matching_rows = ffeat_df[condition]
    
    # Check if matching rows exist and perform the comparison
    if not matching_rows.empty:
        next_start_pos = matching_rows['Start position'].shift(-1)
        if not next_start_pos.empty and current_pos > next_start_pos.iloc[0] + 1000:
            return 'fail'
    
    return 'pass'

# Apply the function to each row in maxcountPA
maxcountPA['status'] = maxcountPA.apply(lambda row: check_status(row, grouped_ffeat), axis=1)

# Display the result
print(maxcountPA)

def check_status(row, ffeat_df):
    current_pos = row['pos']
    
    # Get the rows from ffeat_df that match the condition
    condition = (ffeat_df['chro'] == row['chro']) & (ffeat_df['ID'] == row['ID'])
    matching_rows = ffeat_df[condition]
    
    # Print statements for troubleshooting
    print(f"Current row: {row}")
    print(f"Matching rows: {matching_rows}")
    
    # Check if matching rows exist and perform the comparison
    if not matching_rows.empty:
        next_start_pos = matching_rows['Start position'].shift(-1)
        print(f"Next start pos: {next_start_pos}")
        if not next_start_pos.empty and current_pos > next_start_pos.iloc[0] + 1000:
            return 'fail'
    
    return 'pass'

# Apply the function to each row in maxcountPA
maxcountPA['status'] = maxcountPA.apply(lambda row: check_status(row, grouped_ffeat), axis=1)

# Display the result
print(maxcountPA)


def check_status(row, ffeat_df):
    current_pos = row['pos']
    
    # Get the rows from ffeat_df that match the condition
    condition = (ffeat_df['chro'] == row['chro']) & (ffeat_df['ID'] == row['ID'])
    matching_rows = ffeat_df[condition]
    print(matching_rows)
    
    # Iterate through matching rows and check the condition
    for i in range(len(matching_rows) - 1):
        next_start_pos = matching_rows.iloc[i + 1]['Start position']
        
        # Print general message
        print(f'Checking row {i + 1} in the loop')
        
        if not pd.isna(next_start_pos):
            print(f'Next Start Position: {next_start_pos}')
            if current_pos < next_start_pos + 1000:
                return 'fail'
    
    return 'pass'

# Apply the function to each row in maxcountPA
maxcountPA['status'] = maxcountPA.apply(lambda row: check_status(row, grouped_ffeat), axis=1)

# Display the result
print(maxcountPA)

#%%

#now i need to write something to check if in maxcountPA, 
#if the next gene overlaps with PAS+1000
maxcountPA = maxcountPAy.copy()
maxcountPAy['pos_1000'] = maxcountPAy['pos'] +1000
sns.boxplot(maxcountPAy, y = 'length', x= 'genotype', hue= 'genotype', fill=False)

#quick check that everybody in macountPA is in ffeat
ffeat.rename(columns={'Systematic ID':'ID'}, inplace = True)

ffeat.loc[ffeat['Strand'] == 'forward', 'strand'] = '+'
ffeat.loc[ffeat['Strand'] == 'reverse', 'strand'] = '-'
ffeat.loc[ffeat['Chromosome'] == 'I', 'chro'] = 'chr1'
ffeat.loc[ffeat['Chromosome'] == 'II', 'chro'] = 'chr2'
ffeat.loc[ffeat['Chromosome'] == 'III', 'chro'] = 'chr3'



    maxcountPAy['pos_1000'] = maxcountPAy['pos_1000']/50
    maxcountPAy['pos_1000'] = maxcountPAy['pos_1000'].astype(int)
    maxcountPAy['pos_1000'] = maxcountPAy['pos_1000']*50 +25
    ffeat['Sbinpos'] = ffeat['Start position']/50
    ffeat['Sbinpos'] = ffeat['Sbinpos'].astype(int)
    ffeat['Sbinpos'] = ffeat['Sbinpos']*50 +25

ffeat1 = ffeat[(ffeat['chro'] == 'chr1')]
ffeat2 = ffeat[(ffeat['chro'] == 'chr2')]
ffeat3 = ffeat[(ffeat['chro'] == 'chr3')]

maxcountPA1 = maxcountPAy[(maxcountPAy['chro'] == 'chr1')]
maxcountPA2 = maxcountPAy[(maxcountPAy['chro'] == 'chr2')]
maxcountPA3 = maxcountPAy[(maxcountPAy['chro'] == 'chr3')]

#checkthisbitch = pd.merge(ffeat, maxcountPA[['ID', 'genotype']], on='ID', how='left')
#checkthisbitch.dropna(subset=['genotype'], inplace=True)
#it does, all good

def score(file,maxdf):
    #count=0
    file = file.sort_values(by='Sbinpos')
    maxdf = maxdf.sort_values(by='pos')
    
    for i in range(min(len(file), len(maxdf)) - 1):
        tempID = maxdf.iloc[i]["ID"]
        tempchro = maxdf.iloc[i]["chro"]
        

        tempsubset = file.loc[(file['ID']== tempID) & (file['chro'] == tempchro)]
        if i + 1 < len(file):
            next_row = file.iloc[i + 1]
            maxdf.loc[maxdf['ID'] == tempID, 'Next_Sbinpos'] = next_row['Sbinpos']

            
       # count += 1
        print(next_row)
       # print(count)
        
    return maxdf

#trytrytrytrytry = score(ffeat1, maxcountPA1)




def score(file, maxdf):
    file = file.sort_values(by='Sbinpos')
    maxdf = maxdf.sort_values(by='pos')

    for i in range(len(maxdf) - 1):
        tempID = maxdf.iloc[i]["ID"]
        tempchro = maxdf.iloc[i]["chro"]

        tempsubset = file.loc[(file['ID'] == tempID) & (file['chro'] == tempchro)]
        if not tempsubset.empty:
            next_row = file.iloc[file.index.get_loc(tempsubset.index[-1]) + 1]  # Find the next row after tempsubset
            maxdf.loc[maxdf['ID'] == tempID, 'Next_Sbinpos'] = next_row['Sbinpos']

    return maxdf

# Example usage
trytrytrytrytry1 = score(ffeat1, maxcountPA1)
trytrytrytrytry2 = score(ffeat2, maxcountPA2)
trytrytrytrytry3 = score(ffeat3, maxcountPA3)


notryingtodie = pd.concat([trytrytrytrytry1, trytrytrytrytry2, trytrytrytrytry3])

notryingtodie['distance'] = notryingtodie['Next_Sbinpos'] - notryingtodie['pos_1000']
notryingtodie = notryingtodie.sort_values(by='distance')
notryingtodie["genotype"] = pd.Categorical(notryingtodie["genotype"],
                                                 ordered = True,
                                                 categories = ["sen1D", "dbl8D", "sen1dbl8DD_unique", "control"])

sns.barplot(notryingtodie, x="ID", y="distance", hue="genotype", legend=False)


sns.swarmplot(notryingtodie, x="genotype", y="distance", hue="genotype", legend=False)
flights_wide = notryingtodie.pivot(columns="genotype", values="distance")
sns.barplot(flights_wide)
sns.stripplot(notryingtodie, x="genotype", y="distance", hue="genotype", legend=False)
sns.boxplot(notryingtodie, x="genotype", y="distance", hue="genotype", showfliers=False,legend=False, fill=False)


notryingtodie['relative_pos_log2'] = np.log2(notryingtodie['relativepos'])
notryingtodie['distance_toi_incoming_log2'] = np.log2(notryingtodie['distance'])

sns.lmplot(
    data=notryingtodie, x="relative_pos_log2", y="distance_toi_incoming_log2",
    col="genotype", height=3,
    facet_kws=dict(sharex=False, sharey=False),
)

sns.scatterplot(
    data=notryingtodie, x="relativepos", y="distance", hue="genotype", legend="full"
)


noover_notryingtodie = notryingtodie[notryingtodie['distance'] > 0]


firstPAfor = noover_notryingtodie[(noover_notryingtodie['coding_strand'] == '+')]
firstPArec = noover_notryingtodie[(noover_notryingtodie['coding_strand'] == '-')]

#%%
#let's try to fix the overlap problem take PAS +1000
#say, does this overlap with start of the next gene. 
#rank feat by start pos
#from ggenes, find it in feat df, and ask what is the next gene
#is pas +1000 > start [i+1]
#########
#########
########
#so now I have three dfs to play with:
    #the first is the ayo df which tells me how many PAS in each stall gene
    #the second tells me the first PAS
    #the third tells me the PAS with the max count in the UTR. 
    
    #I shall need to see how far they are from the UTR start stie 

firstPA.loc[firstPA['chro'] == 'I', 'chro'] = 'chr1'
firstPA.loc[firstPA['chro'] == 'II', 'chro'] = 'chr2'
firstPA.loc[firstPA['chro'] == 'III', 'chro'] = 'chr3'

firstPAsen = firstPA[(firstPA['genotype'] == 'sen1D')]
firstPAdbl = firstPA[(firstPA['genotype'] == 'dbl8D')]
firstPAds = firstPA[(firstPA['genotype'] == 'sen1dbl8DD_unique')]


def editforfunc(genesfor):

    genesfor['Sbinpos'] = genesfor['pos']/50
    genesfor['Sbinpos'] = genesfor['Sbinpos'].astype(int)
    genesfor['Sbinpos'] = genesfor['Sbinpos']*50 +25

    return genesfor

firstPA = editforfunc(firstPA)
maxcountPAy = editforfunc(maxcountPAy)

columns_to_drop = ['type_x']
maxcountPAy = maxcountPAy.drop(columns=columns_to_drop)

####so just remeber to hit the [13]
firstPAfor = firstPA[(firstPA['coding_strand'] == '+')]
firstPArec = firstPA[(firstPA['coding_strand'] == '-')]


maxcountPA.loc[maxcountPA['chro'] == 'I', 'chro'] = 'chr1'
maxcountPA.loc[maxcountPA['chro'] == 'II', 'chro'] = 'chr2'
maxcountPA.loc[maxcountPA['chro'] == 'III', 'chro'] = 'chr3'


firstPAfor = maxcountPA[(maxcountPA['strand'] == '+')]
firstPArec = maxcountPA[(maxcountPA['strand'] == '-')]

filtered_df_RF = maxcountPA[(maxcountPA['genotype'] == 'sen1dbl8DD_unique') & (maxcountPA['stalled_fork'] == 'rightward') & (maxcountPA['strand'] == '+')]
count = len(filtered_df_RR)
print(count)

filtered_df_LR = maxcountPA[(maxcountPA['genotype'] == 'sen1dbl8DD_unique') & (maxcountPA['stalled_fork'] == 'leftward') & (maxcountPA['strand'] == '-')]
count = len(filtered_df_LF)
print(count)

def REnd(genesfor, chr1, chr2, chr3, p, k, s):
    xxe =[]
    for ge in genesfor.itertuples():
        if ge[10] == s:
            if ge[11] == k:
                if ge[7] == p:
                    if ge[1] == 'chr1':
                        xe = chr1.index[chr1['pos'] == ge[14]].tolist()
                        xxe.append(xe) 
    
                    if ge[1] == 'chr2':
                        xe2 = chr2.index[chr2['pos'] == ge[14]].tolist()
                        xxe.append(xe2)
    
                    if ge[1] == 'chr3':
                        xe3 = wtchr3.index[wtchr3['pos'] == ge[14]].tolist()
                        xxe.append(xe3)
    return xxe

dsforxxeR = REnd(firstPAfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward', '+')
dsrevxxeR = REnd(firstPArec, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'rightward', '-')
dsforxxeL = REnd(firstPAfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward', '+')
dsrevxxeL = REnd(firstPArec, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique', 'leftward', '-')


#for wt , forward coding strand genes only 
#TSS
def Start(genesfor, chr1, chr2, chr3, p):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[7] == p:
            if g[1] == 'chr1':
                x = chr1.index[chr1['pos'] == g[14]].tolist()
                xx.append(x) 

            if g[1] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[14]].tolist()
                xx.append(x2)

            if g[1] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[14]].tolist()
                xx.append(x3)
                
    return xx

senforxx = Start(firstPAfor, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
senrevxx = Start(firstPArec, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxx = Start(firstPAfor, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
dbl8revxx = Start(firstPArec, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')

controlforxx = Start(firstPAfor, wtchr1, wtchr2, wtchr3, 'control')
controlrevxx = Start(firstPArec, wtchr1, wtchr2, wtchr3, 'control')

#dsforxxeL = Start(firstPAfor, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique')
#dsrevxxeL = Start(firstPArec, sen1chr1, sen1chr2, sen1chr3, 'sen1dbl8DD_unique')




def tssflank(list1,frame):

    yucky=[]
    for z in (list1):
        yy = frame.loc[z[0]-20:z[0]+20,'normalised'].tolist()
        yucky.append(yy)
    tssflanky = pd.DataFrame(yucky)
   # tssflanky = np.array(yucky)

    return tssflanky
senforflankw = tssflank(senforxx, wt)
senforflanks = tssflank(senforxx, sen1)
senforflankb = tssflank(senforxx, dbl8)
senforflankds = tssflank(senforxx, ds)

dbl8forflankw = tssflank(dbl8forxx, wt)
dbl8forflanks = tssflank(dbl8forxx, sen1)
dbl8forflankb = tssflank(dbl8forxx, dbl8)
dbl8forflankds = tssflank(dbl8forxx, ds)

dsforflankRw = tssflank(dsforxxeR, wt)
dsforflankRs = tssflank(dsforxxeR, sen1)
dsforflankRb = tssflank(dsforxxeR, dbl8)
dsforflankRds = tssflank(dsforxxeR, ds)

dsforflankLw = tssflank(dsforxxeL, wt)
dsforflankLs = tssflank(dsforxxeL, sen1)
dsforflankLb = tssflank(dsforxxeL, dbl8)
dsforflankLds = tssflank(dsforxxeL, ds)

senrevflankw = tssflank(senrevxx, wt)
senrevflanks = tssflank(senrevxx, sen1)
senrevflankb = tssflank(senrevxx, dbl8)
senrevflankds = tssflank(senrevxx, ds)


dbl8revflankw = tssflank(dbl8revxx, wt)
dbl8revflanks = tssflank(dbl8revxx, sen1)
dbl8revflankb = tssflank(dbl8revxx, dbl8)
dbl8revflankds = tssflank(dbl8revxx, ds)

dsrevflankRw = tssflank(dsrevxxeR, wt)
dsrevflankRs = tssflank(dsrevxxeR, sen1)
dsrevflankRb = tssflank(dsrevxxeR, dbl8)
dsrevflankRds = tssflank(dsrevxxeR, ds)


dsrevflanklLw = tssflank(dsrevxxeL, wt)
dsrevflanklLs = tssflank(dsrevxxeL, sen1)
dsrevflanklLb = tssflank(dsrevxxeL, dbl8)
dsrevflanklLds = tssflank(dsrevxxeL, ds)

controlflankforw = tssflank(controlforxx, wt)
controlflankfors = tssflank(controlforxx, sen1)
controlflankforb = tssflank(controlforxx, dbl8)
controlflankfords = tssflank(controlforxx, ds)


controlflankrevw = tssflank(controlrevxx, wt)
controlflankrevs = tssflank(controlrevxx, sen1)
controlflankrevb = tssflank(controlrevxx, dbl8)
controlflankrevds = tssflank(controlrevxx, ds)

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
senforflankw = Expand_genes(senforflankw)
senforflanks = Expand_genes(senforflanks)
senforflankb = Expand_genes(senforflankb)
senforflankds = Expand_genes(senforflankds)

dbl8forflankw = Expand_genes(dbl8forflankw)
dbl8forflanks = Expand_genes(dbl8forflanks)
dbl8forflankb = Expand_genes(dbl8forflankb)
dbl8forflankds = Expand_genes(dbl8forflankds)

dsforflankRw = Expand_genes(dsforflankRw)
dsforflankRs = Expand_genes(dsforflankRs)
dsforflankRb = Expand_genes(dsforflankRb)
dsforflankRds = Expand_genes(dsforflankRds)

dsforflankLw = Expand_genes(dsforflankLw)
dsforflankLs = Expand_genes(dsforflankLs)
dsforflankLb = Expand_genes(dsforflankLb)
dsforflankLds = Expand_genes(dsforflankLds)



senrevflankw = Expand_genes(senrevflankw)
senrevflanks = Expand_genes(senrevflanks)
senrevflankb = Expand_genes(senrevflankb)
senrevflankds = Expand_genes(senrevflankds)


dbl8revflankw = Expand_genes(dbl8revflankw)
dbl8revflanks = Expand_genes(dbl8revflanks)
dbl8revflankb = Expand_genes(dbl8revflankb)
dbl8revflankds = Expand_genes(dbl8revflankds)

dsrevflankRw = Expand_genes(dsrevflankRw)
dsrevflankRs = Expand_genes(dsrevflankRs)
dsrevflankRb = Expand_genes(dsrevflankRb)
dsrevflankRds = Expand_genes(dsrevflankRds)


dsrevflanklLw = Expand_genes(dsrevflanklLw)
dsrevflanklLs = Expand_genes(dsrevflanklLs)
dsrevflanklLb = Expand_genes(dsrevflanklLb)
dsrevflanklLds = Expand_genes(dsrevflanklLds)


controlflankforw = Expand_genes(controlflankforw)
controlflankfors = Expand_genes(controlflankfors)
controlflankforb = Expand_genes(controlflankforb)
controlflankfords = Expand_genes(controlflankfords)

controlflankrevw = Expand_genes(controlflankrevw)
controlflankrevs = Expand_genes(controlflankrevs)
controlflankrevb = Expand_genes(controlflankrevb)
controlflankrevds = Expand_genes(controlflankrevds)

def concati(rev, forward):
    forwards = np.stack(forward, axis=0 )
    revs = np.stack(rev, axis=0 )
    rev1 = revs[:, ::-1]
    
    new = np.concatenate([forwards,rev1])
    return new

####TSS only first okay 
#reverse LEFT, forward right 
HTw = concati(dsrevflanklLw, dsforflankRw)
HTs = concati(dsrevflanklLs, dsforflankRs)
HTb = concati(dsrevflanklLb, dsforflankRb)
HTds = concati(dsrevflanklLds,dsforflankRds)

#reverse right + forward leeft 
HHw = concati(dsrevflankRw, dsforflankLw)
HHs = concati(dsrevflankRs, dsforflankLs)
HHb = concati(dsrevflankRb, dsforflankLb)
HHds = concati(dsrevflankRds, dsforflankLds)

sa = concati(senrevflankw, senforflankw)
sas = concati(senrevflanks,senforflanks)
sab = concati(senrevflankb,senforflankb)
sads = concati(senrevflankds,senforflankds)

dbl8w = concati(dbl8revflankw, dbl8forflankw)
dbl8s = concati(dbl8revflanks, dbl8forflanks)
dbl8b = concati(dbl8revflankb, dbl8forflankb)
dbl8ds = concati(dbl8revflankds, dbl8forflankds)


caw = concati(controlflankrevw, controlflankforw)
cas = concati(controlflankrevs, controlflankfors)
cab = concati(controlflankrevb, controlflankforb)
cads = concati(controlflankrevds, controlflankfords)

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
    print(wow)
    wow[0] = wow[0].rolling(window = 50, center=True).mean()   
    
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
#swline = pileup(sew)
#ssline = pileup(ses)
#sbline = pileup(seb)
#sdsline = pileup(seds)


#control gene body
cwline = pileup(caw)
csline = pileup(cas)
cbline = pileup(cab)
cdsline = pileup(cads)



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








figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.tick_params(axis='both', labelsize=24)
#ax1.set_ylim([0, 140])
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


ax1.legend(loc='best')
ax1.set_xticks([0,500,1000])
ax1.set_xticklabels(['-1Kb', 'PAS', '+1Kb'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')



figgy, (ax1) =plt.subplots(1,1, sharey=True)
ax1.tick_params(axis='both', labelsize=24)

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
ax1.set_xticks([0,500,1000])
ax1.set_xticklabels(['-1Kb', 'PAS', '+1Kb'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


figgy, (ax1) =plt.subplots(1, sharey=True)
#ds left
ax1.tick_params(axis='both', labelsize=24)

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
ax1.legend(loc='best')
ax1.set_xticks([0,500,1000])
ax1.set_xticklabels(['-1Kb', 'PAS', '+1Kb'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


figgy, (ax1) =plt.subplots(1, sharey=True)
#ds right
ax1.tick_params(axis='both', labelsize=24)

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
ax1.legend(loc='best')
ax1.set_xticks([0,500,1000])
ax1.set_xticklabels(['-1Kb', 'PAS', '+1Kb'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')


#%%

def Chromosome_plot (data1, data2, data3, data4, featurex, centro, telo, genee, c, con):
    ff, (ax1, ax2, axb, ax3, ax5) = plt.subplots(5,1, sharex=True, sharey=True)
    #ax1.set_title('WT')
    ax1.plot(data1['pos'], data1['norm'], color ='black', alpha=0.8)
    ax2.plot(data2['pos'], data2['norm'], color ='steelblue', alpha=0.8)
    axb.plot(data3['pos'], data3['norm'], color ='orange', alpha=0.8)
    ax3.plot(data4['pos'], data4['norm'], color ='tomato', alpha=0.8)
   # ax1.set_ylim(0,200)
    ax1.set_ylabel('rightward forks (e)')
    

    ax5.set_xlabel('Chromosome position')

                  
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

    for co in con.itertuples(index=False, name=None):
        
        if co[4] == 'forward':
            ax5.axvspan(co[2],co[3],0.5,0.8,color="blue",alpha=0.3)
        elif co[4] == 'reverse':
            ax5.axvspan(co[2],co[3],0.3,0.5,color="blue",alpha=0.3)
    return ff

chromosome1 = Chromosome_plot(wtchr1, sen1chr1, dbl8chr1, dschr1, gene1, centro1, telo1, feat1, 'I', con1)
chromosome2 = Chromosome_plot(ewtchrII, esen1chr2, edbl8chr2, edschr2, gene2, centro2, telo2, feat2, 'II', con2)
chromosome3 = Chromosome_plot(ewtchrIII, esen1chr3, edbl8chr3, edschr3, gene3, centro3, telo3, feat3, 'III', con3)



#%%
def score(file, frame1, frame2, frame3):
    for i in range(len(file)):
        #this is replicon start and stop
        tempstart = file.iloc[i]["Start position"]
        tempend = file.iloc[i]['End position']
        tempchro = file.iloc[i]["chro"]
       # rna = file.iloc[i]['rna']

        tempsubset = frame1.loc[(frame1['pos'] >= tempstart) & (frame1['pos'] <= tempend) & (frame1['chro'] == tempchro)]
        tempsubsett = frame2.loc[(frame2['pos'] >= tempstart) & (frame2['pos'] <= tempend) & (frame2['chro'] == tempchro)]
        tempsubsetm = frame3.loc[(frame3['pos'] >= tempstart) & (frame3['pos'] <= tempend) & (frame3['chro'] == tempchro)]
        

        file.loc[file.index[i], 'thing'] = tempsubset['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_sen'] = tempsubset['d_usage_norm'].sum() / len(tempsubset)
        
        file.loc[file.index[i], 'thingt'] = tempsubsett['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_dbl'] = tempsubsett['d_usage_norm'].sum() / len(tempsubsett)
        
        file.loc[file.index[i], 'thingm'] = tempsubsetm['d_usage_norm'].sum()
        file.loc[file.index[i], 'epb_ds'] = tempsubsetm['d_usage_norm'].sum() / len(tempsubsetm)

        # Add log2 of epb and epbt as new columns
        file.loc[file.index[i], 'log2_epb_sen'] = np.log2(file.loc[file.index[i], 'epb_sen'])
        file.loc[file.index[i], 'log2_epb_dbl'] = np.log2(file.loc[file.index[i], 'epb_dbl'])
        file.loc[file.index[i], 'log2_epb_ds'] = np.log2(file.loc[file.index[i], 'epb_ds'])

        # Calculate the length of the region and add log2_length column
   #     length = tempend - tempstart
    #    file.loc[file.index[i], 'log2_length'] = np.log2(length)
     #   file.loc[file.index[i], 'log2_rna'] = np.log2(rna)

    return file

allthestalls_extra = score(deltaforyourplot, sencut, dblcut, dscut)


#%%


def RStart(genesfor, chr1, chr2, chr3, p, k):
    xx =[]
    
    for g in genesfor.itertuples():
        if g[7] == k:
            if g[8] == p:
                if g[5] == 'chr1':
                    x = chr1.index[chr1['pos'] == g[14]].tolist()
                    xx.append(x) 

                if g[5] == 'chr2':
                    x2 = chr2.index[chr2['pos'] == g[14]].tolist()
                    xx.append(x2)

                if g[5] == 'chr3':
                    x3 = chr3.index[chr3['pos'] == g[14]].tolist()
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
                    xe = chr1.index[chr1['pos'] == ge[15]].tolist()
                    xxe.append(xe) 

                if ge[5] == 'chr2':
                    xe2 = chr2.index[chr2['pos'] == ge[15]].tolist()
                    xxe.append(xe2)

                if ge[5] == 'chr3':
                    xe3 = wtchr3.index[wtchr3['pos'] == ge[15]].tolist()
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
                x = chr1.index[chr1['pos'] == g[14]].tolist()
                xx.append(x) 

            if g[5] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[14]].tolist()
                xx.append(x2)

            if g[5] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[14]].tolist()
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
                x = chr1.index[chr1['pos'] == g[13]].tolist()
                xx.append(x) 

            if g[2] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[13]].tolist()
                xx.append(x2)

            if g[2] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[13]].tolist()
                xx.append(x3)
                
    return xx
controlforxx = ControlStart(controlfor, wtchr1, wtchr2, wtchr3)
controlrevxx = ControlStart(controlrev, wtchr1, wtchr2, wtchr3)


def ControlEnd(genesfor, chr1, chr2, chr3):
    xx =[]
    
    for g in genesfor.itertuples():
            if g[2] == 'chr1':
                x = chr1.index[chr1['pos'] == g[14]].tolist()
                xx.append(x) 

            if g[2] == 'chr2':
                x2 = chr2.index[chr2['pos'] == g[14]].tolist()
                xx.append(x2)

            if g[2] == 'chr3':
                x3 = chr3.index[chr3['pos'] == g[14]].tolist()
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
                xe = chr1.index[chr1['pos'] == ge[15]].tolist()
                xxe.append(xe) 

            if ge[5] == 'chr2':
                xe2 = chr2.index[chr2['pos'] == ge[15]].tolist()
                xxe.append(xe2)

            if ge[5] == 'chr3':
                xe3 = wtchr3.index[wtchr3['pos'] == ge[15]].tolist()
                xxe.append(xe3)
    return xxe

senforxxe = End(genesfor, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
senrevxxe = End(genesrev, sen1chr1, sen1chr2, sen1chr3, 'sen1D')
dbl8forxxe = End(genesfor, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
dbl8revxxe = End(genesrev, dbl8chr1, dbl8chr2, dbl8chr3, 'dbl8D')
#dsforxxe = End(genesfor, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')
#dsrevxxe = End(genesrev, dschr1, dschr2, dschr3, 'sen1dbl8DD_unique')


#selection of random 200 genes 
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
selforxx = selectionlStart(selectionf, wtchr1, wtchr2, wtchr3)
selrevxx = selectionlStart(selectionr, wtchr1, wtchr2, wtchr3)


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
selforxxe = selectionEnd(selectionf, wtchr1, wtchr2, wtchr3)
selrevxxe = selectionEnd(selectionr, wtchr1, wtchr2, wtchr3)


def Gene_bits(list1,list2,frame):

    yucky=[]
    for z, ze in zip(list1, list2):
        yy = frame.loc[z[0]:ze[0],'normalised'].tolist()
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
        yyr = frame.loc[ye[0]:y[0],'normalised'].tolist()
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

def tssflank(list1,frame):

    yucky=[]
    for z in (list1):
        yy = frame.loc[z[0]:z[0]+20,'norm'].tolist()
        yucky.append(yy)
    #tssflanky = pd.DataFrame(yucky)
    tssflanky = np.array(yucky)

    return tssflanky
senforflankw = tssflank(senforxx, wt)
senforflanks = tssflank(senforxx, sen1)
senforflankb = tssflank(senforxx, dbl8)
senforflankds = tssflank(senforxx, ds)

senrevflankw = tssflank(senrevxx, wt)
senrevflanks = tssflank(senrevxx, sen1)
senrevflankb = tssflank(senrevxx, dbl8)
senrevflankds = tssflank(senrevxx, ds)


dbl8forflankw = tssflank(dbl8forxx, wt)
dbl8forflanks = tssflank(dbl8forxx, sen1)
dbl8forflankb = tssflank(dbl8forxx, dbl8)
dbl8forflankds = tssflank(dbl8forxx, ds)

dbl8revflankw = tssflank(dbl8revxx, wt)
dbl8revflanks = tssflank(dbl8revxx, sen1)
dbl8revflankb = tssflank(dbl8revxx, dbl8)
dbl8revflankds = tssflank(dbl8revxx, ds)


dsforflankRw = tssflank(dsforxxR, wt)
dsforflankRs = tssflank(dsforxxR, sen1)
dsforflankRb = tssflank(dsforxxR, dbl8)
dsforflankRds = tssflank(dsforxxR, ds)

dsrevflankRw = tssflank(dsrevxxR, wt)
dsrevflankRs = tssflank(dsrevxxR, sen1)
dsrevflankRb = tssflank(dsrevxxR, dbl8)
dsrevflankRds = tssflank(dsrevxxR, ds)

dsforflankLw = tssflank(dsforxxL, wt)
dsforflankLs = tssflank(dsforxxL, sen1)
dsforflankLb = tssflank(dsforxxL, dbl8)
dsforflankLds = tssflank(dsforxxL, ds)

dsrevflanklLw = tssflank(dsrevxxL, wt)
dsrevflanklLs = tssflank(dsrevxxL, sen1)
dsrevflanklLb = tssflank(dsrevxxL, dbl8)
dsrevflanklLds = tssflank(dsrevxxL, ds)

controlflankforw = tssflank(controlforxx, wt)
controlflankfors = tssflank(controlforxx, sen1)
controlflankforb = tssflank(controlforxx, dbl8)
controlflankfords = tssflank(controlforxx, ds)

controlflankrevw = tssflank(controlrevxx, wt)
controlflankrevs = tssflank(controlrevxx, sen1)
controlflankrevb = tssflank(controlrevxx, dbl8)
controlflankrevds = tssflank(controlrevxx, ds)





def tesflank(list1,frame):

    yucky=[]
    for z in (list1):
        yy = frame.loc[z[0]-20:z[0],'norm'].tolist()
        yucky.append(yy)
    #tssflanky = pd.DataFrame(yucky)
    tssflanky = np.array(yucky)

    return tssflanky
esenforflankw = tesflank(senforxxe, wt)
esenforflanks = tesflank(senforxxe, sen1)
esenforflankb = tesflank(senforxxe, dbl8)
esenforflankds = tesflank(senforxxe, ds)

esenrevflankw = tesflank(senrevxxe, wt)
esenrevflanks = tesflank(senrevxxe, sen1)
esenrevflankb = tesflank(senrevxxe, dbl8)
esenrevflankds = tesflank(senrevxxe, ds)

edbl8forflankw = tesflank(dbl8forxxe, wt)
edbl8forflanks = tesflank(dbl8forxxe, sen1)
edbl8forflankb = tesflank(dbl8forxxe, dbl8)
edbl8forflankds = tesflank(dbl8forxxe, ds)

edbl8revflankw = tesflank(dbl8revxxe, wt)
edbl8revflanks = tesflank(dbl8revxxe, sen1)
edbl8revflankb = tesflank(dbl8revxxe, dbl8)
edbl8revflankds = tesflank(dbl8revxxe, ds)


edsforflankRw = tesflank(dsforxxeR, wt)
edsforflankRs = tesflank(dsforxxeR, sen1)
edsforflankRb = tesflank(dsforxxeR, dbl8)
edsforflankRds = tesflank(dsforxxeR, ds)

edsrevflankRw = tesflank(dsrevxxeR, wt)
edsrevflankRs = tesflank(dsrevxxeR, sen1)
edsrevflankRb = tesflank(dsrevxxeR, dbl8)
edsrevflankRds = tesflank(dsrevxxeR, ds)

edsforflankLw = tesflank(dsforxxeL, wt)
edsforflankLs = tesflank(dsforxxeL, sen1)
edsforflankLb = tesflank(dsforxxeL, dbl8)
edsforflankLds = tesflank(dsforxxeL, ds)

edsrevflanklLw = tesflank(dsrevxxeL, wt)
edsrevflanklLs = tesflank(dsrevxxeL, sen1)
edsrevflanklLb = tesflank(dsrevxxeL, dbl8)
edsrevflanklLds = tesflank(dsrevxxeL, ds)

econtrolflankforw = tesflank(controlforxxe, wt)
econtrolflankfors = tesflank(controlforxxe, sen1)
econtrolflankforb = tesflank(controlforxxe, dbl8)
econtrolflankfords = tesflank(controlforxxe, ds)

econtrolflankrevw = tesflank(controlrevxxe, wt)
econtrolflankrevs = tesflank(controlrevxxe, sen1)
econtrolflankrevb = tesflank(controlrevxxe, dbl8)
econtrolflankrevds = tesflank(controlrevxxe, ds)




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

'''
x = pd.DataFrame(sa)
m = x[0].mean()

s = x[0].std() 
dof = len(x[0])-1 
confidence = 0.95
t_crit = np.abs(t.ppf((1-confidence)/2,dof))
print(t_crit)
(m-s*t_crit/np.sqrt(len(x)), m+s*t_crit/np.sqrt(len(x))) 
'''

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
    print(wow)
    wow[0] = wow[0].rolling(window = 50, center=True).mean()   
    
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

from matplotlib import lines
print(lines.lineStyles.keys())



figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.tick_params(axis='both', labelsize=24)
#ax1.set_ylim([0, 140])
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
#ax1.set_xticklabels(['TES','+1kb'])
ax1.set_xticklabels(['UTR start','UTR END'])
#ax1.set_ylim([-4, 124])
#ax1.set_ylabel('normalised reads')
ax1.legend(loc='best')

ax1.set_xticks([0,1000])
ax1.set_xticklabels(['TSS','TES'])
ax1.legend(loc='best')
ax1.set_ylabel('normalised reads')
#ax1.fill_between(x, y1a, y2a, facecolor='mediumaquamarine', alpha =0.2)
#y2 = linerev
#ax2.fill_between(x, y2, facecolor='green', alpha =0.3)
#filkggy, (ax1) =plt.subplots(1,1, sharey=True)

figgy, (ax1) =plt.subplots(1,1, sharey=True)
ax1.tick_params(axis='both', labelsize=24)
#ax1.set_ylim([0, 140])
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
#ax1.set_xticklabels(['TES','+1kb'])
ax1.set_xticklabels(['-1Kb','TSS'])
#ax1.set_ylim([-4, 124])
ax1.legend(loc='best')
#ax1.set_ylabel('normalised reads')


figgy, (ax1) =plt.subplots(1, sharey=True)
#ds left
ax1.tick_params(axis='both', labelsize=24)
#ax1.set_ylim([0, 140])
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
#ax1.set_xticklabels(['TES','+1kb'])
ax1.set_xticklabels(['-1Kb','TSS'])
ax1.legend(loc='best')
#ax1.set_ylabel('normalised reads')

figgy, (ax1) =plt.subplots(1, sharey=True)
#ds right
ax1.tick_params(axis='both', labelsize=24)
#ax1.set_ylim([0, 140])
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
#ax1.set_xticklabels(['TES','+1kb'])
ax1.set_xticklabels(['-1Kb','TSS'])
ax1.legend(loc='best')
#ax1.set_ylabel('normalised reads')
#ax1.set_ylim([-4, 124])


figgy, (ax1) =plt.subplots(1, sharey=True)
ax1.tick_params(axis='both', labelsize=24)
#ax1.set_ylim([0, 140])
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
ax1.tick_params(axis='both', labelsize=24)
#ax1.set_ylim([0, 140])

ax1.set_xticks([0,1000])
#ax1.set_xticklabels(['TES','+1kb'])
ax1.set_xticklabels(['-1Kb','TSS'])
ax1.legend(loc='best')



#%%%

#DBL8 heatmap !!!!!
fxb, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fxb.tight_layout()
cbar_ax = fxb.add_axes([.91, .3, .03, .4])
sns.heatmap(dbl8w, cmap = 'coolwarm', ax=ax1 ,vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax3,vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax5, vmin= 0, vmax=100,cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=ax7,vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,500,1000])
ax7.set_xticklabels(['-1kb', 'PAS', '+1kb'])


fx, ((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fx.tight_layout()
cbar_ax = fx.add_axes([.91, .3, .03, .4])
sns.heatmap(sa, cmap = 'coolwarm', ax=ax1, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('Forward Strand')
ax1.set_ylabel('WT')
sns.heatmap(sas, cmap = 'coolwarm', ax=ax3, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(sab, cmap = 'coolwarm', ax=ax5, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(sads, cmap = 'coolwarm', ax=ax7, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')

ax7.set_xticks([0,500,1000])
ax7.set_xticklabels(['-1kb', 'PAS', '+1kb'])


fuckkkkkk, (ax2, ax4, ax6, ax8) = plt.subplots(4,1)
fuckkkkkk.tight_layout()
cbar_ax = fuckkkkkk.add_axes([.91, .3, .03, .4])
sns.heatmap(HTw, cmap = 'coolwarm', ax=ax2,vmin= 0, vmax=100,  cbar=True, cbar_ax=cbar_ax)
ax2.set_ylabel('WT')
ax2.axes.yaxis.set_ticks([])
ax2.set_title('S∆D∆ co-directional collisions')
ax2.axes.xaxis.set_ticks([])
sns.heatmap(HTs, cmap = 'coolwarm', ax=ax4, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax4.axes.yaxis.set_ticks([])
ax4.set_ylabel('Sen1∆')
ax4.axes.xaxis.set_ticks([])
sns.heatmap(HTb, cmap = 'coolwarm', ax=ax6, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax6.axes.yaxis.set_ticks([])
ax6.set_ylabel('Dbl8∆')  
ax6.axes.xaxis.set_ticks([])
sns.heatmap(HTds, cmap = 'coolwarm', ax=ax8, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax8.axes.yaxis.set_ticks([])
ax8.set_ylabel('Sen1∆Dbl8∆')

ax8.set_xticks([0,500,1000])
ax8.set_xticklabels(['-1kb', 'PAS', '+1kb'])


fxdsl,((ax1, ax3, ax5,ax7)) = plt.subplots(4,1)  
fxdsl.tight_layout()
cbar_ax = fxdsl.add_axes([.91, .3, .03, .4])
sns.heatmap(HHw, cmap = 'coolwarm', ax=ax1, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax1.axes.yaxis.set_ticks([])
ax1.axes.xaxis.set_ticks([])
ax1.set_title('S∆D∆ Head-Head collisions')
ax1.set_ylabel('WT')
sns.heatmap(HHs, cmap = 'coolwarm', ax=ax3, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax3.axes.yaxis.set_ticks([])
ax3.set_ylabel('Sen1∆')
ax3.axes.xaxis.set_ticks([])
sns.heatmap(HHb, cmap = 'coolwarm', ax=ax5, vmin= 0, vmax=100, cbar=True, cbar_ax=cbar_ax)
ax5.axes.yaxis.set_ticks([])
ax5.set_ylabel('Dbl8∆')
ax5.axes.xaxis.set_ticks([])
sns.heatmap(HHds, cmap = 'coolwarm', ax=ax7, vmin= 0,vmax=100, cbar=True, cbar_ax=cbar_ax)
ax7.axes.yaxis.set_ticks([])
ax7.set_ylabel('Sen1∆Dbl8∆')
ax8.set_xticks([0,500,1000])
ax8.set_xticklabels(['-1kb', 'PAS', '+1kb'])



    fig = plt.figure(constrained_layout=True) #create a figure
    gs = fig.add_gridspec(177, 36)  #sets the figure width and height (can only refer to these as integers)




# Now add as many subplots using the gridspace (gs) to position them precisely
# the first position is the height (ie here from 1 to 3 on the grid)
# the second position is the width (i.e here from 1 to 10 on the grid)
    ax1 = fig.add_subplot(gs[10:19, 5:10])
    sns.heatmap(sa, cmap = 'coolwarm', ax=ax1, cbar=False, vmin= 0, vmax=100)
    ax1.set_xticks([])
    ax1.set_yticks([])

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis

    ax2 = fig.add_subplot(gs[10:19, 12:17]) 
    sns.heatmap(sas, cmap = 'coolwarm', ax=ax2,  cbar=False,vmin= 0, vmax=100)
    ax2.set_xticks([])
    ax2.set_yticks([])

   # ax2.set_title('mutant right forks',loc = 'left', pad = 0) 

    ax3 = fig.add_subplot(gs[10:19, 19:24]) 
    sns.heatmap(sab, cmap = 'coolwarm', ax=ax3, cbar=False,vmin= 0, vmax=100)
    ax3.set_xticks([])
    ax3.set_yticks([])
    
    axa = fig.add_subplot(gs[10:19, 26:31]) 
    sns.heatmap(sads, cmap = 'coolwarm', ax=axa, cbar=False,vmin= 0, vmax=100)
    axa.set_xticks([])
    axa.set_yticks([])

   # ax3.set_title('difference right forks',loc = 'left', pad = 0) 
   
   
   
   
    
    ax4 = fig.add_subplot(gs[24:52, 5:10])
    sns.heatmap(dbl8w, cmap = 'coolwarm', ax=ax4,  cbar=False,vmin= 0, vmax=100,)
    ax4.set_xticks([])
    ax4.set_yticks([])
    
    ax5 = fig.add_subplot(gs[24:52, 12:17])
    sns.heatmap(dbl8s, cmap = 'coolwarm', ax=ax5,  cbar=False,vmin= 0, vmax=100,)
    ax5.set_xticks([])
    ax5.set_yticks([])
    
    ax6 = fig.add_subplot(gs[24:52, 19:24])
    sns.heatmap(dbl8b, cmap = 'coolwarm', ax=ax6,  cbar=False,vmin= 0, vmax=100,)
    ax6.set_xticks([])
    ax6.set_yticks([])
    
    axb = fig.add_subplot(gs[24:52, 26:31])
    sns.heatmap(dbl8ds, cmap = 'coolwarm', ax=axb,  cbar=False,vmin= 0, vmax=100,)
    axb.set_xticks([])
    axb.set_yticks([])    
    
    
    
    ax7 = fig.add_subplot(gs[57:94, 5:10])
    sns.heatmap(HTw, cmap = 'coolwarm', ax=ax7, cbar=False, vmin= 0, vmax=100,)
    #ax7.set_xticks([0,11, 61, 72])
    #ax7.set_xticklabels(['-', 's', 'e', '+'])
    ax7.set_yticks([])
    ax7.set_xticks([])
    
    
    ax8 = fig.add_subplot(gs[57:94, 12:17]) 
    sns.heatmap(HTs, cmap = 'coolwarm', ax=ax8,cbar=False,vmin= 0, vmax=100)
   # ax8.set_xticks([0,11, 61, 72])
  #  ax8.set_xticklabels(['-', 's', 'e', '+'])
    ax8.set_yticks([])
    ax8.set_xticks([])
    
    ax9 = fig.add_subplot(gs[57:94, 19:24]) 
    sns.heatmap(HTb, cmap = 'coolwarm', ax=ax9,cbar=False, vmin= 0, vmax=100,)
   # ax9.set_xticks([0,11, 61, 72])
  #  ax9.set_xticklabels(['-', 's', 'e', '+'])
    ax9.set_yticks([])
    ax9.set_xticks([])

    axc = fig.add_subplot(gs[57:94, 26:31]) 
    sns.heatmap(HTds, cmap = 'coolwarm', ax=axc,cbar=False, vmin= 0, vmax=100,)
   # ax9.set_xticks([0,11, 61, 72])
   # ax9.set_xticklabels(['-', 's', 'e', '+'])
    axc.set_yticks([])   
    axc.set_xticks([]) 


    
    ax7a = fig.add_subplot(gs[99:167, 5:10])
    sns.heatmap(HHw, cmap = 'coolwarm', ax=ax7a, cbar=False, vmin= 0, vmax=100,)
    ax7a.set_xticks([0,500, 1000])
    ax7a.set_xticklabels(['-', 'PAS', '+'])
    ax7a.set_yticks([])
    
    
    ax8a = fig.add_subplot(gs[99:167, 12:17]) 
    sns.heatmap(HHs, cmap = 'coolwarm', ax=ax8a,cbar=False,vmin= 0, vmax=100)
    ax8a.set_xticks([0,500, 1000])
    ax8a.set_xticklabels(['-', 'PAS', '+'])
    ax8a.set_yticks([])
    
    ax9a = fig.add_subplot(gs[99:167, 19:24]) 
    sns.heatmap(HHb, cmap = 'coolwarm', ax=ax9a,cbar=False, vmin= 0, vmax=100,)
    ax9a.set_xticks([0,500, 1000])
    ax9a.set_xticklabels(['-', 'PAS', '+'])
    ax9a.set_yticks([])

    axca = fig.add_subplot(gs[99:167, 26:31]) 
    sns.heatmap(HHds, cmap = 'coolwarm', ax=axca,cbar=False, vmin= 0, vmax=100,)
    axca.set_xticks([0,500, 1000])
    axca.set_xticklabels(['-', 'PAS', '+'])
    axca.set_yticks([])    
 
    
 
    
    fig = plt.figure(constrained_layout=True) #create a figure
    gs = fig.add_gridspec(207, 36)  #sets the figure width and height (can only refer to these as integers)




# Now add as many subplots using the gridspace (gs) to position them precisely
# the first position is the height (ie here from 1 to 3 on the grid)
# the second position is the width (i.e here from 1 to 10 on the grid)
    ax1 = fig.add_subplot(gs[10:197, 5:10])
    sns.heatmap(caw, cmap = 'coolwarm', ax=ax1, cbar=False, vmin= 0, vmax=100)

    ax1.set_xticks([0,500, 1000])
    ax1.set_xticklabels(['-', 'PAS', '+'])
    ax1.set_yticks([])  

   # ax1.set_title('wild type right forks',loc = 'left', pad = 0) # loc(ation) can be left,right,c centre. pad(ding) is space between label and axis

    ax2 = fig.add_subplot(gs[10:197, 12:17]) 
    sns.heatmap(cas, cmap = 'coolwarm', ax=ax2,  cbar=False,vmin= 0, vmax=100)
    ax2.set_xticks([0,500, 1000])
    ax2.set_xticklabels(['-', 'PAS', '+'])
    ax2.set_yticks([])  

   # ax2.set_title('mutant right forks',loc = 'left', pad = 0) 

    ax3 = fig.add_subplot(gs[10:197, 19:24]) 
    sns.heatmap(cab, cmap = 'coolwarm', ax=ax3, cbar=False,vmin= 0, vmax=100)
    ax3.set_xticks([0,500, 1000])
    ax3.set_xticklabels(['-', 'PAS', '+'])
    ax3.set_yticks([])  
    
    axa = fig.add_subplot(gs[10:197, 26:31]) 
    sns.heatmap(cads, cmap = 'coolwarm', ax=axa, cbar=False,vmin= 0, vmax=100)
    axa.set_xticks([0,500, 1000])
    axa.set_xticklabels(['-', 'PAS', '+'])
    axa.set_yticks([])  

   # ax3.set_title('difference right forks',loc = 'left', pad = 0) 
   
   

