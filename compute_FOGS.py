import pandas as pd
import sys, os
import seaborn as sns
import matplotlib.pyplot as plt

pd.set_option("display.precision", 15)



### INPUT
clusterFile = sys.argv[1].strip()
SNPdistance = sys.argv[2].strip()



### MAIN
dClones = pd.read_csv(clusterFile, sep='\t', dtype=str) # Read file with clones


### Load all the indexes
FOGS = pd.read_csv("GainLossEvents.csv", delim_whitespace=True, names = ['genomeA', 'genomeB', 'FOGS'], dtype=str)

SNP = pd.read_csv(SNPdistance, sep='\t', names = ['genomeA', 'genomeB', 'SNPdist'], dtype=str)
SNP['SNPdist'] = SNP['SNPdist'].replace(0, 1)


FOGS['genomeA'] = FOGS['genomeA'].str.replace(r'.gff.fna','')
FOGS['genomeB'] = FOGS['genomeB'].str.replace(r'.gff.fna','')


assignclone_l = FOGS.merge(dClones, left_on='genomeA', right_on='genome', how='left')
assignclone = assignclone_l.merge(dClones, left_on='genomeB', right_on='genome', how='left')

only_intra = assignclone.query('Clone_x == Clone_y')
only_intra.drop(columns=['genome_x', 'genome_y', 'Clone_y'], inplace=True)


addSNP = only_intra.merge(SNP, on=['genomeA','genomeB'], how='left')
addSNP['FOGS'] = addSNP['FOGS'].astype(float)
addSNP['SNPdist'] = addSNP['SNPdist'].astype(float)
addSNP['SNPdist'] = addSNP['SNPdist'].replace(0.0, 1)


addSNP['wFOGS'] = addSNP['FOGS']/addSNP['SNPdist']
final_merge = addSNP[['Clone_x', 'wFOGS']]
final_merge.to_csv("FOGS.tsv", sep='\t', index=False)
print(final_merge)


wFOGSmean = addSNP.groupby('Clone_x')['wFOGS'].mean().reset_index()

#Violin plot on clones or boxplot or boxenplot
sns.set(rc={'figure.figsize':(15,8)}, font_scale=1, style='whitegrid')
sorted_FOGS = wFOGSmean.sort_values(by=['wFOGS'], ascending=False)
sorted_list_Clones_FOGS = list(sorted_FOGS['Clone_x'])
sns.boxplot(data=addSNP, x='Clone_x', y='wFOGS', color='#3295a8', showfliers = False, order=sorted_list_Clones_FOGS, showmeans=True)
plt.set_xlabel("", fontsize=1)
plt.set_ylabel("FOGS", fontsize=10)
plt.show()
plt.savefig("FOGS_boxplots.svg")
