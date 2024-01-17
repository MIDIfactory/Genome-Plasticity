import os, sys, glob
import pandas as pd
import csv
from itertools import combinations
import itertools
import multiprocessing as mp
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components



inFile = sys.argv[1].strip #File with clustered orthology groups


#pandas settings
pd.set_option("display.precision", 32)


#Orthologous groups: keys protein_name: value_group
orthologous_groups={}
with open(inFile, "r") as clustered_proteins:
	c=clustered_proteins.read()
	cluster_content=c[:-1] #remove empty last line
	clustered=cluster_content.split("\n")
	for group in clustered:
		name=group.split(": ")[0]
		proteins=group.strip(": ").split(": ")[1]
		proteins_list=proteins.split("\t")
		for protein in proteins_list:
			p=str(protein)
			orthologous_groups.update({p:name})
all_orthologous_groups = list(set(orthologous_groups.values()))


#read into every prodigal output (or prokka)
#os.chdir("fixed_input_files")
list_roary_files=glob.glob("*.gff")

def GFF_reader(GFF_file):
	with open(GFF_file) as f:
		seqname_list = []
		start_list = []
		end_list = []
		ID_list = []
		for line in f:
			chomped_line = line.rstrip(os.linesep)
			if chomped_line.startswith('##FASTA'):
				break
			if chomped_line.startswith('#'):
				continue
			else:
				elems = chomped_line.split('\t')
				try:
					seqname_list.append(elems[0].lstrip(' ')) # strip leading whitespace
					start_list.append(int(elems[3]))
					end_list.append(int(elems[4]))
					ID = elems[8].split(";")[0].split("=")[1]
					if ID in orthologous_groups.keys(): # match sequence ID with orthologous group
						ID_list.append(orthologous_groups[ID])
				
				except IndexError as ie:
					pass
		data=list(zip(seqname_list, start_list, end_list, ID_list))	

		### Select gene neighbor 
		# two genes are considered neighbor only if the distance between the end of one and the start of the following is less than 1000 bp
		GFF_paralogs = pd.DataFrame(data, columns=['seqname', 'start', 'end', 'ID_og'])
		all_IDs = GFF_paralogs['ID_og'].tolist()
		paralogs = list(set([x for x in all_IDs if all_IDs.count(x) > 1]))
		GFF_df = GFF_paralogs[~GFF_paralogs.ID_og.isin(paralogs)]
		GFF_df['next_seqname'] = GFF_df['seqname'].shift(-1)
		GFF_df['next_start'] = GFF_df['start'].shift(-1)
		GFF_df['next_ID'] = GFF_df['ID_og'].shift(-1)
		GFF_df['distance'] = GFF_df['next_start']-GFF_df['end']
		contig_control = GFF_df[GFF_df['seqname'] == GFF_df['next_seqname']] 
		distance_control = contig_control[contig_control['distance'] < 1000] # deleting neighbor genes distant more than 1000 bp
		neighbors = distance_control[['ID_og', 'next_ID']] # df with only neighbor genes			
		del GFF_df, contig_control, distance_control
		name_file= GFF_file.split(".fna")[0] + ".matrix"
		neighbors.to_csv(name_file, sep='\t', index=False, header=False)
	
				
#GFF_reader(in_file)
if __name__ == '__main__':
	for i in list_roary_files:
		GFF_reader(i)


list_matrix_files = glob.glob("*.matrix")
all_combinations = list(combinations(list_matrix_files, 2))
print(len(all_combinations))

#load couple
def number_of_HGT_events(couple):
	genomeA = pd.read_csv(couple[0], sep='\t', names=['gene1', 'gene2'])
	genomeB = pd.read_csv(couple[1], sep='\t', names=['gene1', 'gene2'])

#find uniq genes (to avoid traslocations)
	genomeA_genes = list(set(genomeA['gene1'].tolist() + genomeA['gene2'].tolist()))
	genomeB_genes = list(set(genomeB['gene1'].tolist() + genomeB['gene2'].tolist()))
	list_genes = list(genomeA_genes + genomeB_genes)
	uniq_genes = pd.DataFrame(list_genes, columns=['uniq_genes'])
	uniq_genes.drop_duplicates(inplace=True, keep=False)
	uniq_genes.reset_index(drop=True, inplace=True)

#find uniq edges for every genome
	uniq_edges_A_1 = uniq_genes.merge(genomeA, left_on='uniq_genes', right_on='gene1')
	uniq_edges_A_2 = uniq_genes.merge(genomeA, left_on='uniq_genes', right_on='gene2')
	uniq_edges_A = pd.concat([uniq_edges_A_1, uniq_edges_A_2])
	del uniq_edges_A_1, uniq_edges_A_2, uniq_edges_A['uniq_genes']

	uniq_edges_B_1 = uniq_genes.merge(genomeB, left_on='uniq_genes', right_on='gene1')
	uniq_edges_B_2 = uniq_genes.merge(genomeB, left_on='uniq_genes', right_on='gene2')
	uniq_edges_B = pd.concat([uniq_edges_B_1, uniq_edges_B_2])
	del uniq_edges_B_1, uniq_edges_B_2, uniq_edges_B['uniq_genes']

	uniq_edges = pd.concat([uniq_edges_A, uniq_edges_B])
	uniq_edges.reset_index(drop=True, inplace=True)

#count events of HGT
#	name_A = []
#	name_B = []
#	events = []
	ev = pd.crosstab(uniq_edges.gene1, uniq_edges.gene2)
	idx = ev.columns.union(ev.index)
	ev = ev.reindex(index = idx, columns=idx, fill_value=0)
	ev = ev.to_numpy()
	graph = csr_matrix(ev)
	n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)

	genomeA_name = couple[0].split(".matrix")[0] + ".fna"
	genomeB_name = couple[1].split(".matrix")[0] + ".fna"
#	name_A.append(genomeA_name)
#	name_B.append(genomeB_name)
#	events.append(n_components)
	return genomeA_name, genomeB_name, n_components
#	data=list(zip(name_A, name_B, events))
#	return pd.DataFrame(data)

def chunked_iterable(iterable, size):
	it = iter(iterable)
	while True:
		chunk = tuple(itertools.islice(it, size))
		if not chunk:
			break
		yield chunk

#parallel
def mp_handler(data):
	with mp.Pool(10) as pool:
		with open('GainLossEvents.csv', 'a') as f:
			writer=csv.writer(f, lineterminator = '\n', delimiter="\t")
			for result in pool.imap_unordered(number_of_HGT_events, data):			
				writer.writerow((result))


if __name__ == '__main__':
	for c in chunked_iterable(all_combinations, size=10):
		mp_handler(c)
