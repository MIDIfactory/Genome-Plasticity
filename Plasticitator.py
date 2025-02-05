## LIBRARIES
import argparse,sys,time,datetime,os,json,pickle,io,re,shutil,glob
import pandas as pd
from multiprocessing import Pool
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from itertools import combinations
from Bio import SeqIO


## FUNCTIONS
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def parse_args():

    parser=MyParser(usage='%(prog)s -gf GROUP_FILE -t NUM_THREADS',
        prog='Plasticitator.py',
        argument_default=argparse.SUPPRESS,
        epilog="FOGS made easy...ish\nREMEMBER! You can use the same genome in different categories")

    parser.add_argument("-gf", help="Group file formatted as follows: genome_name{TAB}genome_path{TAB}group_names (this is also the header).\nGenomes can be in multiple groups, use commas in the last column",dest="groups",required=False)
    parser.add_argument("-snpr", help="Path to reference genome to be used to call coreSNP",dest="snp_ref",required=False)
    parser.add_argument('-t', type=int,help='Number of threads',dest="threads",required=False)

    args = parser.parse_args()
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()
    elif (args.groups is None or args.threads is None or args.snp_ref is None):
        parser.error("\n\nERROR!! If not in TEST mode, the following arguments are required: -gf, -snpr, -t\n\n")
    return args

def rename_prots(i):
    i=i.split("-")
    assert i[0] in samplelist
    i=i[0]+"||"+"-".join(i[1:])
    return i

def GFF_reader(gname,gvalues):
    GFF_file=gvalues[0]
    orthologous_groups=gvalues[1]
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

            elems = chomped_line.split('\t')
            try:
                seqname_list.append(elems[0].lstrip(' ')) # strip leading whitespace
                start_list.append(int(elems[3]))
                end_list.append(int(elems[4]))

                ID = gname+"||"+elems[0]+"-"+elems[8].split(";")[0].split("=")[1]
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
        neighbors=[[row["ID_og"],row["next_ID"]] for index,row in neighbors.iterrows()]
        return gname, neighbors

def number_of_HGT_events(GA,GB):
    genomeA = pd.DataFrame(neighbor_matrices[GA],columns=['gene1', 'gene2'])
    genomeB = pd.DataFrame(neighbor_matrices[GB],columns=['gene1', 'gene2'])

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
    ev = pd.crosstab(uniq_edges.gene1, uniq_edges.gene2)
    idx = ev.columns.union(ev.index)
    ev = ev.reindex(index = idx, columns=idx, fill_value=0)
    ev = ev.to_numpy()
    graph = csr_matrix(ev)
    n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)

    return [GA, GB, n_components]


def core_snps(list_genomes, window,reference_fasta):
    diz_pos={}
    diz_pos_ref={}
    for gg in list_genomes:
        p=open("snps/%s.snp" %gg,'r')
        for x in p:
            if not re.search(r'\.',x):
                pos=int(x.split('\t')[0])
                base_ref=x.split('\t')[1]
                base_org=x.split('\t')[2]
                letters = set('atcgATCG')
                try:
                    diz_pos[pos]
                except:
                    diz_pos.update({int(pos):'0'})
                    diz_pos_ref.update({int(pos):base_ref})
            if not ((letters & set(base_org)) and (letters & set(base_ref))):
                diz_pos.update({int(pos):'1'})
        p.close()
    keys_sort=sorted(diz_pos.keys())
    core=[]
    if ((keys_sort[1]-keys_sort[0]) > int(window)) and (diz_pos[keys_sort[0]] == "0"):
        core.append(keys_sort[0])
    for i in range(1,len(keys_sort)-2):
        if (keys_sort[i]-keys_sort[i-1]) > int(window) and (keys_sort[i+1]-keys_sort[i]) > int(window) and (diz_pos[keys_sort[i]] == "0"):
            core.append(keys_sort[i])
    if (keys_sort[-1]-keys_sort[-2]) > int(window) and (diz_pos[keys_sort[-1]] == "0"):
        core.append(keys_sort[-1])

    count=0
    for seq_record in SeqIO.parse(reference_fasta,"fasta"):
        list_pos_ref=(list(seq_record.seq.lower()))
        reference=seq_record.id
        count=count+1
        if count>1:
            print('Error')
            break
    outname_snp_fasta="coreSNPs.fna"
    outfasta=open(outname_snp_fasta,"w")
    s=''
    for i in core:
        list_pos_ref[int(i)-1].strip()
        s+=list_pos_ref[int(i)-1]
    outfasta.write(">REF\n%s\n" %(s.strip().upper()))
    for gg in list_genomes:
        basi={}
        p2=open("snps/%s.snp" %gg,'r')
        for x2 in p2:
            if not re.search(r'\.',x2):
                pos2=int(x2.split('\t')[0])
                base_org2=x2.split('\t')[2]
                basi.update({int(pos2):base_org2.strip()})
        s=''
        for i in core:
            try:
                basi[int(i)]
                s+=basi[int(i)]
            except:
                list_pos_ref[int(i)-1].strip(),
                s+=list_pos_ref[int(i)-1]
        outfasta.write(">%s\n%s\n" %(gg,s.strip().upper()))
        p2.close()
    outfasta.close()



## STATICS
RAM_tmp="/dev/shm" ## temporary folder on ramdisk [default=RAM tmp in ubuntu]
###KEEP IT EFFING TIDY!!! (i.e empty after use)


####BUGFIX for conda env within env
condabase="/mnt/bunnydisk/greta/miniconda3/envs/panta/bin/panta"


##MAIN
args = parse_args()
time_now=datetime.datetime.now()
datestamp=time_now.strftime("%Y-%m-%d_%H-%M-%S")
Results_folder_name= str(str(args.groups) + "_RUN_Plasticity")
print("Resuls will be saved to %s" %Results_folder_name)
os.mkdir(Results_folder_name)
os.mkdir("%s/genomes" %Results_folder_name)
THIS_DIR=os.path.abspath( os.path.dirname( sys.argv[0]))


df=pd.read_csv(str(args.groups),sep="\t", dtype=str)
df = df.set_index("genome_name")

for index, row in df.iterrows():
    os.system("cp %s %s/genomes/%s.fna" %(row["genome_path"],Results_folder_name,index))

samplelist=df.index.tolist()

REF=os.path.abspath(args.snp_ref)
os.mkdir("%s/tmp_plasticity" %(RAM_tmp))
os.system("cp -r %s/genomes %s/tmp_plasticity" %(Results_folder_name,RAM_tmp))
os.chdir("%s/tmp_plasticity" %RAM_tmp)
os.mkdir("gff")
os.system("echo '##FASTA' > mid")
os.system("ls genomes/* | parallel -j %i 'prodigal -i {} -f gff -o {}.gff 2>/dev/null ; cat {}.gff mid {} > {}.GFF'" %args.threads)
os.system("rm genomes/*.gff")

### CHANGED FROM STE

# Define the source and destination directories
source_directory = 'genomes'
destination_directory = 'gff'

# Create the destination directory if it doesn't exist
os.makedirs(destination_directory, exist_ok=True)

# Use glob to find all files matching the pattern
source_files = glob.glob(os.path.join(source_directory, '*.fna.GFF'))

# Loop through each file and rename/move it
for file_path in source_files:
    # Get the base name of the file (without the directory)
    base_name = os.path.basename(file_path)
    
    # Change the extension and create the new file name
    new_file_name = base_name.replace('.fna.GFF', '.gff')
    
    # Create the full destination path
    destination_path = os.path.join(destination_directory, new_file_name)
    
    # Move (rename) the file to the new location
    shutil.move(file_path, destination_path)

print("Files have been renamed and moved.")

###
os.system("/mnt/bunnydisk/greta/miniconda3/envs/panta/bin/panta main -g gff/*.gff -o panta -t %i " %(args.threads))
os.system("cp -r panta %s/%s/panta; cp -r gff %s/%s/gff" %(THIS_DIR,Results_folder_name,THIS_DIR,Results_folder_name))
clusters = json.load(open('panta/clusters.json'))
cluster_assign={}
c=0
for i in clusters:
    c+=1
    others=clusters[i]
    i=rename_prots(i)
    cluster_assign[i]="Group%i" %c
    for o in others:
        cluster_assign[rename_prots(o)]="Group%i" %c

with open('%s/%s/cluster_assign.pkl' %(THIS_DIR,Results_folder_name), 'wb') as f:
    pickle.dump(cluster_assign, f)


mpinput={g:["gff/%s.gff" %(g),cluster_assign] for g in samplelist}

with Pool(args.threads) as p: #####DEBUG WITH 1
    results=p.starmap(GFF_reader, mpinput.items())
neighbor_matrices = {key: value for key, value in results}

with open('%s/%s/neighbor_matrices.pkl' %(THIS_DIR,Results_folder_name), 'wb') as f:
    pickle.dump(neighbor_matrices, f)

categories={}
for index,row in df.iterrows():
    cats=row["group_names"].strip().split(",")
    for cat in cats:
        if cat not in categories.keys():
            categories[cat]=[]
        categories[cat].append(index)

os.mkdir("snps")
os.system('ls genomes/*.fna | parallel -j %i "nucmer --maxgap=500 -p {} %s {}; delta-filter -1 {}.delta > {}_filtered.delta; show-snps -H -C -I -r -T {}_filtered.delta | cut -f1,2,3 >{}.snp"' %(args.threads,REF))
#os.system("mmv 'genomes/*.fna.snp' 'snps/#1.snp'")
input_dir = 'genomes'
output_dir = 'snps'

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Find all files matching the *.fna.snp pattern in the input directory
input_files = glob.glob(os.path.join(input_dir, '*.fna.snp'))

# Loop through each file and rename/move it
for snp_file in input_files:
    # Get the base name of the file (e.g., file1.fna.snp)
    base_name = os.path.basename(snp_file)
    
    # Create the new file name by replacing the extension .fna.snp with .snp
    new_file_name = base_name.replace('.fna.snp', '.snp')
    
    # Construct the full path for the new file in the output directory
    destination_path = os.path.join(output_dir, new_file_name)
    
    # Move and rename the file
    shutil.move(snp_file, destination_path)



os.system("cp -r snps %s/%s/snps" %(THIS_DIR,Results_folder_name))
core_snps(samplelist,10,REF)
os.system("cp -r coreSNPs.fna %s/%s/coreSNPs.fna" %(THIS_DIR,Results_folder_name))
SNP_distances=os.popen("snp-dists -m coreSNPs.fna").read().strip("\n").split("\n")
SNP_distances={"%s||%s" %(r.split("\t")[0],r.split("\t")[1]) : int(r.split("\t")[2]) for r in SNP_distances}

with open('%s/%s/SNP_distances.pkl' %(THIS_DIR,Results_folder_name), 'wb') as f:
    pickle.dump(SNP_distances, f)


GS_distances={}
for cat in categories:
    all_combinations = list(combinations(categories[cat], 2))
    all_combinations = [[g[0],g[1]] for g in all_combinations]
    with Pool(args.threads) as p: #####DEBUG WITH 1 ||||| args.threads
        results=p.starmap(number_of_HGT_events, all_combinations)
    GS_distances[cat]={"%s||%s" %(r[0],r[1]):int(r[2]) for r in results}

with open('%s/%s/GS_distances.pkl' %(THIS_DIR,Results_folder_name), 'wb') as f:
    pickle.dump(GS_distances, f)

def get_wGS(row):
    snpvalue=float(row['snps'])
    if snpvalue==0:
        snpvalue=1.0
    wGS=float(row['GS'])/snpvalue
    return wGS


FinalDF=pd.DataFrame()
FOGSfile=open("%s/%s/FOGS.txt" %(THIS_DIR,Results_folder_name),"w")
for cat in categories:
    diz={k:{"GS":GS_distances[cat][k],"snps":SNP_distances[k],"group":cat} for k in GS_distances[cat].keys()}
    diz=pd.DataFrame.from_dict(diz,orient="index")
    diz["wGS"]=diz.apply(get_wGS, axis=1)
    diz["pair"]=diz.index
    wGS=diz["wGS"].tolist()
    avg=sum(wGS)/len(wGS)
    FOGSfile.write("%s\t%f\n" %(cat,avg))
    diz=diz.reset_index(drop=True)
    if FinalDF.empty:
        FinalDF=diz
    else:
        FinalDF=pd.concat([FinalDF, diz], ignore_index=True)
FOGSfile.close()

FinalDF.to_csv("%s/%s/FOGStable.tsv" %(THIS_DIR,Results_folder_name), sep="\t")


os.chdir(THIS_DIR)
os.system("rm -rf %s/tmp_plasticity" %RAM_tmp)
