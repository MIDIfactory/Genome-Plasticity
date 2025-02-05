# Genome-Plasticity

Scripts used for computing Flux of Genes Segments (FOGS), which enable to quantify the rate of gene acquisition or loss events within a collection of bacterial genomes.


![Schermata del 2025-02-05 09-57-01](https://github.com/user-attachments/assets/64b14d15-0b41-4f1a-9fb3-f262ba98e926)
where 𝛸 is the number of gene gain/loss events, d is the SNP distance between genome A and genome B, N is the number of genomes considered). The higher the value of FOGS, the higher is the genome plasticity. 

## Installation

### Requirments
Plasticitator requires **conda**, [Python 3](https://www.python.org/) and the following dependencies:
* argparse
* pandas
* multiprocessing
* scipy
* BioPython
* glob, shutil, re, json, pickle
* prodigal (for gene annotation, ensure it is installed and available in PATH)

### Setup
```
git clone https://github.com/yourusername/Plasticitator.git
cd Plasticitator
```

## Usage
### Command-line Arguments

```
python Plasticitator.py -gf GROUP_FILE -snpr REFERENCE_GENOME -t NUM_THREADS
```

Required Arguments:

```-gf```: Path to the group file, formatted as: genome_name \t genome_path \t group_names.

```-snpr```: Path to the reference genome for SNP calling.

```-t```: Number of threads for parallel execution.


### Example

```
python Plasticitator.py -gf genomes.tsv -snpr ref_genome.fna -t 8
```
Input Files:

Group File (-gf): A tab-separated file listing genomes with their paths and associated groups.

Reference Genome (-snpr): A FASTA file to be used for core SNP identification.


### Output

Results are saved in a directory named groupname_RUN_Plasticity, containing:

```cluster_assign.pkl```: A pickle file storing clustering assignments.

```coreSNPs.fna```: A FASTA file containing core SNP sequences identified across genomes, obtained using [P-DOR](https://github.com/gtonkinhill/fastbaps).

```FOGStable.tsv```: A tab-separated file containing the number of gene gain/loss events, core SNP distance, cluster, and gene gain/loss events weighted by SNP distance for each genome pair

```FOGS.txt```: A text file containing the FOGS value for each identified cluster.

```GS_distances.pkl```: A pickle file containing FOGS calculations for each genome pair.

```neighbor_matrices.pkl```: A pickle file with gene neighborhood matrices for adjacency analysis.

```SNP_distances.pkl```: A pickle file containing SNP-based genetic distance calculations.

```snps/```: A directory storing raw SNP information for further downstream analysis.

```genomes/```: Copies of input genomes.

```panta/```: A directory contaning [panta](https://github.com/amromics/panta) results.


### Warning

Plasticitator creates and stores temporary files during execution in the ```/dev/shm`` directory. If the execution is interrupted, temporary files may remain in ```/dev/shm```, consuming memory. It is advisable to clean up any leftover files manually to free up system resources using the following command: 
```
rm -rf /dev/shm/tmp_plasticity
```

