# Bacterial_Pathogens_GENRE_Database

## fluxsampling.py usage

fluxsampling.py takes all .sbml files in a folder, samples flux for all sbml files, performs NMDS, and clusters data based on gram status, pathogen name, aerobe status.

NOTE: .sbml files MUST take on the name format as follows: genus_species.feature_protien.sbml


### required inputs:

#### --cluster_type
What you would like to color clusters by. options include: 'gram', 'pathogen', 'aerobic'

NOTE: the --cluster_type pathogen input is TEMPORARILY non-functional if you do not have a specific subset of bacterial species

#### --sample_size
How many flux samples would you like to take for a given reconstruction (default is 500)

## fastaextractor.py usage

fastaextractor.py takes a file called output.fasta that contains protein fasta sequences of all strains. Then, it separates the output.fasta into individual fasta files based on genomeid in the form of genomeid.fa. 

## jobcreator.py usage

jobcreator.py creates a .txt file with all srun commands needed to create all .sbml GENRE reconstructions of selected strains. this .txt file can be copied and pasted into a .slurm script to create reconstructions in parallel on rivanna. 

## 16S_extractor.py usage

16S_extractor.py uses annotated rna FASTA sequences as an input (can be multiple rna FASTA sequences from multiple species), and extracts the 16S sequence from the annotated rna FASTA input file. If a given record does not contain a 16S sequence, no error is raised, just nothing is added to the ouptut file.

## Genesvsrxns.py usage

Genesvsrxns.py reads in a taxonomy info file that has a column with all species names, another column with the number of genes in the GENERE corresponding to the species, and a third column wih the number of reactions in the GENRE corresponding to the species. This program then creates a simple plot of genes vs reactions for each model.

## GetPairwiseEssentialGenes.py usage

GetPairwiseEssentialGenes.py performs a gene essentiality screen on each of 914 GENREs in the PATHGENN database. The gene essentiality screen is performed using the cobra.flux_analysis.variability.find_essential_genes(model) function, and adds a list of essential genes to a column corresponding to a given taxid for a given GENRE. So, there will be 914 genres, each with an essential gene list. 

## NMDS.py usage

NMDS.py uses isolation source data, phylum data, class data, and taxid data, and reaction presence data or flux data and performs NMDS dimensionsionality reduction on either reaction presence or flux data, coloring on either isolation source, phylum data, or class data, which can be altered in the code to the user's liking. 

## Reaction_annotations.py usage

Reaction_annotations.py takes a list of species specific genes and converts them into KEGG orthologs by interfacing with KEGG. If a given essential gene is not found in the KEGG database, NA will be inserted in its place.



























