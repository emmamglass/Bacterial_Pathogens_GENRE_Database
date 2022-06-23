# PATHGENN_project

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

## essentialgenes.py usage

This is an UNFINISHED program to determine the essential genes across all pathobionts in PATHGENN. Will work on uploading a newer version.

## fastaextractor.py usage

This program takes the input as a .fasta file with ALL geome annotations in the singular .fasta file. This program then separates each individual annotated genome .fasta file into its own file with the title as it's NCBI taxid. 


## flux_sampling_rivanna.py usage

This program performs flux sampling for all models in a given folder. This version was used on the UVA supercomputer system Rivanna.

## fluxjobcreator.py usage

This program writes a .slurm job to be used on the supercomputer. It creates paralell processes for flux sampling of each model. This can be added to a .slurm script. 


## fluxsampling_pseudomonas.py usage

This program was adapted from fluxsampling.py above to be used specifically on pseudomonas aeruginosa isolates. This was used as part of a project in one of my classes.


## gapsplit.py usage

This program was NOT written by me. This program was created in Paul Jensen's group and is used for more efficient and better flux sampling. It is necessary to contain in this folder since many programs in this folder utilize this tool.

## gene_to_KO.py usage

gene_to_KO.py usage converts KEGG genes to KO values. 

## genesrxnmetabolitesfig.py usage

genesrxnsmetabolitesfig.py takes an input called phyluminfo.csv which contains taxonomic class information, and the number of genes, metabolites, and reactions in each model. This program then creates a sereies of boxplots that describe the number and range of genes reactions and metabolites in each PATHGENN model.

## genome_features.py usage

genome_features.py takes an input of PATRIC genome ids in a .txt file and annotates each genome using the RAST annotation server. The result is one .fasta file containing annotations for each genome, which then can be separated using fastaextractor.py

## genome_features_dna.py usage

genome_features_dna.py takes an input of PATRIC genome ids in a .txt file and annotates the rRNA sequences of each genome using the RAST annotation server. This can be used to annotate the genome for 16s rRNA sequences

## gettaxonomy.py usage

geettaxonomy.py takes an input of NCBI taxids and gets their specified taxonomy information. Then it saves that information to a .csv file

## k20.py usage 

k20.py is not a very useful program. pay it no mind

## metabolicphylogeny.py usage

metabolicphylogeny.py reads all .sbml in the folder and then can perform a variety of tasks. It will make a list of all present reactions across the models. The 'panreactome'. Then, it will make an empty dataframe based on the columns being the panreactome, and the rows being the taxids for each strain. Additionally, it will populate the dataframe with 1s or 0s based on if a given reaction is present in a given strain. Then, it can calculate the jaccard distance between two pairs of pathobionts, and create a pairwise similarity dataframe. Additionally, it can create a reaction presence heatmap, perform nmds on reaction presence, and finally, perform kmenas clustering.

## quickflux.py usage

This program is used to perform Flux sampling on all models in a folder. The sampling method incorporates gapsplit.py function by Paul Jensen's group. This program is meant to run on the UVA supercomputer rivanna.

## reactionannotationfig.py usage

reactionannotationfig.py creates a figure based on the reaction annotations for each model in PATHGENN. Additionally, it bins reaction annotations according to peripheral and universal metabolism categories. This creates figure 2e in the PATHGENN paper
























