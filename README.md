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


