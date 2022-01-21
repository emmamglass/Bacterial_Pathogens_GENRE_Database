# Bacterial_Pathogens_GENRE_Database
fluxsampling.py usage
fluxsampling.py takes all .sbml files in a folder, samples flux for all sbml files, performs NMDS,
#and clusters data based on gram status, pathogen name, aerobe status.

NOTE: .sbml files MUST take on the name format as follows
      genus_species.feature_protien.sbml

required inputs: --cluster_type, --sample_size
    --cluster_type
          what you would like to color clusters by. options include: 'gram', 'pathogen', 'aerobic'
    --sample_size
          how mnay flux samples would you like to take for a given reconstruction (default is 500)

NOTE: the --cluster_type pathogen input is TEMPORARILY non-functional if you do not have a specific 
      subset of bacterial species
