#!/usr/bin/bash


#Performing the initial clustering
../clust complexes_allgrids_c2c_nspec p2_noh.pdb 5 'a' 500

#Doing the treecutting 

../cut cluster_average_linkage_full_scoring.out complexes_allgrids_c2c_nspec p2_noh.pdb 5 'a' > stats.out

#Getting increments, avd. ClSize, min. ClSize and max. ClSize columns
../aux/complete_scoring.awk cluster_average_linkage.out > cluster_average_linkage_full_scoring.out

#Plotting the results
../aux/plot_cycles.gnu cluster_average_linkage_full_scoring.out plot.pdf "bnbs docking" 500
