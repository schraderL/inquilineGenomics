# Calculate divergence time with MCMC tree
######################
#http://www.fish-evol.com/mcmctreeExampleVert6/text1Eng.html

#Calibration from:
# Nygaard et al 2011:														Atta-Acro split 10 MYA
# Mikheyev et al 2010:														Atta-Acro split 10-16 or 8-12 MYA
# Hoelldobler & Wilson The Leafcutter Ants: Civilization by Instinct:		Atta-Acro split 8-12 MYA
# Schultz & Brady 2008 used fossil data										Origin of leaf-cutting 8-12 MYA

# !!! But Nygaard et al 2016													Atta-Acro split ~11-20 MYA
# From the supplement of the Nygaard 2016, the split between Atta and Acro occured 16.2 mya (mean) 5%-95% (12.6 mya - 19.7 mya)
#Input:
#alignment
#/Users/lukas/CSE/inquilineGenomics/Phylogeny/input/cat.sequential.phy
#tree with specific point calibration (12.6 - 19.7)
#/Users/lukas/CSE/inquilineGenomics/Phylogeny/tree.baseml
#tree with range calibration (16.2 mya)
#/Users/lukas/CSE/inquilineGenomics/Phylogeny/MCMCtree_Input.tre

#0. baseml
#########

# run baseml to get estimate for substitution rate per time unit
# For this, I need a specific calibration date. I used the mean for the Atta-Acro split as calculated from Haofu's alignment.
# Split between Acro and Atta at mean 16.2 mya
cd /Users/lukas/sciebo/inquilineGenomics18/phylogeny/TreeDating/MCMCtree/
baseml ./ctl/baseml.ctl

#Check file mlb for substitution rate and alpha
#########

#Substitution rate is per time unit
#    0.313239 +- 0.001812

#Alpha
# alpha (gamma, K=5) =  1.54319
