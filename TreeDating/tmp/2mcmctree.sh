#2. Second run of MCMCtree
#########
base=~/data/inqGen18/phylogeny/treeDating/MCMCtree
bd=~/data/inqGen18/phylogeny/
software=/usr/local/home/lschrader/software
scripts=~/data/inqGen18/phylogeny/scripts/

SEQFILE=$base/data/4dAll.fa.phy
TOPOLOGYFILE=$base/data/4dAll.topo.phy
ALPHA=1.47227
RGENE1=1
RGENE2=2.173
MCMCtreeInput=$base/data/4dAll.mcmcTreeInput.phy

cp out.BV in.BV
cp mcmctree.ctl mcmctree2ndrun.ctl

cat $base/ctl/mcmctree1.ctl | \
sed "s/usedata = 3/usedata = 2/g" | \
sed "s/outfile = out/outfile = out_usedata2/g" > $base/ctl/mcmctree2.ctl

# change usedata to 2 in mcmctree2ndrun.ctl
# change outfile to out_usedata2

mcmctree $base/ctl/mcmctree2.ctl


#Adjust the fine tuning in case it is not between 0.2-0.4
# Now let’s look at the first line of the MCMC proper:# example:     -4% 0.16 0.64 0.33 0.00 0.69  0.178 0.149 0.091 0.067 0.032 - 1.438 -34999.3
# my results:  -4% 0.01 0.83 0.08 0.99 0.00  0.237 0.127 0.054 0.047 0.013 - 0.392 -3.8
# my results2: -4% 0.28 0.40 0.46 0.52 0.00  0.603 0.354 0.258 0.209 0.013 - 0.771 -5.2 #after manual adjustment of finetune
# The negative percentage (−4%) indicates that we are in the burn-in stage of the MCMC. The next 5 numbers are the acceptance proportions. They are printed in the order times, rates, mixing,
# substitution model parameters, and rate param- eters. For example, 16% of all proposed times were accepted during this stage of the MCMC (i.e. 84% were rejected), while 64% of the rates propos ed
# were accepted. A good MCMC analysis should have acceptance proportions close to 30% (20-40% being a good range and 15-70% being acceptable). You can see that the program goes through various
# rounds of finetune improvement until the acceptance proportions get very close to 30%:

# A few seconds after you start the program, the screen output will look like the above. Scroll back to check that the tree, fossil calibrations and the sequence alignments are read correctly. The output here is generated from a
# run under the JC model and global clock (clock = 1). The percentage % indicates the progress of the run, with negative values for the burn-in. Then the five ratios (e.g., 0.33 0.01 0.25 0.00 0.00 on the first line) are the
# acceptance proportions for the corresponding proposals. The optimal acceptance proportions are around 0.3, and you should try to make them fall in the interval (0.2, 0.4) or at least (0.15, 0.7). If the acceptance proportion is
# too small (say, <0.10), you decrease the corresponding finetune parameter. If the acceptance proportion is too large (say, >0.80), you increase the corresponding finetune parameter.

#Output
# ~/CSE/inquilineGenomics/Phylogeny/divergenceTime/FigTree.tre
