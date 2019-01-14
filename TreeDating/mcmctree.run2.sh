############################################################################################################
# 2nd MCMCtree RUN
############################################################################################################

############################################################################################################
#0. setup environment
############################################################################################################

base=~/data/inqGen18/phylogeny/treeDating/MCMCtree
bd=~/data/inqGen18/phylogeny/
software=/usr/local/home/lschrader/software
scripts=~/data/inqGen18/phylogeny/scripts/
SEQFILE=$base/data/4dAll.fa.phy
TOPOLOGYFILE=$base/data/4dAll.topo.phy
MCMCtreeInput=$base/data/4dAll.mcmcTreeInput.phy

############################################################################################################
# Run mcmctree
#     sampFreq = 5000
#     substitution_rate = 0.394531 estimated from @17.5 calibration
#     sigma2_gamma=5.5
# use default priors (gamma-dirichlet prior for locus rates (RGENE) and conditional iid for sigma2)
# use default alphas (1)
############################################################################################################

SAMPFREQ=5000
ALPHA=$(cat $base/baseml.17.5calibration/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml.17.5calibration/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.394531/0.394531)^2
#  RGENE2 = 0.394531/0.394531^2
RGENE1=$(echo "($SR/$SR)^2" | bc -l)
RGENE2=$(echo "$SR/($SR)^2" | bc -l)
SIGMA1=1
SIGMA2=5.5 #tested different values, 0.9, 1, 4.5, 10 and 5.5 seems well fitting
#############################
#1.1 First run of MCMCtree
#############################

#The branch lengths are estimated by maximum likelihood together with the gradient and Hessian
#(i.e. vector of first derivatives and matrix of second derivatives) of the likelihood function at the maximum likelihood estimates.
#The gradient and Hessian contain information about the curvature of the likelihood surface (dos Reis and Yang, 2013, P10)．

mkdir $base/mcmctree.sf500.sr17.5cal.sigma2.5.5.run2/
cd $base/mcmctree.sf500.sr17.5cal.sigma2.5.5.run2/

cat $base/ctl/mcmctree.template.ctl | \
sed "s|<SIGMA1>|$SIGMA1|g" | \
sed "s|<SIGMA2>|$SIGMA2|g" | \
sed "s|<SEQFILE>|$SEQFILE|g" | \
sed "s|<SAMPFREQ>|$SAMPFREQ|g" | \
sed "s|<MCMCtreeInput>|$MCMCtreeInput|g" | \
sed "s|<ALPHA>|$ALPHA|g" | \
sed "s|<RGENE1>|$RGENE1|g" | \
sed "s|<RGENE2>|$RGENE2|g" > ./mcmctree1.ctl

mcmctree ./mcmctree1.ctl

#############################
#1.2 Second run of MCMCtree
#############################
cp out.BV in.BV
cat mcmctree1.ctl | \
sed "s/usedata = 3/usedata = 2/g" | \
sed "s/outfile = out/outfile = out_usedata2/g" > mcmctree2.ctl

# change usedata to 2 in mcmctree2ndrun.ctl
# change outfile to out_usedata2
mcmctree mcmctree2.ctl
mv mcmc.txt p.mcmc.run2.txt
mv FigTree.tre FigTree.run2.tre
############################################################################################################
