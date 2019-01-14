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
# 1. Run mcmctree with sampFreq500
############################################################################################################

SAMPFREQ=500
ALPHA=$(cat $base/baseml1/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml1/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.460286/0.460286)^2
#  RGENE2 = 0.460286/0.460286^2
RGENE1=1
RGENE2=2.173

#############################
#1.1 First run of MCMCtree
#############################

#The branch lengths are estimated by maximum likelihood together with the gradient and Hessian
#(i.e. vector of first derivatives and matrix of second derivatives) of the likelihood function at the maximum likelihood estimates.
#The gradient and Hessian contain information about the curvature of the likelihood surface (dos Reis and Yang, 2013, P10)．

mkdir $base/mcmctree.run.sampfreq500/
cd $base/mcmctree.run.sampfreq500/

cat $base/ctl/mcmctree.template.ctl | \
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
cd $base/mcmctree.run.sampfreq500/
cp out.BV in.BV
cat mcmctree1.ctl | \
sed "s/usedata = 3/usedata = 2/g" | \
sed "s/outfile = out/outfile = out_usedata2/g" > mcmctree2.ctl

# change usedata to 2 in mcmctree2ndrun.ctl
# change outfile to out_usedata2

mcmctree mcmctree2.ctl
############################################################################################################


############################################################################################################
# 2. Run mcmctree on posterior only (1st run)
############################################################################################################
SAMPFREQ=500
ALPHA=$(cat $base/baseml1/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml1/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.460286/0.460286)^2
#  RGENE2 = 0.460286/0.460286^2
RGENE1=$(echo "($SR/$SR)^2" | bc -l)
RGENE2=$(echo "$SR/($SR)^2" | bc -l)

#############################
#1.1 First run of MCMCtree
#############################

mkdir $base/mcmctree.p1/
cd $base/mcmctree.p1/

cat $base/ctl/mcmctree.template.ctl | \
sed "s|usedata = 3|usedata = 0|g" | \
sed "s|<SEQFILE>|$SEQFILE|g" | \
sed "s|<SAMPFREQ>|$SAMPFREQ|g" | \
sed "s|<MCMCtreeInput>|$MCMCtreeInput|g" | \
sed "s|<ALPHA>|$ALPHA|g" | \
sed "s|<RGENE1>|$RGENE1|g" | \
sed "s|<RGENE2>|$RGENE2|g" > ./mcmctree.p1.ctl

mcmctree ./mcmctree.p1.ctl

############################################################################################################

############################################################################################################
# 3. Run mcmctree on posterior only (2nd run)
############################################################################################################
SAMPFREQ=500
ALPHA=$(cat $base/baseml1/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml1/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.460286/0.460286)^2
#  RGENE2 = 0.460286/0.460286^2
RGENE1=$(echo "($SR/$SR)^2" | bc -l)
RGENE2=$(echo "$SR/($SR)^2" | bc -l)

#############################
#1.1 First run of MCMCtree
#############################

mkdir $base/mcmctree.p2/
cd $base/mcmctree.p2/

cat $base/ctl/mcmctree.template.ctl | \
sed "s|usedata = 3|usedata = 0|g" | \
sed "s|<SEQFILE>|$SEQFILE|g" | \
sed "s|<SAMPFREQ>|$SAMPFREQ|g" | \
sed "s|<MCMCtreeInput>|$MCMCtreeInput|g" | \
sed "s|<ALPHA>|$ALPHA|g" | \
sed "s|<RGENE1>|$RGENE1|g" | \
sed "s|<RGENE2>|$RGENE2|g" > ./mcmctree.p2.ctl

mcmctree ./mcmctree.p2.ctl

############################################################################################################

############################################################################################################
# 4. Run mcmctree with sampFreq50 and substitution rate 0.394531 estimated from @17.5 calibration
############################################################################################################

SAMPFREQ=50
ALPHA=$(cat $base/baseml2/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml2/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.394531/0.394531)^2
#  RGENE2 = 0.394531/0.394531^2
RGENE1=$(echo "($SR/$SR)^2" | bc -l)
RGENE2=$(echo "$SR/($SR)^2" | bc -l)

#############################
#1.1 First run of MCMCtree
#############################

#The branch lengths are estimated by maximum likelihood together with the gradient and Hessian
#(i.e. vector of first derivatives and matrix of second derivatives) of the likelihood function at the maximum likelihood estimates.
#The gradient and Hessian contain information about the curvature of the likelihood surface (dos Reis and Yang, 2013, P10)．

mkdir $base/mcmctree.run.sampfreq50.17.5cal/
cd $base/mcmctree.run.sampfreq50.17.5cal/

cat $base/ctl/mcmctree.template.ctl | \
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
cd $base/mcmctree.run.sampfreq50.17.5cal/
cp out.BV in.BV
cat mcmctree1.ctl | \
sed "s/usedata = 3/usedata = 2/g" | \
sed "s/outfile = out/outfile = out_usedata2/g" > mcmctree2.ctl

# change usedata to 2 in mcmctree2ndrun.ctl
# change outfile to out_usedata2
mcmctree mcmctree2.ctl
############################################################################################################


############################################################################################################
# 5. Run mcmctree with sampFreq500 and substitution rate 0.394531 estimated from @17.5 calibration and sigma=2.5
############################################################################################################

SAMPFREQ=500
ALPHA=$(cat $base/baseml.17.5calibration/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml.17.5calibration/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.394531/0.394531)^2
#  RGENE2 = 0.394531/0.394531^2
RGENE1=$(echo "($SR/$SR)^2" | bc -l)
RGENE2=$(echo "$SR/($SR)^2" | bc -l)
SIGMA1=1
SIGMA2=2.5 #formerly 4.5
#############################
#1.1 First run of MCMCtree
#############################

#The branch lengths are estimated by maximum likelihood together with the gradient and Hessian
#(i.e. vector of first derivatives and matrix of second derivatives) of the likelihood function at the maximum likelihood estimates.
#The gradient and Hessian contain information about the curvature of the likelihood surface (dos Reis and Yang, 2013, P10)．

mkdir $base/mcmctree.run.sampfreq500.17.5cal.sigma.2.5/
cd $base/mcmctree.run.sampfreq500.17.5cal.sigma.2.5/

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
############################################################################################################

############################################################################################################
# 6. Run mcmctree with sampFreq500 and substitution rate 0.394531 estimated from @17.5 calibration and sigma=2.5,
# but use default priors (gamma-dirichlet prior for locus rates (RGENE) and conditional iid for sigma2 (I think)) and alphas (1)
############################################################################################################

SAMPFREQ=500
ALPHA=$(cat $base/baseml.17.5calibration/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml.17.5calibration/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.394531/0.394531)^2
#  RGENE2 = 0.394531/0.394531^2
RGENE1=$(echo "($SR/$SR)^2" | bc -l)
RGENE2=$(echo "$SR/($SR)^2" | bc -l)
SIGMA1=1
SIGMA2=2.5 #formerly 4.5
#############################
#1.1 First run of MCMCtree
#############################

#The branch lengths are estimated by maximum likelihood together with the gradient and Hessian
#(i.e. vector of first derivatives and matrix of second derivatives) of the likelihood function at the maximum likelihood estimates.
#The gradient and Hessian contain information about the curvature of the likelihood surface (dos Reis and Yang, 2013, P10)．

mkdir $base/mcmctree.run.sampfreq500.17.5cal.sigma.2.5.defaultPriors/
cd $base/mcmctree.run.sampfreq500.17.5cal.sigma.2.5.defaultPriors/

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
############################################################################################################


############################################################################################################
# 6. Run mcmctree with sampFreq50 and substitution rate 0.394531 estimated from @17.5 calibration and sigma=0.9,
# but use default priors (gamma-dirichlet prior for locus rates (RGENE) and conditional iid for sigma2 (I think)) and alphas (1)
############################################################################################################

SAMPFREQ=50
ALPHA=$(cat $base/baseml.17.5calibration/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml.17.5calibration/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.394531/0.394531)^2
#  RGENE2 = 0.394531/0.394531^2
RGENE1=$(echo "($SR/$SR)^2" | bc -l)
RGENE2=$(echo "$SR/($SR)^2" | bc -l)
SIGMA1=1
SIGMA2=0.9 #formerly 4.5
#############################
#1.1 First run of MCMCtree
#############################

#The branch lengths are estimated by maximum likelihood together with the gradient and Hessian
#(i.e. vector of first derivatives and matrix of second derivatives) of the likelihood function at the maximum likelihood estimates.
#The gradient and Hessian contain information about the curvature of the likelihood surface (dos Reis and Yang, 2013, P10)．

mkdir $base/mcmctree.run.sampfreq500.17.5cal.sigma.0.9.defaultPriors/
cd $base/mcmctree.run.sampfreq500.17.5cal.sigma.0.9.defaultPriors/

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
mv mcmc.txt mcmc.run.sampfreq500.17.5cal.sigma.0.9.defaultPriors.txt
############################################################################################################


############################################################################################################
# 8. Run mcmctree with sampFreq50 and substitution rate 0.394531 estimated from @17.5 calibration and sigma=10,
# but use default priors (gamma-dirichlet prior for locus rates (RGENE) and conditional iid for sigma2 (I think)) and alphas (1)
############################################################################################################

SAMPFREQ=50
ALPHA=$(cat $base/baseml.17.5calibration/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')
SR=$(sed 1d $base/baseml.17.5calibration/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
# calculate RGENE1 and RGENE2 as described in 1mcmctree.readme.txt
#  RGENE1 = (0.394531/0.394531)^2
#  RGENE2 = 0.394531/0.394531^2
RGENE1=$(echo "($SR/$SR)^2" | bc -l)
RGENE2=$(echo "$SR/($SR)^2" | bc -l)
SIGMA1=1
SIGMA2=10 #formerly 4.5
#############################
#1.1 First run of MCMCtree
#############################

#The branch lengths are estimated by maximum likelihood together with the gradient and Hessian
#(i.e. vector of first derivatives and matrix of second derivatives) of the likelihood function at the maximum likelihood estimates.
#The gradient and Hessian contain information about the curvature of the likelihood surface (dos Reis and Yang, 2013, P10)．

mkdir $base/mcmctree.run.sampfreq500.17.5cal.sigma.10.defaultPriors/
cd $base/mcmctree.run.sampfreq500.17.5cal.sigma.10.defaultPriors/

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
mv mcmc.txt mcmc.run.sampfreq500.17.5cal.sigma.10.defaultPriors.txt
############################################################################################################
