#################################################################################
# Calculate divergence time with MCMC tree
#################################################################################
# Notes:
##############
## good tutorials
#http://www.fish-evol.com/mcmctreeExampleVert6/text1Eng.html
# https://github.com/mariodosreis/divtime

## Calibration from:
  # Nygaard et al 2011:Atta-Acro split 10 MYA
  # Mikheyev et al 2010:Atta-Acro split 10-16 or 8-12 MYA
  # Hoelldobler & Wilson The Leafcutter Ants: Civilization by Instinct:Atta-Acro split 8-12 MYA
  # Schultz & Brady 2008 used fossil data:Origin of leaf-cutting 8-12 MYA
  # Nygaard et al 2016:Atta-Acro split ~11-20 MYA
  # From the supplement of the Nygaard 2016, the split between Atta and Acro occured 16.2 mya (mean) 5%-95% (12.6 mya - 19.7 mya)
  ### Acol/Acep 7.05
  ### Atta/Acro 16.2
  ### Tsep/leafcutters 17.78
  ### Tcor/lc+Tsep 19.18
  ### Tzet/lc+Tset+Tcor 22.87
  ### Ccos/lc+Trachymyrmex 26.6

  # Ješovnik et al 2016 PlosOne
  ### Atta/Acro 19.9 (22.5 - 17.7)
  ### Tsep/leafcutters 22.3 (25.1 - 20.1)
  ### Tcor/lc+Tsep 29.2 (31.2 - 27.4)

  # Brandstetter et al 2017:
  ### Aech/Aoct 2.99
  ### Aech/Ahey 9.25
  ### Atta/Acro 16.93
  ### Tsep/leafcutters 19.31
  ### Tcor/lc+Tsep 21.59
  ### Tzet/lc+Tset+Tcor 27.29
  ### Ccos/lc+Trachymyrmex 34.59

  # Li et al 2018
  ### Aech/Aoct 2.26
  ### Aech/Ahey 6.29
  ### Atta/Acro 12.64
  ### Tsep/leafcutters 14.31
  ### Tcor/lc+Tsep 16.11
  ### Tzet/lc+Tset+Tcor 21.89 (fossil-calibrated)
  ### Ccos/lc+Trachymyrmex 28.73

  ### Minimum age bound for Tzet/lc+Tset+Tcor 15 (see supplement Table S7) or 22 (18-27) (see supplement Table S8)
  ### from Li et al 2018
  # Except for the root age, all four fossil calibration points were specified as a truncated Cauchy distribution indicated
  #  by L (tL, p, c), where tL = minimum-age bound (set as 15 Ma), P = offset value (default value of 0.1),
  #  and c = scale parameter value (default value of 1) representing a heavy-tailed density (50).

##################################################################
# First Calibration point: Split between Tzet and others
##################################################################

#################################
#> Thus, I chose for baseml control file a point calibration at 15 MA
#################################
(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)))),TZET)'@0.175',CCOS);
#################################
# So, the split is set to be 17.5 MA old, as Li wrote that the amber is between 15 MYA and 20 MYA.
# This strict setting is just a starting point for the prior.

#################################
#> Thus, I chose for mcmctree control file with a prior for the root at 'L(0.15, 0.1, 1, 0.1)'
#################################
(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)))),TZET)'L(0.15, 0.1, 1, 0.1)',CCOS);
#################################
# So, the root is set to be at least 15 MA, following the same distribution as Li used in their 2018 paper

##################################################################
# Second Calibration point: Atta/Acro split
##################################################################
# I will use the "S2N (skew 2 normals)" distribution (see PAML doc) for the second calibration point
#'SN2(p1, loc1, scale1, shape1, loc2, scale2, shape2)'

# First distribution includes 16.2, 19.9, 16.93, which are the estimates from Nygaard et al 2016, Jesovnik et al 2016 and Brandstetter et al 2017.
location = 17.68
scale = 2
shape = -1

# Second distribution includes the 12.64 from Li et al 2018
location = 12.64
scale = 2
shape = 1

p = value that defines which of the 2 normal distributions are more likely.
I set p=0.6
'SN2(p1, loc1, scale1, shape1, loc2, scale2, shape2)'


# plot in R
Rs
#plot(density(rskewnorm(1000000, location = 17.67667, scale = 2, shape = -3)),xlim=c(5,22))
#points(density(rskewnorm(1000000, location = 12.64, scale = 2, shape = 3)),type="l")
p<-0.6
hist(rskewnorm(p*1000000, location = 17.67667, scale = 2, shape = -1),xlim=c(5,22),200,col="red",border="red",main="",xlab="S2N(0.6,17.68,2,-1,12.64,2,1)")
hist(rskewnorm((1-p)*1000000, location = 12.64, scale = 2, shape = 1),add=T,200,col="red",border="red")

S2N(0.6,17.68,2,-1,12.64,2,1)

#################################
#> Thus, I chose for mcmctree control file with a prior for the Tzet split at 'L(0.15, 0.1, 1, 0.1)' and a
#>  prior following two skewed normal distributions capturing previous results
#################################
(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL))'S2N(0.6,17.68,2,-1,12.64,2,1)')),TZET)'L(0.15, 0.1, 1, 0.1)',CCOS)
#################################

#Input:
#alignment
#tree with lower bound calibration of root
#~/data/inqGen18/phylogeny/treeDating/MCMCtree/data/4dAll.topo.phy

#tree with range calibration prior distribution
#~/inqGen18/phylogeny/treeDating/MCMCtree/data/4dAll.mcmcTreeInput.phy


#################################################################################
#0. baseml
#################################################################################

# run baseml to get estimate for substitution rate per time unit
# For this, I need a specific calibration date. I used the mean for the Atta-Acro split as calculated from Haofu's alignment.
# Split between Acro and Atta at mean 16.2 mya

###############################
# Setup environment
###############################
base=~/data/inqGen18/phylogeny/treeDating/MCMCtree
bd=~/data/inqGen18/phylogeny/
software=/usr/local/home/lschrader/software
scripts=~/data/inqGen18/phylogeny/scripts/
mkdir $base
mkdir $base/data
cd $base/data
ln -s $bd/SCO/4d/4dAll.fa.phy .


#################################################################################
# 1st Run: 15.0 MY calibration point at Tzet/lc+Trachy split
##############################################################
mkdir $base/baseml1
cd $base/baseml1

SEQFILE=$base/data/4dAll.fa.phy
TOPOLOGYFILE=$base/data/4dAll.topo.15.0MYcalibration.phy

echo "11 1" > $TOPOLOGYFILE
echo "(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)))),TZET)'@0.15',CCOS);" >> $TOPOLOGYFILE

sed "s|<SEQFILE>|$SEQFILE|g" $base/ctl/baseml.template.ctl|sed "s|<TOPOLOGYFILE>|$TOPOLOGYFILE|g" > ./baseml1.ctl

baseml ./baseml1.ctl
#Check file mlb for substitution rate and alpha
#########
grep "Substitution rate is per time unit" -A 1 mlb > substitution.rate.txt
SR=$(sed 1d $base/baseml1/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
cat substitution.rate.txt
#Substitution rate is per time unit
#    0.460286 +- 0.001115

grep "alpha ("  mlb > alpha.txt
#Alpha
# alpha (gamma, K=5) =  1.47227
ALPHA=$(cat $base/baseml1/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')

#################################################################################

#################################################################################
# 2nd Run: 17.5 MY calibration point at Tzet/lc+Trachy split
##############################################################
mkdir $base/baseml2
cd $base/baseml2

SEQFILE=$base/data/4dAll.fa.phy
TOPOLOGYFILE=$base/data/4dAll.topo.17.5MYcalibration.phy

echo "11 1" > $TOPOLOGYFILE
echo "(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)))),TZET)'@0.175',CCOS);" >> $TOPOLOGYFILE


sed "s|<SEQFILE>|$SEQFILE|g" $base/ctl/baseml.template.ctl|sed "s|<TOPOLOGYFILE>|$TOPOLOGYFILE|g" > ./baseml1.ctl

baseml ./baseml1.ctl
#Check file mlb for substitution rate and alpha
#########
grep "Substitution rate is per time unit" -A 1 mlb > substitution.rate.txt
SR=$(sed 1d $base/baseml1/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
cat substitution.rate.txt
#Substitution rate is per time unit
#    0.394531 +- 0.000955

grep "alpha ("  mlb > alpha.txt
#Alpha
# alpha (gamma, K=5) =  1.47227
ALPHA=$(cat $base/baseml2/alpha.txt |perl -pe 's/.* = +(.*)$/$1/g')

#################################################################################

#################################################################################
# 3rd Run: 15-20 MY calibration point at Tzet/lc+Trachy split
##############################################################
mkdir $base/baseml3.fail
cd $base/baseml3.fail

SEQFILE=$base/data/4dAll.fa.phy
TOPOLOGYFILE=$base/data/4dAll.topo.15-20MYcalibration.phy

echo "11 1" > $TOPOLOGYFILE
echo "(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)))),TZET)''>0.15<0.2'',CCOS);" >> $TOPOLOGYFILE


sed "s|<SEQFILE>|$SEQFILE|g" $base/ctl/baseml.template.ctl|sed "s|<TOPOLOGYFILE>|$TOPOLOGYFILE|g" > ./baseml1.ctl

baseml ./baseml1.ctl
#Check file mlb for substitution rate and alpha
#########
grep "Substitution rate is per time unit" -A 1 mlb > substitution.rate.txt
SR=$(sed 1d $base/baseml1/substitution.rate.txt|perl -pe 's/  +(.*)? \+\-.*/$1/g')
cat substitution.rate.txt

## EMPTY! NOT WORKING!



#################################################################################
