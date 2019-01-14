#################################################################################
# TAG: Calculate divergence time with MCMC tree
#################################################################################

############################
# TAG: General overview
############################
: '
I use MCMCtree to infer species divergence dates. MCMCtree uses a bayesian approximate likelihood approach
to infer divergence dates using a calibrated topological tree as input. I here use the 4d site alignment to
infer divergence dates between the different attine species.
'

#####################
#TAG INPUT files
#####################

#alignment
: '~/data/inqGen18/phylogeny/SCO/4d/4dAll.fa.phy'

#tree with lower bound calibration of split at Tzet node. (for baseml)
: '~/data/inqGen18/phylogeny/treeDating/MCMCtree/data/4dAll.topo.17.5MYcalibration.phy'

#tree with soft range calibration prior distribution. (for mcmctree)
: '~/data/inqGen18/phylogeny/treeDating/MCMCtree/data/4dAll.mcmcTreeInput.phy'

# template control file for baseml
: '~/data/inqGen18/phylogeny/treeDating/MCMCtree/ctl/baseml.template.ctl'

# template control file for mcmctree
: '~/data/inqGen18/phylogeny/treeDating/MCMCtree/ctl/mcmctree.template.ctl'


############################
# TAG Tutorials & Resources
############################

# http://www.fish-evol.com/mcmctreeExampleVert6/text1Eng.html
# https://github.com/mariodosreis/divtime


############################
# TAG Calibration
############################
:'

Here is an overview of various different results for different splits in the attine phylogeny.

Nygaard et al 2011
  ## Atta-Acro split 10 MYA

Mikheyev et al 2010
  ## Atta-Acro split 10-16 or 8-12 MYA

Hoelldobler & Wilson The Leafcutter Ants: Civilization by Instinct:
  ## Atta-Acro split 8-12 MYA

Schultz & Brady 2008 used fossil data:
  ## Origin of leaf-cutting 8-12 MYA

Nygaard et al 2016:
  ## Acol/Acep 7.05
  ## Atta/Acro 16.2 [Nygaard et al 2016 Supplement, 16.2 mya (mean) 5%-95% (12.6 mya - 19.7 mya)]
  ## Tsep/leafcutters 17.78
  ## Tcor/lc+Tsep 19.18
  ## Tzet/lc+Tset+Tcor 22.87
  ## Ccos/lc+Trachymyrmex 26.6

Ješovnik et al 2016 PlosOne
  ## Atta/Acro 19.9 (22.5 - 17.7)
  ## Tsep/leafcutters 22.3 (25.1 - 20.1)
  ## Tcor/lc+Tsep 29.2 (31.2 - 27.4)

Brandstetter et al 2017:
  ## Aech/Aoct 2.99
  ## Aech/Ahey 9.25
  ## Atta/Acro 16.93
  ## Tsep/leafcutters 19.31
  ## Tcor/lc+Tsep 21.59
  ## Tzet/lc+Tset+Tcor 27.29
  ## Ccos/lc+Trachymyrmex 34.59

Li et al 2018
  ## Aech/Aoct 2.26
  ## Aech/Ahey 6.29
  ## Atta/Acro 12.64
  ## Tsep/leafcutters 14.31
  ## Tcor/lc+Tsep 16.11
  ## Tzet/lc+Tset+Tcor 21.89 (based on fossil-calibrated prior)
  ## Ccos/lc+Trachymyrmex 28.73

  ## Minimum age bound for Tzet/lc+Tset+Tcor 15 (see supplement Table S7) or 22 (18-27) (see supplement Table S8)
Li et al 2018
  Except for the root age, all four fossil calibration points were specified as a truncated Cauchy distribution indicated
   by L (tL, p, c), where tL = minimum-age bound (set as 15 Ma), P = offset value (default value of 0.1),
   and c = scale parameter value (default value of 1) representing a heavy-tailed density (50).


##################################################################
# TAG First Calibration: Point calibration for estimating substitution rate and alpha
##################################################################

Using baseml, I estimated an overall substitution rate and alpha that will subsequently
be used in the mcmctree runs. For this overall estimates, I need a point calibrated phylogeny as input.

Li et al 2018 have a fossil (N4) for the split between Tzet/lc+Tset+Tcor. The fossil is from Dominican amber,
which is estimated to be between 15 to 20 mya.

Li et al 2018:
  "the dating of Dominican amber is ambiguous, ranging from 15 to 20 Ma"
#> Thus, I chose for baseml control file a point calibration at 17.5 MA

#> My point calibrated phylogeny for baseml is:

(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)))),TZET)'@0.175',CCOS);

  # This strict setting is just a starting point for the prior.
'
# TAG Run baseML with the point calibrated phylogeny.
~/sciebo/inquilineGenomics18/phylogeny/TreeDating/MCMCtree/scripts/0baseml.sh


##################################################################
# TAG  Second Calibration: Soft calibration for split between Tzet and others
##################################################################
: '

Li et al 2018:
  "all four fossil calibration points were specified as a truncated Cauchy distribution indicated
  by L (tL, p, c), where tL = minimum-age bound (set as 15 Ma), P = offset value (default value of 0.1),
  and c = scale parameter value (default value of 1) representing a heavy-tailed density (50)."
#> Thus, I chose for mcmctree control file with a prior for the root at 'L(0.15, 0.1, 1, 0.1)'

#################################
(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL)))),TZET)'L(0.15, 0.1, 1, 0.1)',CCOS);
#################################
# So, the Tzet node is set to be at least 15 MA, following the same distribution as Li used in their 2018 paper.

##################################################################
# TAG Third Calibration: Soft calibration for Atta/Acro split.
##################################################################
I will use the "S2N (skew 2 normals)" distribution (see PAML doc) for the second calibration point.
This calibration is based on previous estimates of the divergence dates between Atta and Acro.
Some papers have suggested a split ~10 MYA (Nygaard et al 2011, Mikheyev et al 2010, Schultz & Brady 2008, Li et al 2018),
others suggest an older origin of ~16 MY (Nygaard et al 2016, Ješovnik et al 2016, Brandstetter et al 2017, Mikheyev et al 2010).

In mcmctree a "skew 2 normal" prior distribution can be set using the following tag in the input tree.
SN2(p1, loc1, scale1, shape1, loc2, scale2, shape2)

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
"SN2(p1, loc1, scale1, shape1, loc2, scale2, shape2)"

So, the final S2N distribution for the Atta/Acro split in MCMCtree is:

S2N(0.6,17.68,2,-1,12.64,2,1)
'

#To plot this S2N distribution in R, use:
:'
library(VGAM)
#plot(density(rskewnorm(1000000, location = 17.67667, scale = 2, shape = -3)),xlim=c(5,22))
#points(density(rskewnorm(1000000, location = 12.64, scale = 2, shape = 3)),type="l")
p<-0.1
plot(density(c(rskewnorm((1-p)*1000000, location = 12.64, scale = 2, shape = 1),rskewnorm(p*1000000, location = 17.67667, scale = 2, shape = -1))))
'
:'
#################################
Together, I chose for mcmctree control file with a prior for the Tzet split at "L(0.15, 0.1, 1, 0.1)" and a
prior following two skewed normal distributions for the Atta/Acro split capturing previous results: S2N(0.6,17.68,2,-1,12.64,2,1).
#################################

(((TCOR,(TSEP,(((AHEY,(ACHA,PARG)),(AECH,AINS)),(ACEP,ACOL))'S2N(0.6,17.68,2,-1,12.64,2,1)')),TZET)'L(0.15, 0.1, 1, 0.1)',CCOS)

'
##################################################################
# TAG Run MCMCtree
##################################################################

# two runs to see if they converge on the same posterior.
# two runs sampling only from the prior for comparison.
# See ./mcmctree.readme.txt for more details.
:'
Final settings used:
  sampFreq = 5000
  substitution_rate = 0.394531 estimated from @17.5 calibration
  sigma2_gamma=5.5
'

# RUN 1st run
/Users/lukas/sciebo/inquilineGenomics18/phylogeny/TreeDating/MCMCtree/scripts/mcmctree.run1.sh
# RUN 2nd run
/Users/lukas/sciebo/inquilineGenomics18/phylogeny/TreeDating/MCMCtree/scripts/mcmctree.run2.sh
# RUN 1st run sampling from prior only
/Users/lukas/sciebo/inquilineGenomics18/phylogeny/TreeDating/MCMCtree/scripts/p.mcmctree.run1.sh
# RUN 2nd run sampling from prior only
/Users/lukas/sciebo/inquilineGenomics18/phylogeny/TreeDating/MCMCtree/scripts/p.mcmctree.run2.sh

# RUN 3rd run with even longer chain
/Users/lukas/sciebo/inquilineGenomics18/phylogeny/TreeDating/MCMCtree/scripts/mcmctree.run3.sh
