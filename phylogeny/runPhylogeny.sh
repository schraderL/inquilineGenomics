
################################################################################
# 4-fold degenerate sites
################################################################################

################################################################################
# 3 Run phylogeny with 4d sites
################################################################################

base=~/data/inqGen18/phylogeny
data=~/data/inqGen18/phylogeny/data
software=/usr/local/home/lschrader/software
scripts=~/data/inqGen18/phylogeny/scripts/

cd $base/

################################################################################
# 3.1 Convert to different formats
################################################################################

cd $base/SCO/4d
# convert to sequential phylip (Fasta2Phylip.pl path/to/input/fastafile)
$software/Sequence-manipulation/Fasta2Phylip.pl $base/SCO/4d/4dAll.fa
# convert fa to nexus
prank -convert -d=./4dAll.fa -o=./4dAll -f=nexus

################################################################################
# 3.2 Run model selection with Jmodeltest
################################################################################
java -jar $software/jmodeltest-2.1.10/jModelTest.jar -d $base/SCO/4d/4dAll.fa -g 4 -i -f -AIC -BIC -a > $base/SCO/4d/jmodeltest.out
# GTR+I+G is the best fitting model


################################################################################
# 3.2 Run Fasttree
################################################################################
cd $base/fasttree
ln -s $base/SCO/4d/4dAll.fa.phy $base/fasttree
ln -s $base/SCO/4d/4dAll.fa $base/fasttree
#compute tree with FastTree
nice FastTree -nt -gtr < 4dAll.fa > 4d.fasttree.tre

################################################################################
# 3.3 Run Paup for starting tree
################################################################################
mkdir $base/Paup
ln -s $base/SCO/4d/4dAll.fa.phy $base/Paup
ln -s $base/SCO/4d/4dAll.fa $base/Paup
ln -s $base/SCO/4d/4dAll.nex $base/Paup

#compute starting tree with Paup (Maximum Parsomony cf Nygaard 2016)
#http://www.peter.unmack.net/molecular/programs/paup.command.blocks.html

cd $base/Paup
$software/PAUP/paup4a163_ubuntu64

#Parsimony Bootstrap Analysis
    execute 4dAll.nex
    set autoclose=yes;
    set criterion=parsimony;
    set root=outgroup;
    set storebrlens=yes;
    set increase=auto;
    outgroup Ccos;
    bootstrap nreps=1000 search=heuristic/ addseq=random nreps=10 swap=tbr hold=1;
    savetrees from=1 to=1 file=attines.1000bootstrap.MP.tree.nex format=altnex brlens=yes savebootp=NodeLabels MaxDecimals=0;

egrep "^tree" attines.1000bootstrap.MP.tree.nex |perl -pe 's/.*?(\(.*\;$)/$1/g' > attines.1000bootstrap.MP.tree.tre


################################################################################
# 3.4 Run PhyML
################################################################################

#run PhyML with b -1 option
# Run PhyML without aLRT tests

# -m GTR = GTR model
# -b -1 = approximate likelihood ratio test returning aLRT statistics.
# -s BEST = Can be either NNI (default, fast) or SPR (a bit slower than NNI) or BEST (best of NNI and SPR search).
# -c 8 = nb_subst_cat : number of relative substitution rate categories. Default : nb_subst_cat=4.
# -a e = gamma : distribution of the gamma distribution shape parameter. Can be a fixed positive value or e to get the maximum likelihood estimate.
# -u = user_tree_file : starting tree filename. The tree must be in Newick format.


# phyml1
# Simple run with approximate likelihood ratio tests and no bootstrapping
mkdir $base/phyml/1
cd $base/phyml/1
cp $base/SCO/4d/4dAll.fa.phy .
startingTree=$base/phyml/4d.fasttree.tre
$software/PhyML-3.1/PhyML-3.1_linux64 -i ./4dAll.fa.phy -d nt -m GTR -s BEST -b -1 -v 0 -c 8 -a e -u $startingTree
~/mpich-3.2/bin/mpirun -n 4 ~/software/PhyML-3.1/phyml-mpi -i ./4dAll.fa.phy -d nt -m GTR -s BEST -b -1 -v 0 -c 8 -a e -u ../4d.fasttree.tre

# phyml2
# Run with 10 bootstraps, runs already for 2 h.
mkdir $base/phyml/2
cd $base/phyml/2
cp $base/SCO/4d/4dAll.fa.phy .
#perl $folder/scripts/convertAln.pl ../4dAll.fa fasta 4dAll.phy phylip
nice mpirun -n 30 phyml-mpi -i 4dAll.fa.phy -d nt -m GTR -s BEST -b 2 -v 0 -c 8 -a e > phyml2.log
# killed after >90 h of runtime without progress

#phyml 3
# PhyML with -b -3 option and Paup MP Starting Tree
mkdir $base/phyml/3
cd $base/phyml/3
cp $base/SCO/4d/4dAll.fa.phy .
startingTree=$base/phyml/Paup/attines.1000bootstrap.MP.tree.tre
$software/PhyML-3.1/PhyML-3.1_linux64 -i ./4dAll.fa.phy -d nt -m GTR -s BEST -b -3 -v 0 -c 8 -a e -u $startingTree

# phyml 4
mkdir $base/phyml/4
cd $base/phyml/4
cp $base/SCO/4d/4dAll.fa.phy .
startingTree=$base/phyml/Paup/attines.1000bootstrap.MP.tree.tre
nice mpirun -n 10 phyml-mpi -i ./4dAll.fa.phy -d nt -m GTR -s BEST -b 10 -v 0 -c 4 -a e -u $startingTree
# crashes with segmentation fault
$software/PhyML-3.1/PhyML-3.1_linux64 -i ./4dAll.fa.phy -d nt -m GTR -s BEST -b 10 -v 0 -c 4 -a e -u $startingTree
#Number of substitution rate categories
#The default is having all the sites evolving at the same rate, hence having one substitution rate category. A discrete-gamma distribution can be used to account for variable substitution rates among sites, in which case the number of categories that defines this distribution is supplied by the user. The higher this number, the better is the goodness-of-fit regarding the continuous distribution. The default is to use four categories, in this case the likelihood of the phylogeny at one site is #averaged over four conditional likelihoods corresponding to four rates and the computation of the likelihood is four times slower than with a unique rate. Number of categories less than four or higher than eight are not recommended. In the first case, the discrete distribution is a poor approximation of the continuous one. In the second case, the computational burden becomes high and an higher number of categories is not likely to enhance the accuracy of phylogeny estimation.


################################################################################
# 3.5 Run phylogeny with MrBayes
################################################################################
base=~/data/inqGen18/phylogeny
mkdir $base/mrbayes
mkdir $base/mrbayes/1

# convert fa to nex
cd $base/SCO/4d/
#prank -convert -d=./4dAll.fa -o=./4dAll -f=nexus
perl $base/scripts/convertAln.pl 4dAll.fa fasta 4dAll.nex nexus

#http://mrbayes.sourceforge.net/wiki/index.php/Tutorial_3.2
#At the end of the run, MrBayes asks whether or not you want to continue with the analysis.
#Before answering that question, examine the average standard deviation of split frequencies.
#As the two runs converge onto the stationary distribution, we expect the average standard
#deviation of split frequencies to approach zero, reflecting the fact that the two tree samples
#become increasingly similar. In our case, the average standard deviation is about 0.07 after
#1,000 generations and then drops to less than 0.000001 towards the end of the run.
#Your values can differ slightly because of stochastic effects. Given the extremely
#low value of the average standard deviation at the end of the run, there appears to
# be no need to continue the analysis beyond 10,000 generations so when MrBayes asks
# "Continue with analysis? (yes/no):", stop the analysis by typing no.

#By default, MrBayes will run two simultaneous, completely independent analyses
#starting from different random trees (Nruns = 2). Running more than one analysis
#simultaneously is very helpful in determining when you have a good sample from
#the posterior probability distribution, so we suggest that you leave this setting as is.

#http://mrbayes.sourceforge.net/wiki/index.php/Tutorial_3.2
cd $base/mrbayes/1
ln -s $base/SCO/4d/4dAll.nex .
mb
> execute 4dAll.nex
# We want the GTR model (Nst=6) with gamma distribution and a proportion of invariant sites(= GTR + G + I)
> lset nst=6 rates=invgamma
> mcmc ngen = 200000 samplefreq = 10
# If the standard deviation of split frequencies is below 0.01 after 100,000 generations, stop the run by answering no when the program asks "Continue the analysis? (yes/no)". Otherwise, keep adding generations until the value falls below 0.01.
> help lset
> showmodel
# Summarize the parameter values by typing sump burnin=250 (or whatever value corresponds to 25 % of your samples)
> sumt burnin = 50000
# Check the output and make sure that the potential scale reduction factor (PSRF) is reasonably close to 1.0 for all parameters (ideally below 1.02); if not, you need to run the analysis longer.

# summarize parameter values
> sump burnin = 50000

mb runMrB.txt > log.txt
# begin mrbayes;
#   set autoclose=yes nowarn=yes;
#   execute 4dAll.nex;
#   lset nst=6 rates=invgamma;
#   mcmc nruns=2 ngen=200000 samplefreq=10 file=4dAll.nex1;
#   mcmc file=4dAll.nex2;
#   mcmc file=4dAll.nex3;
# end;


cd $base/mrbayes/2
ln -s $base/SCO/4d/4dAll.nex .
mb
> execute 4dAll.nex
# We want the GTR model (Nst=6) with gamma distribution and a proportion of invariant sites(= GTR + G + I)
> lset nst=6 rates=invgamma
> help lset
> showmodel
> mcmc ngen = 200000 samplefreq = 10
# If the standard deviation of split frequencies is below 0.01 after 100,000 generations, stop the run by answering no when the program asks "Continue the analysis? (yes/no)". Otherwise, keep adding generations until the value falls below 0.01.
> sumt burnin = 50000
# Check the output and make sure that the potential scale reduction factor (PSRF) is reasonably close to 1.0 for all parameters (ideally below 1.02); if not, you need to run the analysis longer.
# summarize parameter values
> sump burnin = 50000

#mb runMrB.txt > log.txt
## begin mrbayes;
##   set autoclose=yes nowarn=yes;
##   execute 4dAll.nex;
##   lset nst=6 rates=invgamma;
##   mcmc nruns=2 ngen=20000 samplefreq=10 file=4dAll.nex1;
##   mcmc file=4dAll.nex2;
##   mcmc file=4dAll.nex3;
## end;

ln -s 4dAll.nex.con.tre 4dAll.mrbayes.tre

################################################################################
# 3.6 Run RaxML
################################################################################
mkdir $base/raxml
mkdir $base/raxml/1
cd $base/raxml/1

ln -s $base/SCO/4d/4dAll.fa.phy .

#nice raxmlHPC-PTHREADS –f i -s 4dAll.fa.phy -n boot -m GTRGAMMA -b 12345 -p 123456 -# 2

# run with RAxML version 8.2.8 released by Alexandros Stamatakis on March 23 2016.
nice raxmlHPC-PTHREADS -f a -m GTRGAMMA -p 12345 -x 12345 -# 500 -s 4dAll.fa.phy -n run1 -T 20

###########################
# Run with rapid bootstraps
# Run GTRGAMMAI
###########################
mkdir $base/raxml/2
cd $base/raxml/2
# run with RAxML version 8.2.12 released by Alexandros Stamatakis on May 2018.
ln -s $base/SCO/4d/4dAll.fa.phy .
#https://mctavishlab.github.io/MLsearchLab
#To set the evolutionary model to GTR+I+G, as we did in the Garli run, we will use '-m GTRGAMMAI'. RAxML also implements a similar, faster model, GTRCAT, that uses a different approximation to capture rate heterogeneity across sites. GTRCAT is faster and can yield trees with slightly better likelihood values (see Stamatakis 2006). It is not a good idea to use the CAT approximation of rate heterogeneity on datasets with less than 50 taxa. In general there will not be enough data per alignment column available to reliably estimate the per-site rate parameters. We only have 29 taxa here, so we will stick with GTRGAMMA.
nice raxmlHPC-PTHREADS -f a -m GTRGAMMAI -p 12345 -x 12345 -# 500 -s 4dAll.fa.phy -n T20 -T 20

###########################
# Run with rapid bootstraps
# Run GTRGAMMA
###########################
mkdir $base/raxml/3
cd $base/raxml/3
# run with RAxML version 8.2.12 released by Alexandros Stamatakis on May 2018.
ln -s $base/SCO/4d/4dAll.fa.phy .
nice raxmlHPC-PTHREADS -f a -m GTRGAMMA -p 12345 -x 12345 -# 500 -s 4dAll.fa.phy -n T20 -T 20

###########################
# Run 50 ML trees with 100 regular bootstraps
# Run GTRGAMMA
###########################
mkdir $base/raxml/4
cd $base/raxml/4
ln -s $base/SCO/4d/4dAll.fa.phy .
nice raxmlHPC-PTHREADS -m GTRGAMMA -p 12345 -# 50 -s 4dAll.fa.phy -n 4 -T 10
nice raxmlHPC-PTHREADS -m GTRGAMMA -p 12345 -b 12345 -# 100 -s 4dAll.fa.phy -n 4boot -T 10
nice raxmlHPC-PTHREADS -m GTRCAT -p 12345 -f b -t RAxML_bestTree.4 -z RAxML_bootstrap.4boot -n 4bp
nice raxmlHPC-PTHREADS -m GTRCAT -J STRICT -z RAxML_bootstrap.4boot -n 4strict
#https://sco.h-its.org/exelixis/web/software/raxml/hands_on.html

###########################
# Run 50 ML trees with 500 regular bootstraps
# Run GTRCAT
# Use Ccos as outgroup
###########################
mkdir $base/raxml/5
cd $base/raxml/5
ln -s $base/SCO/4d/4dAll.fa.phy .

#http://evomics.org/learning/phylogenetics/raxml/
# don't use I according to Stamakakis
# GTRCAT and GTRGAMMA very similar, but CAT much faster
nice raxmlHPC-PTHREADS -m GTRCAT -p 12345 -# 50 -s 4dAll.fa.phy -n 5 -T 10 -o CCOS
nice raxmlHPC-PTHREADS -m GTRCAT -p 12345 -b 12345 -# 500 -s 4dAll.fa.phy -n 5boot500 -T 10 -o CCOS
nice raxmlHPC-PTHREADS -m GTRCAT -p 12345 -f b -t RAxML_bestTree.5 -z RAxML_bootstrap.5boot500 -n 5bp
nice raxmlHPC-PTHREADS -m GTRCAT -J STRICT -z RAxML_bootstrap.5boot500 -n 5strict
nice raxmlHPC-PTHREADS -m GTRCAT -J MR -z RAxML_bootstrap.5boot500 -n 5majrule
nice raxmlHPC-PTHREADS -m GTRCAT -J MRE -z RAxML_bootstrap.5boot500 -n 5extmajrule
#https://sco.h-its.org/exelixis/web/software/raxml/hands_on.html

#see /corefac/cse/lukas/inqGen18/phylogeny/raxml/5/RAxML_bipartitionsBranchLabels.5bp
