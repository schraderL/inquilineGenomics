#!cafe
# nice /corefac/cse/lukas/software/CAFE/release/cafe
#load -i filtered_cafe_input.txt -t 4 -l reports/log_run4.txt
# mkdir /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/lhtest/
load -i /usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.5.cafe.filter100.tsv -t 20 -l /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/lhtest/lhtest.out
tree (((AHEY:0024.983,(ACHA:0016.335,PARG:0016.335):0008.648):0027.567,(AECH:0009.605,AINS:0009.605):0042.945):0065.115,ACOL:0117.665);

#assign error models

errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp/cafe_errormodel_0.0526635742187.txt -sp AECH
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp/cafe_errormodel_0.0478759765625.txt -sp AHEY
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp/cafe_errormodel_0.0430883789063.txt -sp AINS
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp/cafe_errormodel_0.0287255859375.txt -sp ACHA
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp/cafe_errormodel_0.0.txt -sp ACOL
errormodel -model /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/CAFE.errorprediction/errorTest.smallSet.filter100/tmp/cafe_errormodel_0.129265136719.txt -sp PARG

# set null model lambda estimate
lambda -l 0.000574758
# could also use lambdamu
## lambdamu -l 0.000574758 -m 0.00041658305593
#The product of λ or μ and the depth of the tree should not exceed one (i.e., λ × t < 1 and μ × t < 1 must be true; where t is the time from tips to the root). See the Known limitations section for details as to how to diagnose this problem.

# create randomized gene set
genfamily /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/lhtest/rnd -t 100
# test against best clustered tree
lhtest -d /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/lhtest/ -t (((3,(3,1)5)2,(4,5)2)2,2) -l 0.000574758 -o /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/lhtest/lhtest_result.txt


# cut -f 2,4 /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/tests/lhtest_result.txt > /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/tests/lk_diffs.txt
# cd /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/tests/
# Rscript lhtest.R lk_diffs.txt -27399.1 -27947.2

# testing
# create randomized gene set
#genfamily /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/lhtest/rnd -t 1
# test against best clustered tree
#lhtest -d /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/lhtest/ -t (((3,(3,1)5)2,(4,5)2)2,2) -m 0.00041658305593 -l 0.000574758 -o /usr/local/home/lschrader/data/inqGen18/geneFamilies/CAFE/lhtest/lhtest_result.tmp.txt
