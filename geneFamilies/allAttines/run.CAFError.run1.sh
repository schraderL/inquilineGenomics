#############################################
#Run1
#############################################

#!/corefac/cse/lukas/software/CAFE/release/cafe
# -filter = filter clusters that have an ancestral size of 0
# -p 0.05 = report gene families with a significant deviation from the other clusters (not important here)
# -t = number of threads

#load -i filtered_cafe_input.txt -t 4 -l reports/log_run6.txt
#load -i /usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.0.cafe.tsv -t 20 -l ErrSim1Run.txt -p 0.05 -filter
load -i /usr/local/home/lschrader/data/inqGen18/geneFamilies/results/mcl.clean/largeQC.TEpso.nofilter.MCL.I01.0.cafe.tsv -t 20 -l ErrSim1Run.txt
tree (CCOS:8.0,(TZET:7.0,(TCOR:6.0,(TSEP:5.0,((((ACHA:1.0,PARG:1.0)PA:1.0,AHEY:2.0)AH:1.0,(AECH:1.0,AINS:1.0)AI:2.0)AC:1.0,(ACEP:1.0,ACOL:1.0)AT:3.0)TS:1.0)TC:1.0)TZ:1.0)CC:1.0);


#Estimate lambda for all species using a single rate
# Error model estimate restricted to lambda. Mu is not included.
lambda -s -t (1,(1,(1,(1,((((1,1)1,1)1,(1,1)1)1,(1,1)1)1)1)1)1);
report reports/report_run1


#############################################
#Run1
#############################################
#Parg	0.28361328125
#10.4	0.23998046875
#Ains	0.1308984375
#5.9	0.23998046875
#Tsep	0.0
#2.6	0.2181640625
#Ahey	0.28361328125
#3.4	0.23998046875
#16.2	0.23998046875
#17.8	0.261796875
#Aech	0.10908203125
#Acha	0.261796875
#Acol	1.0
#This final error estimate has a score of: 133501.460892
#and a lambda of: 0.0207925496659
#############################################

#############################################
#Run2
#############################################
#Parg	0.333764648438
#10.4	0.262243652344
#Ains	0.214562988281
#5.9		0.286083984375
#Tsep	0.119201660156
#2.6		0.286083984375
#Ahey	0.357604980469
#3.4		0.262243652344
#16.2	0.262243652344
#17.8	0.286083984375
#Aech	0.166882324219
#Acha	0.309924316406
#Acol	0.0715209960938
#This final error estimate has a score of: 87865.408697
#and a lambda of: 0.00597081513724
#############################################
