
We used CAFE to model gene family size evolution across the phylogeny of leaf-cutting ants. We only included Acromyrmex species and Atta colombica in the analysis as these genomes were sequenced using the same approach (long-range Illumina libraries). 


0. Combine the MCL clusterd gene families with the MCMCtree computed ultrametric tree
# run.geneFamilies.mcl.allAttines.sh

1. First run error estimation for each genome: run.CAFError.sh.
# run.CAFError.Acro.Acol.sh

2. Then run CAFE with 2 lambda/mu parameters on each species (run.CAFE.sh)
# run.CAFE.Acro.Acol.sh

3. Cluster these 2-parameter models to get a model that converges.
# plotCAFE.AcroAtta.Rmd

4. Run kmeans-clustered model using run.CAFE.kmeansModels.sh
# run.CAFE.kmeansModels.Acro.Acol.sh

5. Run LRT test to compare 1-parameter model with 5-parameter model
# run.CAFE.Acro.Acol.lrTests.sh
