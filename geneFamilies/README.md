
We used CAFE to model gene family size evolution across the phylogeny of leaf-cutting ants. We only included Acromyrmex species and Atta colombica in the analysis as these genomes were sequenced using the same approach (long-range Illumina libraries).


0. Combine the MCL clustered gene families with the MCMCtree computed ultrametric tree
    ```
    run.geneFamilies.mcl.sh
    ```
1. First run error estimation for each genome
    ```
    run.CAFError.sh
    ```
2. Then run CAFE with 2 lambda/mu parameters on each species
    ```
    run.CAFE.sh
    ```
3. Cluster these 2-parameter models to find models that can converge.
    ```
    analyseCAFEmodels.Rmd
    ```
4. Run LRT tests to compare models
    ```
    CAFE.lhtest.R
    ```
5. Analyse final CAFE output
    ```
    analyseCAFEoutput.Rmd
    ```
