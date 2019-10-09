          seed = -1
       seqfile = <SEQFILE>
      treefile = <MCMCtreeInput>
       outfile = out

         ndata = 1
       seqtype = 0  * 0: nucleotides; 1:codons; 2:AAs
       usedata = 3    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
         clock = 2    * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = <0.4  * safe constraint on root age, used if no fossil for root.

         model = 7    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 7: GTR+G
         alpha = <ALPHA>    * alpha for gamma rates at sites (collected from baseml results in file mlb)
         ncatG = 5    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0    * birth, death, sampling
   kappa_gamma = 6 2      * gamma prior for kappa. 
   alpha_gamma = 1 1      * gamma prior for alpha

   rgene_gamma = <RGENE1> <RGENE2>   * gamma prior for overall rates for genes. In the form rgene_gamma =  au bu a prior, where prior is either 0 (conditional iid) or 1 (gamma-dirichlet, default)
  sigma2_gamma = <SIGMA1> <SIGMA2>   * gamma prior for sigma^2     (for clock=2 or 3). In the form sigma2_gamma = au bu a prior, where prior is either 0 (conditional iid, default) or 1 (gamma-dirichlet)

     *finetune = 0: 0.06 0.5 0.008 0.12 0.4 * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr
     finetune = 1: 0.05  0.5  0.008  0.21 0.4  * auto (0 or 1) : times, musigma2, rates, mixing, paras, FossilErr

         print = 1
        burnin = 50000
      sampfreq = <SAMPFREQ>
       nsample = 20000

*** Note: Make your window wider (100 columns) before running the program.
