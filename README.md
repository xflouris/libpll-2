# Libpll-2

libpll-2 is the new official fork of libpll (https://github.com/xflouris/libpll/). It implements site repeats to speed up computations.


Please read the wiki for more information.



# Projects that are already using libpll-2
 List of projects already using libpll-2 and site repeats, and reported speedups compared with the tip pattern optimization:
 * [RAxML-NG](https://github.com/amkozlov/raxml-ng): speedup ranges between 1.2 and 1.5 
 * [ModelTest-NG](https://github.com/ddarriba/modeltest): speedup around 2
 * [EPA-ng](https://github.com/Pbdas/epa-ng): no speedup (likelihood computation is not the main bottleneck) but memory footprint reduced by 30%.
