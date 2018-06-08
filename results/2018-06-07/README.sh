#!/bin/bash
#
#				2018-06-07
#				==========
#
# I wrote the script mutationSelectionBalance.py, which simulates mutations
# on loci that affect the probability of becoming targeted for natural
# death. Mutant alleles all increase the parameter 'a' of the two-phases
# model of aging, and are therefore deleterious. They are recessive, and
# expressed in males if present in the X chromosome. The model includes
# 28 loci, with recombination. The gene effects are all equal and additive.
#
# The challenge is to adjust the parameter values to reproduce some real-life
# features. Mainly, a low but positive equilibrium allele frequency, and a
# longer lifespan for females (unguarded X hypothesis).
#
# The number of loci should be more than enough. The gene effects should be
# small for the equilibrium frequencies to be positive, but large enough to
# cause a difference in lifespan between males and females.

for min_a in .003 .005 .010 .015 .020; do
   for max_a in .05 .10 .20 .40; do
      for mu in .00001 .0001 .001 .01; do
         if [ ! -e MSB_$min_a.$max_a.$mu.txt ]; then
            python mutationSelectionBalance.py \
               -m $min_a -M $max_a -u $mu -N 2000 -G 10000 \
               -o MSB_$min_a.$max_a.$mu.txt > MSB_$min_a.$max_a.$mu.log
         fi
      done
   done
done
