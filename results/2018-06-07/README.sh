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
      for mu in .00001 .0001 .001 .005 .01; do
         if [ ! -e MSB$min_a$max_a$mu.txt ]; then
            python mutationSelectionBalance.py \
               -m $min_a -M $max_a -u $mu -N 2000 -G 500000 \
               -o MSB$min_a$max_a$mu.txt > MSB$min_a$max_a$mu.log &
         fi
      done
   done
done
wait

if [ ! -e A_difference.png ]; then
   gnuplot < plotA.gnp
fi

if [ ! -e age_difference.png ]; then
   gnuplot < plotAges.gnp
fi

if [ ! -e frequencies.png ]; then
   if [ ! -e frequencies.txt ]; then
      echo "#MIN_a" > frequencies.txt
      echo "#MAX_a" >> frequencies.txt
      echo "#MU"    >> frequencies.txt
      cut -f 1 MSB.003.05.00001.txt >> frequencies.txt
      for min_a in .003 .005 .010 .015 .020; do
         for max_a in .05 .10 .20 .40; do
            for mu in .00001 .0001 .001 .005 .01; do
               gawk -v MIN_A=$min_a -v MAX_A=$max_a -v MU=$mu 'BEGIN{
                  print MIN_A "\n" MAX_A "\n" MU
               }{
                  S=0
                  for (i=2; i<=29; i++) S += $i
                  printf("%.4f\n", S/28)
               }' MSB$min_a$max_a$mu.txt > z1
               paste frequencies.txt z1 > z2
               mv z2 frequencies.txt
               rm z1
            done
         done
      done
   fi
   gnuplot < plotFrequencies.gnp
fi

# CONCLUSIONS
# ===========
#
# Under mutation-selection balance, mutant allele frequencies vary along the
# generations, due to drift and mutation. Natural selection is better able to
# keep frequencies low when mutation rate is lower and the gene effect on aging
# rate is higher, as expected. This happens when the difference between the
# minimum and the maximum values of the 'a' parameter are further appart, which
# makes sense, given the model of gene effects. Under low effect and high mutation
# rate, eventually all loci become fixed for the mutant allele. At that point,
# the rate of aging becomes constant and equal for males and females.
#
# The rate of aging, 'a', ends up differring between males and females very little,
# and only if mutation rates are relatively high. Of course, males tend to have a
# higher value of 'a', because of the heterozygosity on the X chromosome, as
# expected. However, the differences in 'a' are not reflected on the maximum
# lifespan. This may be due to a bad choice of the quantile (100th) to represent
# the age distribution in males and females.
#
# In any case, it is clear that the effect of the heterozygosity at X-linked loci
# is quite mild. In what follows, I should limit simulations to high gene effects,
# relatively high mutation rates, and probably a higher proportion of QTL on the
# X chromosome to represent the mutation-selection balance.

