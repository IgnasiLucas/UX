#!/bin/bash
#
#				2018-06-07
#				==========
#
# I wrote the script mutationSelectionBalance.py, which simulates mutations
# on loci that affect the probability of becoming targeted for natural
# death. Mutant alleles all increase the parameter 'a' of the two-phases
# model of aging, and are therefore deleterious. They are recessive, and
# expressed in males if present in the X chromosome. The gene effects are
# all equal and additive.
#
# The challenge is to adjust the parameter values to reproduce some real-life
# features. Mainly, a low but positive equilibrium allele frequency, and a
# longer lifespan for females (unguarded X hypothesis).
#
# Some preliminar simulations suggested that the largest effect of the mutations
# produced a larger, more realistic difference between males and females. Thus,
# I let the 'a' parameter to vary between 0.003 and 0.400, which is the widest
# range tested so far. I will focus on the effect of the number of genes, and its
# distribution among chromosome types.

for XLoci in 20 100 200; do
   for ALoci in 0 20 100 200 1000; do
      for mu in .000001 .00001 .0001 .001; do
         python mutationSelectionBalance.py \
            -m 0.003 -M 0.400 -N 50000 -G 500000 --step 100 \
            -X $XLoci -A $ALoci -u $mu -o MSB${XLoci}.${ALoci}${mu}.txt > MSB${XLoci}.${ALoci}${mu}.log &
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

