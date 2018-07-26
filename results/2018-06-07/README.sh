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
# Previous simulations showed that it is not easy to produce a noticeable
# difference in lifespan between males and females. One mistake was to have
# mutation effects too small. A wide range of values of the 'a' parameter
# makes gene effects be larger. Another mistake was to fix the relationship
# between 'a' and 'b'. The rationale was to prevent larvae from becoming
# targeted for natural death before reaching adulthood (10 days). However, it
# is clear that such misfortune does happen in nature. The model can stay
# simple by fixing the value of 'b'.
#
# Calculations of survival curves on 2018-06-21 suggest the following range
# of parameters: 0.003 < a < 0.0100, b ~ -0.019.

for XLoci in 20 100 200; do
   for ALoci in 0 20 100 200; do
      for mu in .000001 .00001 .0001 .001; do
         if [ ! -e MSB${XLoci}.${ALoci}${mu}.txt ]; then
            python mutationSelectionBalance.py \
               -m 0.003 -M 0.010 -N 50000 -G 500000 --step 100 \
               -X $XLoci -A $ALoci -u $mu -o MSB${XLoci}.${ALoci}${mu}.txt > MSB${XLoci}.${ALoci}${mu}.log &
         fi
      done
   done
done
wait

# The script mutationSelectionBalance.py wrote a first line for generation 0,
# which gnuplot complains about, because of the logaritmic scale of the X axis.
# I need to remove the first line.

for file in `ls -1 MSB*.txt`; do
   if head -n 1 $file | grep -q "^0"; then
      tail -n +2 $file > z1
      mv z1 $file
   fi
done

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
# I killed the simulations after 8 days and more than 200000 generations. The
# average mutant frequencies seemed to be at equilibrium, althogh the log-log
# plots do not make it apparent.
#
# 1. The combination of X-linked and autosomal genes affecting parameter 'a'
#    makes it difficult to interprete the results. I should have reported the
#    average allele frequencies for each kind of gene separately, because they
#    have different dynamics.
# 2. Only the highest mutation rates simulated (0.001 and 0.0001) produce significant
#    age differences between males and females. The difference is higher (~5 days)
#    when only X-linked genes are included.
# 3. With 200 genes in the X chromosome and a large enough number of autosomal
#    genes, the age difference drops after having reached a peak around generation
#    20000 or 30000. Why?
# 4. Equilibrium frequencies are higher when autosomal genes are included in the
#    simulation. The more autosomal genes, the higher the equilibrium frequency.
#
