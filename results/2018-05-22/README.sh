#!/bin/bash
#
#				2018-05-22
#				==========
#
# This is the first attempt to implement a demographic model of Drosophila
# melanogaster with the Python package simuPOP. It is the basis for the
# simulations of the unguarded X hypothesis and its alternative, the sexually
# antagonistic pleiotropy model (see 2018-06-07 and 2018-06-21, respectively).
#
# Here, only a two-phases aging model is implemented (Tricoire and Rera, 2015).
# See 2018-07-20 for the classic Gompertz and Weibull models. Up to now, the
# goal is simply to convince myself that the demographic model works as
# expected.
#
# Some simplifications of the demographic model needs justification. First,
# the number of chromosomes does not have to be larger than two: an X chromosome
# and one autosome. The distribution of genes affecting the parameters of the
# aging model among the two types of chromosomes is to be studied in the future,
# and does not need to be fixed. Thus, I purposefully ignore the lengths of
# the real Drosophila chromosomes.
#
# Second, from the Wikipedia: "the last male to mate with a female sires about
# 80% of her offspring". They also say that females lay around 400 eggs. However,
# simuPOP populates the next generation up to a certain population size. There
# is no need to represent fecundity in accurate terms when keeping population
# size constant. Setting the number of offsprings per mating equal to 1 maximizes
# the effective population size, and therefore, the effect of natural selection.
#
# In the two-stages model of aging, there are two types of individuals: normal
# and those targeted to die soon. The latter are called 'smurfs' in the literature.
# Normal individuals do not die, while smurfs have a constant mortality rate
# (exponential model of mortality). Normal individuals randomly become smurfs
# at a rate linearly increasing with age since the first age at which they can
# become smurfs. That is, the hazard function of the event of becoming a smurf
# is 'a·x + b', where 'x' is the age, and 't0 = -b/a' is the first age at which
# the fly can become a smurf.
#
# To implement this model, instead of the rate we need the actual probability of
# becoming a smurf at a certain moment of the simulation. Time is continuous, even
# though simulations proceed stepwise, in discrete cycles. Thus, at the time the
# probability of death is applied in the simulation, it should be equivalent to
# the total probability of death along one cycle, which is one day, the unit of
# time in Drosophila models of aging. Depending on the sequence of events in the
# cycle, the probability of death is calculated differently. There are three events
# in the life cycle and two possible orders: aging -> reproducing -> surviving (ARS),
# or aging -> surviving -> reproducing (ASR). In the two phases model, becoming a
# smurf is just part of the chance to die. The meaning of surviving to age 'x'
# is to be able to reproduce at age 'x'. In an ASR cycle, those who survive also
# reproduce before further aging. Thus, survival to age 'x' is challenged when age
# is already set to 'x'. Thus, the probability of death is that of dying between
# 'x-1' (last time they reproduced) and 'x'.
#
# Conversely, in an ARS cycle, death is a filter to aging. When survival is challenged
# at age 'x', they had already reproduced at that age, and it is reproduction at
# age 'x+1' that is challenged. Thus, the probability of death applied at age 'x'
# should represent the chance of death at any point between 'x' and 'x+1'.
#
# Finally, I decided for an ARS cycle, where survival is challenged between
# reproduction and aging. In the two-phases model this involves a chance to become
# a smurf, if it is not one yet, and then a chance of dying if one is a smurf. The
# probability of becoming a smurf at any point between ages 'x' and 'x+1', conditional
# on not having become one before age 'x' is expressed as follows in terms of the
# function N(t) = exp(-a·(t-t0)^2 / 2), which expresses the probability of not having
# become a smurf by age 't':
#
#                              N(t) - N(t+1)
#    P(t < T < t+1 | T > t) = --------------- = 1 - exp( -a·t + a·t0 - a/2)
#                                  N(t)
#
# Then, if the individual is a smurf, the probability of dying between 'x' and
# 'x+1' is calculated similarly, from an exponential survival function S(t)=exp(-r·t):
# Note that 't' is time since having become a smurf. Luckily, the hazard function
# is constant, and I don't need to keep track of when an individual became a
# smurf:
#
#                              S(t) - S(t+1)
#    P(t < T < t+1 | T > t) = --------------- = 1 - exp(-r)
#                                  S(t)
#
# Of course, there is a small error here: the day a fly turns into a smurf,
# I will be killing it more probably than I should, because it could have
# turned into a smurf just before 'x+1', and I just assume it did right
# after 'x'.
#
# Only after being able to generate these two types of source populations,
# I will be in conditions to calculate the quantitative genetics parameters
# that John Kelly attempts to estimate with his experiments, and simulate the
# experiments themselves.
#

N=5000 # Population size.
k=0.4
if [ ! -e ageStructures.png ]; then
   b=-0.019
   for a in 0.0025 0.0050 0.0075 0.0100 0.0125 0.0150 0.0175 0.0200; do
      if [ ! -e ages_$a.txt ]; then
         #b=$(echo "scale=4; -$a * 10" | bc -l)
         # Python 3 is assumed.
         python age_structure.py -a $a -N $N -b $b -G 300 -o ages_$a.txt 1> ages_$a.log
      fi
   done
   if [ ! -e default.txt ]; then
      python age_structure.py -N $N -G 300 -o default.txt 1> default.log
   fi
   gnuplot < plotAgeStructures.gnp
fi

if [ ! -e default.png ]; then
   gnuplot < plotDefault.gnp
fi

# See Tricoire and Rera, 2015, 'A new, discontinuous 2 phases of aging model: lessons
# from Drosophila melanogaster', PLoS ONE ( https://doi.org/10.1371/journal.pone.0141920).
# The authors obtain the following estimates for the parameters:
#
# k = 0.1911 (IC_(95) [0.1694, 0.2129]
# a = 0.0039
# b = -0.019
#
# This implies 5-days old larvae to start becoming smurfs. This proved difficult to
# model with simuPOP, because the age-structured population is obtained by 'cloning'
# all individuals every generation, and then producing additional offspring through
# mating. But in a cohort, if they start dying before reaching adulthood, they cannot
# keep population size constant. I added a function to provide the right population
# size each generation.
