#!/bin/bash
#
#				2018-05-22
#				==========
#
# I need to simulate several types of populations, and I want to go step
# by step. First, I should simulate a source population of flies with the
# following properties:
#   * At least 2000 individuals.
#   * With sexual conflict: at least one locus with opposite fitness effects
#     between sexes.
#   * Age-structured.
#   * Random mating.
#
# From the Wikipedia: "the last male to mate with a female sires about 80%
# of her offspring". They also say that females lay around 400 eggs.
#
# +-----------------------------------------+
# | Chromosome | Length (Mbp) | Length (cM) |
# +------------+--------------+-------------|
# |     X      |     22.4     |     ~75     |
# |     Y      |     40.0     |       0     |
# |     2      |     21.1     |    ~107     |
# |     3      |     27.9     |    ~110     |
# |     4      |      1.4     |       0     |
# +------------+--------------+-------------+
#
# After gaining some familiarity with simuPOP, I have a better idea of what is
# feasible. So far, I got an age-structured population which implements the
# model of the two stages of aging: adults with zero mortality rate and smurfs
# with a constant mortaility rate. Two parameters, a and b control the rate of
# transformation of adults in smurfs, which increases linearly with age. Those
# two parameters could be reduced to one, assuming that only adults (at least
# 10 days old flies) and not larvae can be targeted by natural death. The aging
# rate parameter, thus, could be affected by either several deleterious
# spontaneous mutations in mutation-selection equilibrium, or by sexually
# antagonistic variation, maintained in the X chromosome (see James D. Fry,
# 2009, Evolution 64-5:1510-1516.).
#
# Some details of the life cycle are difficult to implement and can be left
# aside for later. For example, the polygamy.
#
# Only after being able to generate these two types of source populations,
# I will be in conditions to calculate the quantitative genetics parameters
# that John Kelly attempts to estimate with his experiments, and simulate the
# experiments themselves.
#
# The discontinuous two-phases model of aging says that during the first phase
# the probability of transitioning into the second ("smurf") phase increases
# linearly with age: y = ax + b, where x is the age (in days) and -b/a is the
# age at which flies start to become smurfs. The second phase is described by
# an exponential survival function: S(x) = exp(-kx), where k is the rate constant.
# Thus, the probability of dying during a certain day x is 1 - exp(-k).

N=5000 # Population size.
k=0.4
if [ ! -e ageStructures.png ]; then
   for a in 0.05 0.04 0.03 0.02 0.01 0.005 0.002 0.001; do
      if [ ! -e ages_$a.txt ]; then
         b=$(echo "scale=4; -$a * 10" | bc -l)
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
