#!/bin/bash
#
#				2018-06-21
#				----------
#
# Here I implement a sexually antagonistic pleiotropy model, with one locus
# on the X chromosome affecting both male lifespan and female reproduction
# in a sexually antagonistic way. That is, the "1" allele gives females a
# selective advantage during reproduction. There is a coefficient of dominance
# and the fitness values are the following:
#
#                +--------------------+-------------+
#                |        Females     |    Males    |
#    +-----------+------+------+------+------+------+
#    | genotype  |  00  |  01  |  11  |  0Y  |  1Y  |
#    +-----------+------+------+------+------+------+
#    | fitness   | 1-s  | 1-hs | 1.0  | 1-s  | 1-s  |
#    +-----------+------+------+------+------+------+
#
# Where the selection coefficient "s" is set by the "-s" option, and the
# coefficient of dominance, "h" is set by the "-h" option.
#
# At the same time, allele "1" is deleterious for males by increasing their
# aging rate, "a", by a certain amount, which is set by the option "-e" in
# the SexChromSelectionBalance.py script.
#
# It is time to find an expression for the survival function under the smurf
# model of aging. In the simulations, aging, becoming a smurf or not, and dying
# or not all happen right before mating, in that order. To reach age "x" should
# mean to be able to mate at age "x" (if old enough to reproduce). Thus the
# probability of reaching age 1 is equal to the probability of not becoming a
# smurf or becoming a smurf and not dying at age 1. The probability of reaching
# age 2 is equal to the sum of the probabilities of three alternative events:
# becoming a smurf at age 1 and not dying neither at age 1 nor at age 2; not
# becoming a smurf at age 1, but at age 2 and not dying at age 2; and not becoming
# a smurf neither at age 1 nor at age 2.
#
# Let B(x) be the probability of becoming a smurf (blue) at age x, which is ax+b.
# The survival function among smurfs is an exponential decay, exp(-dx), where d
# is the constant death rate experienced by smurfs. In symbols, the probability
# of surviving to age x, S(x) is:
#
#          i=x              i=x                           j=i-1
# S(x) = PRODUCT (1-B(i)) + SUM [B(i) * exp((i-x-1)*d) * PRODUCT (1-B(j))]
#          i=1              i=1                            j=1
#
#
# See the formula in the survival.pdf file for a nicer format. I need the
# survival function to determine the fitness of each fenotype, in order to have
# a theoretical expectation of the allele frequency at equilibrium. I take
# fitness to be the lifetime expected reproduction, which in discrete ages is:
#
#       inf
#   R = SUM S(x)F(x)
#       x=1
#
# Where S(x) is the survival function, and F(x) the fecundity. In the model,
# adults (x >= 10) reproduce at a constant rate, limited by the constant population
# size. That is, there is no age-dependency in fecundity. Thus, F(x) can be
# substituted by the reproductive fitness, or probability of being chosen as
# a parent.
#
# The script survival.py outputs a list of survival probabilities for each
# combination of 'a' and 'b' values.

if [ ! -e survival.png ]; then
   if [ ! -e survival.txt ]; then
      python survival.py --min-a 0.001 --max-a 0.050 --min-b -0.010 --max-b -0.500 --min-d 0.1911 --max-d 0.1911
   fi
fi