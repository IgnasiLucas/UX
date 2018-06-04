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
