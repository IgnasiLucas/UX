#!/bin/bash
#
#				2018-06-21
#				----------
#
# Here I implement a sexually antagonistic pleiotropy model, with one locus
# on the X chromosome affecting both male lifespan and female reproduction
# in a sexually antagonistic way. That is, the "1" allele gives females a
# selective advantage during reproduction. There is a coefficient of dominance
# and that the fitness values are the following:
#
#                +--------------------+-------------+
#                |        Females     |    Males    |
#    +-----------+------+------+------+------+------+
#    | genotype  |  00  |  01  |  11  |  0Y  |  1Y  |
#    +-----------+------+------+------+------+------+
#    | fitness   | 1-s  | 1-hs | 1.0  | 1-s  | 1-s  |
#    +-----------+------+------+------+------+------+
#
# Where the selection coefficient "s" is set by the "-f" option, and the
# coefficient of dominance, "h" is set by the "-d" option.
#
# At the same time, allele "1" is deleterious for males by increasing their
# aging rate, "a", by a certain amount, which is set by the option "-e" in
# the SexChromSelectionBalance.py script.
#
# It is time to find an expression for the survival function under the smurf
# model of aging. In the simulations, aging, becoming a smurf or not, and dying
# or not all happen right before mating, in that order. To reach age "x"
# means to be able to mate at age "x".
