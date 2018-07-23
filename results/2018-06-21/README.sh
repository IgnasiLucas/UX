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
#    | fitness   | 1-s  | 1-ds |   1  |   1  | 1-r  |
#    +-----------+------+------+------+------+------+
#
# Where the selection coefficient "s" is set by the "-s" option, and the
# coefficient of dominance, "d" is set by the "-d" option.
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

MIN_A=0.0025
MAX_A=0.0250
MEAN_A=0.0039
MIN_B=-0.030
MAX_B=-0.0001
MEAN_B=-0.019


if [ ! -e default.txt ]; then
   python survival.py -a 0.0039 -b -0.019 -d 0.1911 -n 1 -x 150 -o default.txt
fi

if [ ! -e constant_b.png ]; then
   if [ ! -e constant_b.txt ]; then
      python survival.py --min_a $MIN_A --max_a $MAX_A -b $MEAN_B -n 10 \
                         -x 150 --death_rate 0.1911 -o constant_b.txt
   fi
   gnuplot -e "min_a=$MIN_A; max_a=$MAX_A; b=$MEAN_B" plot_constant_b.gnp
fi

if [ ! -e constant_a.png ]; then
   if [ ! -e constant_a.txt ]; then
      python survival.py -a $MEAN_A --min_b $MIN_B --max_b $MAX_B -n 10 \
                         -x 150 --death_rate 0.1911 -o constant_a.txt
   fi
   gnuplot -e "min_b=$MIN_B; max_b=$MAX_B; a=$MEAN_A" plot_constant_a.gnp
fi

if [ ! -e constant_t.png ]; then
   if [ ! -e constant_t.txt ]; then
      python survival.py --min_a $MIN_A --max_a $MAX_A -t 4.87 -n 10 \
                         -x 150 --death_rate 0.1911 -o constant_t.txt
   fi
   gnuplot -e "min_a=$MIN_A; max_a=$MAX_A; t0=4.87" plot_constant_t.gnp
fi

# Say that the 'a' parameter of the two-stages aging model is 0.0039 for
# either females or males with the 0 allele. That is the empirical value
# estimated by Tricoire and Rera (2015). Allele 1 increses 'a'. Below, I
# represent the fitness and the selection coefficient against allele 1 as
# a function of the increase it produces in 'a'. I assume reasonable values
# of 'a' may go up to 0.025.

N=50
if [ ! -e fitness_a.png ]; then
   if [ ! -e fitness_a.txt ]; then
      if [ ! -e constant_b_$N.txt ]; then
         python survival.py --min_a $MEAN_A --max_a $MAX_A -b $MEAN_B -n $N \
                            -x 150 --death_rate 0.1911 -o constant_b_50.txt
      fi
      gawk -v SIZE=$N '(/^#a/){
         for (i=2; i<=NF; i++) {
            A[i] = $i
            W[$i] = 0
         }
      }(/^[^#]/){
         for (i=2; i<=NF; i++) {
            W[A[i]] += $i
         }
      }END{
         for (i=2; i<=SIZE + 1; i++) {
            print A[i] "\t" W[A[i]] / W[A[2]]
         }
      }' constant_b_50.txt > fitness_a.txt
   fi
   gnuplot plot_fitness_a.gnp
fi

# I see that fitness drops fast with small increments of 'a'. Now, we can
# choose a range for the effects of allele 1 on 'a' and a range for the fitness
# effects on females of allele '0', 's' and 'd', in the table below:
#
#                +--------------------+-------------+
#                |        Females     |    Males    |
#    +-----------+------+------+------+------+------+
#    | genotype  |  00  |  01  |  11  |  0Y  |  1Y  |
#    +-----------+------+------+------+------+------+
#    | fitness   | 1-s  | 1-ds |   1  |   1  | 1-r  |
#    +-----------+------+------+------+------+------+
#
# According to Fry (2009), who cites Rice (1984) and Hedrick (2000), the
# polimorphism will be stable if:
#
#        2·d         r      2·(1 - d)
#     ---------  <  --- <  -----------
#      1 + d·s       s       1 - d·s
#
# To make things simple, we can set 's' = 'r', and 'd' = 0.5, which does
# fulfill the stability condition. According to file fitness_a.txt, the
# following increments of 'a' ('e') produce such selection coefficients
# ('r').
#
#   +---------+----------+
#   |    e    |    r     |
#   +---------+----------+
#   | 0.0004  | 0.05336  |
#   | 0.0009  | 0.098107 |
#   | 0.0013  | 0.137243 |
#   | 0.0017  | 0.1711   |
#   | 0.0022  | 0.200737 |
#   | 0.0026  | 0.227182 |
#   | 0.003   | 0.251203 |
#   | 0.0034  | 0.272798 |
#   | 0.0039  | 0.292342 |
#   . ...     . ...      .
#
#
