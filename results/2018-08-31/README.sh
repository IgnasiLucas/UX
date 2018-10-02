#!/bin/bash
#
#				2018-08-31
#				==========
#
# Here I implement a sexually antagonistic pleiotropy model, with one locus
# on the X chromosome affecting both male lifespan and female reproduction
# in a sexually antagonistic way. That is, the "1" allele gives females a
# selective advantage during reproduction. There is a coefficient of dominance,
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
# coefficient of dominance, "d" is set by the "-d" option of the SexChromSelectionBalance.py
# script.
#
# At the same time, allele "1" is deleterious for males by increasing their
# aging rate, "a", by a certain amount, which is set by the option "-e" in
# the SexChromSelectionBalance.py script.
#
# According to Fry (2009), who cites Rice (1984) and Hedrick (2000), the
# polimorphism will be stable if:
#
#        2·d         r      2·(1 - d)
#     ---------  <  --- <  -----------
#      1 + d·s       s       1 - d·s
#
# To make things simple, we can set 's' = 'r', and 'd' = 0.5, which does
# fulfill the stability condition. According to the last results from 2018-06-21,
# an increase of 0.0013 of the 'a' parameter of the two-phases model, from
# 0.0039 to 0.0052 causes a reduction in relative fitness of 0.10827. Thus,
# using those values, the selective coefficient against the 0 allele in females
# should be 0.10827. This particular implementation may serve just to prove
# that a stable equilibrium can be reached, and that the model is well implemented.

for s in 0.108 0.109 0.110 0.112 0.114 0.116; do
   for i in 1 2 3 4; do
      if [ ! -e `LC_ALL=C printf "s_%.3f_sim%02i.txt" $s $i` ]; then
         python SexChromSelectionBalance.py -N 20000 -G 500000 \
                                            -e 0.0013 -s $s \
                                            -o `LC_ALL=C printf "s_%.3f_sim%02i.txt" $s $i` > `LC_ALL=C printf "s_%.3f_sim%02i.log" $s $i` &
      fi
   done
done
wait

if [ ! -e frequencies.png ]; then
   gnuplot plot_freqs.gnp
fi

if [ ! -e age_diff.png ]; then
   for s in 108 109 110 112 114 116; do
      for i in 1 2 3 4; do
         if [ ! -e smooth_agediff_$s\_$i.txt ]; then
            gawk '{
                     DIFF[$1 % 1000] = $4 - $3
                     if (length(DIFF) == 10) {
                        S = 0
                        for (diff in DIFF) {
                           S += DIFF[diff]
                        } print ($1 - 1000) + 500 "\t" S / 10
                     }
                  }' s_0.$s\_sim0$i.txt > smooth_agediff_$s\_$i.txt
         fi
      done
   done
   gnuplot plot_agediff.gnp
   rm smooth_agediff*
fi

# CONCLUSIONS
# ===========
#
# 1. Using the theoretical conditions of stability, the allele that benefits females
#    and damages males is lost. This means that the damage to males must actually be
#    larger than measured in 2018-06-21. A slightly higher selection coefficient against
#    the 'wild-type' is required in females to reach stability.
#
# 2. The 'stability' is actually very variable. Even in a population of 20000, fluctuations
#    of allele frequency are extraordinary.
#
# 3. The model is well implemented, in as much it does show a 'stable' equilibrium over
#    500000 days with 0.112 < s < 0.116. Actually, s > 0.116 seem to be also compatible
#    with the equilibrium.
#
# 4. Under polymorphism, females are on average 1 day older than males. Occasional peaks
#    of allele frequency cause the difference to reach up to 2 days or to drop below 0.
