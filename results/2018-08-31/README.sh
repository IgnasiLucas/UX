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
# fulfill the stability condition. According to file 2018-06-21/fitness_a.txt,
# the following increments of 'a' ('e') produce such selection coefficients:
# 'r1' for constant fecundity, and 'r2' for age-dependent fecundity according
# to: 1 - ((x-12)^2)/100.
#
#   +---------+----------+----------+
#   |    e    |    r1    |    r2    |
#   +---------+----------+----------+
#   | 0.0004  | 0.051169 | 0.020016 |
#   | 0.0009  | 0.095492 | 0.041183 |
#   | 0.0013  | 0.133012 | 0.061168 |
#   | 0.0017  | 0.165454 | 0.080400 |
#   | 0.0022  | 0.193837 | 0.098916 |
#   | 0.0026  | 0.220821 | 0.119469 |
#   | 0.0030  | 0.243837 | 0.137667 |
#   | 0.0034  | 0.264521 | 0.155195 |
#   | 0.0039  | 0.283232 | 0.172082 |
#   . ...     . ...      .

for i in `seq 1 20`; do
   if [ ! -e  `printf "sim%02i.txt" $i` ]; then
      python SexChromSelectionBalance.py -N 10000 -G 30000 -o `printf "sim%02i.txt" $i` > `printf "sim%02i.log" $i` &
   fi
done
wait
