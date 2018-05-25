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
#
# 
