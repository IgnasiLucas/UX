#!/bin/bash
#
#				2018-08-08
#				==========
#
# Once I can simulate populations of flies with specific demographic models
# and with mutations affecting male longevity in an equilibrium frequency,
# I would like to start simulating the experiment proposed by Kelly (1999).
#
# The python script should evolve the population to equilibrium and then apply
# the experiment, which consists on an inbreeding experiment, a short selection
# experiment, and a second inbreeding experiment. It will accept options to
# determine what kind of equilibrium to simulate, etc. For the moment, I will
# test Kelly's original idea, with only autosomal loci affecting a component
# of fitness. Once this works, I will test the experiment in the case of X-linked
# variation affecting male longevity.
#
