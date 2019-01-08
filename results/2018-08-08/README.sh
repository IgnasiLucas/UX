#!/bin/bash
#
#				2018-08-08
#				==========
#
# Once I can simulate populations of flies with specific demographic models
# and with mutations affecting male longevity in an equilibrium frequency,
# I would like to start simulating the experiment proposed by Kelly (1999).
#
# I need to split the process in several scripts, for the sake of modularity.
# First, I can use previous scripts to bring a population to an equilibrium,
# and then save it. The equilibrium will be either an antagonistic pleiotropy or
# a mutation-selection equilibrium. Then, the saved population will be the inpiut
# to the script that runs Kelly's experiment. Actually, the experiment can be split
# in two scripts, one for the selection steps and one for the inbreeding. A last
# script could be used to determine the statistics. All scripts should accept and/or
# produce populations, that need to be saved and loaded.
#
# 
