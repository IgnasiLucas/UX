###############################################################
#                           MODULES                           #
###############################################################

import simuPOP as sim
import random
import math
import argparse
import re

###############################################################
#                 ARGUMENTS AND  VARIABLES                    #
###############################################################

parser = argparse.ArgumentParser(description = "Evolves a population to equilibrium and then runs Kelly's experiment.")
parser.add_argument('-d', '--demography', default='two_phases', type=str, help='Aging model: "two_phases", "gompertz", or "weibull". Default: "two_phases".')
parser.add_argument('-q', '--equilibrium', default='mutationselection', type=str, help='Type of equilibrium: "mutation_selection" or "antagonistic_pleiotropy".')
parser.add_argument('-a', default=0.0039, type=float, help='a parameter. Default: 0.0039.')
parser.add_argument('-b', default=-0.019, type=float, help='b parameter. Default: -0.019')
parser.add_argument('-k', default=0.1911, type=float, help='rate of mortality of smurfs, for two-phases model. Default: 0.1911')

parser.add_argument('-e', '--meffect', type=float, default=0.1, help='Male-specific selective coefficient against mutant allele in antagonistic pleiotropy model. Default: 0.1.')
parser.add_argument('-s', '--feffect', type=float, default=0.1, help='Female selective coefficient against wild type allele in antagonistic pleiotropy model. Default: 0.1.')
parser.add_argument('-d', '--female-dominance', type=float, default=0.5, help='Coefficient of dominance of deleterious allele in females in antagonistic pleiotropy model. Default: 0.5')
parser.add_argument('-H', '--male-dominance', type=float, default=0.5, help='Coefficient of dominance of deleterious allele in males in antagonistic pleiotropy model. Default: 0.5')

parser.add_argument('-M', default=0.0500, type=float, help='Maximum value of the "a" parameter in the mutation-selection balance model. Default: 0.050.')
parser.add_argument('-u', '--mutation', default=0.0001, type=float, help='Mutation rate in the mutation-selection balance model. Default 1.0E-04.')
parser.add_argument('-A', '--autosomal', default=50, type=int, help='Number of autosomal loci with effect on fitness in the mutation-selection balance model. Default 50.')
parser.add_argument('-X', '--xlinked', default=0, type=int, help='Number of X-linked loci with effect on longevity in the mutation-selection balance model. Default 0.')

parser.add_argument('-N', default=5000, type=int, help='Population size.')
parser.add_argument('-G', default=200, type=int, help='Number of days to evolve before experiment. Default: 200000.')
parser.add_argument('-o', '--output', default='z1.txt', type=argparse.FileType('w'))
parser.add_argument('-r', '--replicates', default=1, type=int, help='Number of replicates. Default: 1.')
parser.add_argument('-S', '--step', default=100, type=int, help='Periodicity of statistics output. Default: 100 (every 100 days).')
args = parser.parse_args()

if re.search("two|phases|smurf", args.demography, re.I):
   args.demography = "two_phases"
   print("Demographic model set to the two-phases model.")

if re.search("gompert?z", args.demography, re.I):
   args.demography = "gompertz"
   print("Demographic model set to Gompertz.")

if re.search("weibull?", args.demography, re.I):
   args.demography = "weibull"
   print("Demographic model set to Weibull.")

if re.search("mut|sel", args.equilibrium, re.I):
   args.equilibrium = "mutation_selection"
   print("Population under mutation-selection balance.")

if re.search("ant|plei", args.equilibrium, re.I):
   args.equilibrium = "antagonistic_pleiotropy"
   print("Population under antagonistic pleiotropy.")


###############################################################
#                    FUNCTIONS AND CLASSES                    #
###############################################################

def aging_model(model):
   def two_phases(smurf, pop):
      k = pop.dvars().k
      if smurf == 1 and random.random() < 1.0 - math.exp(-k):
         return True
      else:
         return False
   def weibull(age, a, b):
      if random.random() < 1.0 - math.exp(-(a/b) * (age + 1) ** b + (a/b) * age ** b):
         return True
      else:
         return False
   def gompertz(age, a, b):
      if random.random() < 1.0 - math.exp( (a * math.exp(b * age) / b) * (1.0 - math.exp(b)) ):
         return True
      else:
         return False
   if model == "two_phases":
      return two_phases
   if model == "gompertz":
      return gompertz
   if model == "weibull":
      return weibull

def fitness_func(geno, ind, pop):
   if ind.sex() == 2:
      if geno[0] + geno[1] == 0:
         value = 1.0 - args.feffect
      if geno[0] + geno[1] == 1:
         value = 1.0 - args.feffect * args.dominance
      if geno[0] + geno[1] == 2:
         value = 1.0
   else:
      value = 1.0
   return value

def MaleEffect(geno, ind):
   if ind.sex() == 2:
      a = min_a
   else:
      # allele "1" makes male aging rate increase.
      a = min_a + args.meffect * geno[0]
   b = args.b
   t0 = -b / a
   return (a, b, t0)

def demo(gen, pop):
   return pop.popSize()

def AdditiveCodominant(geno, ind):
   if ind.totNumLoci() > 0:
      a = min_a + sum(geno) * (max_a - min_a) / ind.totNumLoci()
   else:
      a = args.a
   b = args.b
   return (a, b)

def AdditiveRecessive(geno, ind):
   # males have a shorter 'geno' tuple, because there is only one
   # allele in each locus of the X chromosome. Assuming that X loci
   # are listed before the autosomal ones.
   a = min_a
   gene_effect = (max_a - min_a) / ind.totNumLoci()
   if ind.sex() == 1:
      for X_locus in range(X_loci):
         if geno[X_locus] == 1:
            a += gene_effect
      for A_locus in range(ind.totNumLoci() - X_loci):
         if geno[X_loci + A_locus * 2] + geno[X_loci + A_locus * 2 + 1] == 2:
            a += gene_effect
   elif ind.sex() == 2:
      for locus in range(ind.totNumLoci()):
         if geno[locus * 2] + geno[locus * 2 + 1] == 2:
            a += gene_effect
   b = args.b
   return (a, b)

class sexSpecificRecombinator(sim.PyOperator):
    def __init__(self, intensity=0, rates=0, loci=sim.ALL_AVAIL, convMode=sim.NO_CONVERSION,
            maleIntensity=0, maleRates=0, maleLoci=sim.ALL_AVAIL, maleConvMode=sim.NO_CONVERSION,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.Recombinator = sim.Recombinator(rates, intensity, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = sim.Recombinator(maleRates, maleIntensity,
            maleLoci, maleConvMode)
        sim.PyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)
    def transmitGenotype(self, pop, off, dad, mom):
        # Form the first homologous copy of offspring.
        self.Recombinator.transmitGenotype(mom, off, 0)
        # Form the second homologous copy of offspring.
        self.maleRecombinator.transmitGenotype(dad, off, 1)
        return True

def OutputStats_MutSelBal(pop):
   X = 0
   for locus in range(X_loci):
      X += pop.dvars().alleleFreq[locus][1]
   A = 0
   for locus in range(X_loci, X_loci + A_loci):
      A += pop.dvars().alleleFreq[locus][1]
   if pop.dvars().gen == 0:
      args.output.write("#Gen. \tMeanFreqX \tMeanFreqAu\tMaleAge\tFemAge \tMale_a\tFem_a \n")
   outstring = "{:6d}".format(pop.dvars().gen)
   if X_loci > 0:
      outstring += "\t{:.8f}".format(X/X_loci)
   else:
      outstring += "\t          "
   if A_loci > 0:
      outstring += "\t{:.8f}".format(A/A_loci)
   else:
      outstring += "\t          "
   outstring += "\t{:7.4f}".format(pop.dvars().meanOfInfo_malesAge['age'])
   outstring += "\t{:7.4f}".format(pop.dvars().meanOfInfo_femalesAge['age'])
   outstring += "\t{:.4f}".format(pop.dvars().meanOfInfo_males['a'])
   outstring += "\t{:.4f}\n".format(pop.dvars().meanOfInfo_females['a'])
   args.output.write(outstring)
   return True

def OutputStats_AntaPlei(pop):
   sim.stat(pop, alleleFreq=sim.ALL_AVAIL)
   sim.stat(pop, meanOfInfo='age', subPops=[(0,3)], suffix='_males')
   sim.stat(pop, meanOfInfo='age', subPops=[(0,4)], suffix='_females')
   sim.stat(pop, meanOfInfo='a', subPops=[(0,3)], suffix='_amales')
   sim.stat(pop, meanOfInfo='a', subPops=[(0,4)], suffix='_afemales')
   outstring = str(pop.dvars().gen)
   for locus in range(pop.totNumLoci()):
      outstring += "\t%.3f" % pop.dvars().alleleFreq[locus][1]
   outstring += "\t%3d" % pop.dvars().meanOfInfo_males['age']
   outstring += "\t%3d" % pop.dvars().meanOfInfo_females['age']
   outstring += "\t%4f" % pop.dvars().meanOfInfo_amales['a']
   outstring += "\t%4f\n" % pop.dvars().meanOfInfo_afemales['a']
   args.output.write(outstring)
   return True

