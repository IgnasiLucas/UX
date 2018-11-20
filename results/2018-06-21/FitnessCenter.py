###############################################################
#                           MODULES                           #
###############################################################

import simuPOP as sim
import random
import math
import argparse
import re
import numpy

###############################################################
#                 ARGUMENTS AND  VARIABLES                    #
###############################################################

parser = argparse.ArgumentParser(description = 'Simulates a population and reports the average relative fitness of genotypes.')
parser.add_argument('-m', '--model', default='two_phases', type=str, help='Aging model: "two_phases", "gompertz", or "weibull". Default: "two_phases".')
parser.add_argument('-a', default=0.0039, type=float, help='Minimum value of the "a" parameter of the two phases model of aging. Default: 0.0039.')
parser.add_argument('-b', default=-0.0190, type=float, help='Parameter "b" of the two phases model o aging. Default: -0.019.')
parser.add_argument('-k', default=0.1911, type=float, help='Rate of mortality of smurfs. Default: 0.1911')
parser.add_argument('-N', default=50000, type=int, help='Population size. Default: 50000.')
parser.add_argument('-G', default=500, type=int, help='Number of generations. Default: 500.')
parser.add_argument('-e', '--meffect', type=float, default=0.0013, help='Mutant effect on male rate of aging. Default: 0.0013.')
parser.add_argument('-s', '--feffect', type=float, default=0.114, help='Female selective coefficient against wild type allele. Default: 0.114.')
parser.add_argument('-d', '--dominance', type=float, default=0.5, help='Coefficient of dominance of deleterious allele in females. Default: 0.5')
parser.add_argument('-f', '--initfreq', type=float, default=0.5, help='Initial frequency of allele 1. Default: 0.5.')
parser.add_argument('-o', '--output', default='z1', help='Ouput file prefix. Default: z1.')
args = parser.parse_args()
min_a = args.a
if re.search("two|phases|smurf", args.model, re.I):
   args.model = "two_phases"
if re.search("gompert?z", args.model, re.I):
   args.model = "gompertz"
if re.search("weibull?", args.model, re.I):
   args.model = "weibull"

fh = open(args.output + '.ped', 'w')

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

def demo(gen, pop):
   # This prevents running out of reproducers in a cohort of larvae where one
   # or more died.
   if gen < 10:
      return pop.popSize()
   else:
      return args.N

def fitness_func(geno, ind, pop, age):
   # First, I determine age-specific fecundity (between 0 and 1), acording
   # to an arbitrary function. Then, I multiply female's fecundity by (1-s)
   # or (1-hs). This should produce relative fitness values of 1-s and 1-hs
   # if fitness is the expected total number of descendants in a lifetime.
   # At least, among females, because they all share the same survival function.
   # The actual fecundity will just be proportional to these results.
   # Note that I am ignoring the fact that smurfs should not reproduce.
   if age < 10:
      value = 0.0
   else:
      value = 1.0 - ((age - 12) ** 2) / 100
   if value < 0:
      value = 0
   if ind.sex() == 2:
      if geno[0] + geno[1] == 0:
         value *= 1.0 - args.feffect
      if geno[0] + geno[1] == 1:
         value *= 1.0 - args.feffect * args.dominance
      if geno[0] + geno[1] == 2:
         value *= 1.0
   else:
      value *= 1.0
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

###############################################################
#                         POPULATION                          #
###############################################################

pop = sim.Population(args.N, loci = [1], ploidy = 2,
   chromTypes = [sim.CHROMOSOME_X],
   infoFields = ['age', 'a', 'b', 'smurf', 'ind_id',  'father_id', 'mother_id', 'luck', 'fitness', 't0', 'birthday'])
#pop.dvars().seed = args.seed
pop.dvars().min_a = min_a
pop.dvars().meffect = args.meffect
pop.dvars().feffect = args.feffect
pop.dvars().dominance = args.dominance
pop.dvars().k = args.k
pop.dvars().maxAge = 200

# Below the population is split first by an age cutoff of 10 days
# (larvae and adults), and then by the aging stage (normal, smurf).
# Both splitters are combined in a product, so that they generate
# VSP 0 (normal larvae), 1 (smurf larvae), 2 (normal adults), and
# 3 (smurf adults). Then I create two additional and independent
# VSPs: males (VSP 4), and females (VSP 5). Finally, I join 1 and 3.
pop.setVirtualSplitter(
   sim.CombinedSplitter(
      splitters = [
         sim.ProductSplitter(
            splitters = [
               sim.InfoSplitter(field = 'age', cutoff = 10),
               sim.InfoSplitter(field = 'smurf', values = [0, 1])
            ]
         ),
         sim.SexSplitter()
      ],
      vspMap = [(0), (2), (1,3), (4), (5)],
      names = ['larvae', 'adults', 'smurfs', 'males', 'females']
   )
)

# This is to be able to call random and math from InfoExec:
exec("import random", pop.vars(), pop.vars())
exec("import math", pop.vars(), pop.vars())

###############################################################
#                         SIMULATION                          #
###############################################################

simu = sim.Simulator(pop, 1)
simu.evolve(
   initOps = [
      sim.InitSex(),
      sim.InitGenotype(freq = [1.0 - args.initfreq, args.initfreq]),
      sim.InitInfo([0], infoFields = 'age'),
      sim.InitInfo([1], infoFields = 'fitness'),
      sim.InitInfo([0], infoFields = 'birthday'),
      # At this point, even males are diploids! Only 1/1 and 1/0 (but not 0/1) males get 'a' increased.
      sim.InfoExec("a = min_a + meffect if ind.sex() == 1 and ind.allele(0,0) == 1 else min_a", exposeInd = 'ind'),
      sim.InitInfo([ args.b ], infoFields = 'b'),
      sim.InitInfo(lambda: random.random(), infoFields = 'luck'),
      sim.InfoExec("t0 = -ind.b / ind.a", exposeInd = 'ind'),
      sim.InfoExec("smurf = 1.0 if (ind.smurf == 1 or (ind.age > ind.t0 and ind.luck < 1.0 - math.exp(-ind.a * ind.age + ind.a * ind.t0 - ind.a / 2.0))) else 0.0", exposeInd = 'ind'),
      sim.IdTagger(),
      sim.PedigreeTagger(output = '>>{}.ped'.format(args.output), outputFields = ['a', 'birthday'], outputLoci = [0]),
      sim.PyExec("AccumAges = {1: {0: {x: 0 for x in range(maxAge)}, 1: {x: 0 for x in range(maxAge)}, 2: {x: 0 for x in range(maxAge)}}, 2: {0: {x: 0 for x in range (maxAge)}, 1: {x: 0 for x in range(maxAge)}, 2: {x: 0 for x in range(maxAge)}}}")
   ],
   preOps = [
      sim.InfoExec("luck = random.random()"),
      sim.InfoExec("smurf = 1.0 if (ind.smurf == 1 or (ind.age > ind.t0 and ind.luck < 1.0 - math.exp(-ind.a * ind.age + ind.a * ind.t0 - ind.a / 2.0))) else 0.0", exposeInd='ind'),
      sim.DiscardIf(aging_model(args.model)),
      sim.InfoExec("age += 1"),
      # Here, ind.allele(0,1) is 0 for all males, except in first generation.
      sim.InfoExec("AccumAges[ind.sex()][ind.allele(0,0) + ind.allele(0,1)][int(ind.age)] += 1", usePopVars=True, exposeInd='ind'),
      sim.PySelector(loci=[0], func=fitness_func)
   ],
   matingScheme = sim.HeteroMating(
      [
         sim.CloneMating(subPops = [(0,0), (0,1), (0,2)], weight = -1),
         sim.RandomMating(
            ops = [
               sim.IdTagger(),
               sim.InfoExec("smurf = 0.0"),
               sim.InfoExec("birthday = gen"),
               sim.MendelianGenoTransmitter(),
               sim.PyQuanTrait(loci = sim.ALL_AVAIL, func = MaleEffect, infoFields = ['a', 'b', 't0']),
               sim.PedigreeTagger(output = '>>{}.ped'.format(args.output), outputFields = ['a', 'birthday'], outputLoci=[0])
            ],
            weight = 1,
            subPops = [(0,1)],
            numOffspring = 1
         )
      ],
      subPopSize = demo
   ),
   gen = args.G
)

fh.close()

pop = simu.extract(0)
AccumAges = pop.dvars().AccumAges
OffspringNumber = {}
Sex = {}
Genotype = {}
Birthday = {}
with open(args.output + '.ped', 'r') as ped:
   for line in ped:
      Fields = line.split()
      # Initialize the number of children that focal individual will have at each age. Now, all zeros.
      OffspringNumber[int(Fields[0])] = numpy.zeros(50, dtype=numpy.uint64)
      Sex[int(Fields[0])] = Fields[3]
      Genotype[int(Fields[0])] = int(Fields[7]) + int(Fields[8])
      Birthday[int(Fields[0])] = int(Fields[6])
      # The first generation, genotypes are initialized as if all chromosomes were
      # autosomes, and some males are assigned genotype '1 1'. I change that to have
      # consistent hemizygous genotypes in all males.
      if Fields[3] == 'M' and Genotype[int(Fields[0])] == 2:
         Genotype[int(Fields[0])] = 1
      try:
         # This counts the individual as offspring of its _father_ at the father's age.
         OffspringNumber[int(Fields[1])][int(Fields[6]) - Birthday[int(Fields[1])]] += 1
      except KeyError:
         assert Fields[1] == '0'
      try:
         # And now of its mother.
         OffspringNumber[int(Fields[2])][int(Fields[6]) - Birthday[int(Fields[2])]] += 1
      except KeyError:
         assert Fields[2] == '0'

# This will record the numbers of offspring from age-specific parents, across all generations.
# To compare the age-specific fecundities between genotypes, we need to know the overall
# (across generations) numbers of parents at each age. Note the use of a numpy array.
TotalMaleFitness   = {0: numpy.zeros(50, dtype=numpy.uint64), 1: numpy.zeros(50, dtype=numpy.uint64)}
TotalFemaleFitness = {0: numpy.zeros(50, dtype=numpy.uint64), 1: numpy.zeros(50, dtype=numpy.uint64), 2: numpy.zeros(50, dtype=numpy.uint64)}
NumMales = {0: 0, 1: 0}
NumFemales = {0: 0, 1: 0, 2: 0}
AverageMaleFitness = {0: 0, 1: 0}
AverageFemaleFitness = {0: 0, 1: 0, 2: 0}
# This iterates over all individuals simulated, including those who didn't reach adulthood by the end of the simulation...
for ind in OffspringNumber.keys():
   if Sex[ind] == 'M':
      TotalMaleFitness[Genotype[ind]] += OffspringNumber[ind]
      NumMales[Genotype[ind]] += 1
   else:
      TotalFemaleFitness[Genotype[ind]] += OffspringNumber[ind]
      NumFemales[Genotype[ind]] += 1
assert sum(TotalMaleFitness[0] + TotalMaleFitness[1]) == sum(TotalFemaleFitness[0] + TotalFemaleFitness[1] + TotalFemaleFitness[2])
for genotype in range(2):
   try:
      AverageMaleFitness[genotype] = TotalMaleFitness[genotype].sum() / NumMales[genotype]
   except ZeroDivisionError:
      AverageMaleFitness[genotype] = 'nan'
MaxMaleFitness = max([x for x in AverageMaleFitness.values() if x != 'nan'])
# Because the keys (genotypes) are like an index (0, 1...) I can use an array instead of a dictionary.
# Here, I trust that dictionaries preserve the insertion order of both keys and values (new in python 3.7).
RelativeMaleFitness = [ x / MaxMaleFitness if x != 'nan' and MaxMaleFitness > 0 else 'nan' for x in AverageMaleFitness.values()]
for genotype in range(3):
   try:
      AverageFemaleFitness[genotype] = TotalFemaleFitness[genotype].sum() / NumFemales[genotype]
   except ZeroDivisionError:
      AverageFemaleFitness[genotype] = 'nan'
MaxFemaleFitness = max([x for x in AverageFemaleFitness.values() if x != 'nan'])
RelativeFemaleFitness = [ x / MaxFemaleFitness if x != 'nan' and MaxFemaleFitness > 0 else 'nan' for x in AverageFemaleFitness.values()]

with open(args.output + '.fit', 'w') as FitnessFile:
   print("#Model    \tBasic_a\tDelta_a\tFem_s\tDomin.\tFreq.\tSex\tGenot.\tOffs.\tParents\tAve_fit.\tRel_fit.", file=FitnessFile)
   for genotype in range(2):
      print("{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.2f}\tMales\t{:d}\t{:d}\t{:d}\t{}\t{}".format(args.model, args.a, args.meffect,
         args.feffect, args.dominance, args.initfreq, genotype, int(TotalMaleFitness[genotype].sum()), NumMales[genotype],
         AverageMaleFitness[genotype], RelativeMaleFitness[genotype]), file = FitnessFile)
   for genotype in range(3):
      print("{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.2f}\tFemales\t{:d}\t{:d}\t{:d}\t{}\t{}".format(args.model, args.a, args.meffect,
         args.feffect, args.dominance, args.initfreq, genotype, int(TotalFemaleFitness[genotype].sum()), NumFemales[genotype],
         AverageFemaleFitness[genotype], RelativeFemaleFitness[genotype]), file = FitnessFile)
   
with open(args.output + '.age', 'w') as AgeFile:
   # Recall, AccumAges is dictionary created along simulation that records the accumulated number of individuals ever
   # alive right before mating with a specific sex (1 or 2), genotype (0, 1 or 2) and age (up to MaxAge).
   print("#Age\tM0_Abs\tM1_Abs\tF0_Abs\tF1_Abs\tF2_Abs\t\tM0_Num\tM1_Num\tF0_Num\tF1_Num\tF2_Num", file = AgeFile)
   for age in range(1,50):
      print("{: 3}\t{: 6}\t{: 6}\t{: 6}\t{: 6}\t{: 6}\t\t{: 6}\t{: 6}\t{: 6}\t{: 6}\t{: 6}".format(age,
         TotalMaleFitness[0][age], TotalMaleFitness[1][age], TotalFemaleFitness[0][age], TotalFemaleFitness[1][age], TotalFemaleFitness[2][age],
         AccumAges[1][0][age], AccumAges[1][1][age], AccumAges[2][0][age], AccumAges[2][1][age], AccumAges[2][2][age]), file = AgeFile)
