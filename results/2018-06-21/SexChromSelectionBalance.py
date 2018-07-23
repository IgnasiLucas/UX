###############################################################
#                           MODULES                           #
###############################################################

import simuPOP as sim
import random
import math
import argparse

###############################################################
#                 ARGUMENTS AND  VARIABLES                    #
###############################################################

parser = argparse.ArgumentParser(description = 'Simulates a population and prints allele frequencies out.')
parser.add_argument('-m', default=0.0039, type=float, help='Minimum value of the "a" parameter of the two phases model of aging. Default: 0.003.')
parser.add_argument('-k', default=0.1911, type=float, help='Rate of mortality of smurfs. Default: 0.1911')
parser.add_argument('-b', default=-0.0190, type=float, help='Parameter "b" of the two phases model o aging. Default: -0.019.')
parser.add_argument('-N', default=50000, type=int, help='Population size. Default: 50000.')
parser.add_argument('-G', default=500, type=int, help='Number of generations. Default: 500.')
parser.add_argument('-e', '--meffect', type=float, default=0.0013, help='Mutant effect on male rate of aging. Default: 0.0013.')
parser.add_argument('-s', '--feffect', type=float, default=0.137243, help='Female selective coefficient against wild type allele. Default: 0.137243.')
parser.add_argument('-d', '--dominance', type=float, default=0.5, help='Coefficient of dominance of deleterious allele in females. Default: 0.5')
parser.add_argument('-o', '--output', default='z1.txt', type=argparse.FileType('w'), help='Ouput file name. Default: z1.txt.')
parser.add_argument('-r', '--seed', default=115, type=int, help='Random number generator seed. Default: 115.')
args = parser.parse_args()

min_a = args.m

###############################################################
#                    FUNCTIONS AND CLASSES                    #
###############################################################

def natural_death(smurf):
   if smurf == 1.0 and random.random() < 1.0 - math.exp(-args.k):
      return True
   else:
      return False

def demo(gen, pop):
   if gen < 10:
      return pop.popSize()
   else:
      return args.N

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
   return (a, b)

def OutputStats(pop):
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

###############################################################
#                         POPULATION                          #
###############################################################

pop = sim.Population(args.N, loci = [1], ploidy = 2,
   chromTypes = [sim.CHROMOSOME_X],
   infoFields = ['age', 'a', 'b', 'smurf', 'ind_id',  'father_id', 'mother_id', 'luck', 'fitness'])
pop.dvars().seed = args.seed
pop.dvars().min_a = min_a
pop.dvars().meffect = args.meffect
pop.dvars().feffect = args.feffect
pop.dvars().dominance = args.dominance

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

# This is to be able to call random from InfoExec:
exec("import random\nrandom.seed(seed)", pop.vars(), pop.vars())

###############################################################
#                         SIMULATION                          #
###############################################################

simu = sim.Simulator(pop, rep=1)

simu.evolve(
   initOps = [
      sim.InitSex(),
      sim.InitGenotype(freq = [0.5, 0.5]),
      sim.InitInfo([0], infoFields = 'age'),
      sim.InitInfo([1], infoFields = 'fitness'),
#      sim.InfoExec("fitness = {0: 0.98, 1: 0.99, 2:1.0}[ind.allele(0,0) + ind.allele(0,1)] if ind.sex() == 2 else 0.98", exposeInd = 'ind'),
      sim.InfoExec("a = min_a + meffect if ind.sex() == 1 and ind.allele(0,0) == 1 else min_a", exposeInd = 'ind'),
      sim.InitInfo([ args.b ], infoFields = 'b'),
      sim.InitInfo(lambda: random.random(), infoFields = 'luck'),
      sim.InfoExec("smurf = 1.0 if ind.luck < ind.age * ind.a + ind.b else 0.0", exposeInd = 'ind'),
      sim.IdTagger()
   ],
   preOps = [
      sim.InfoExec("age += 1"),
      sim.InfoExec("luck = random.random()"),
      sim.InfoExec("smurf = 1.0 if ind.luck < (ind.age * ind.a + ind.b) else 0.0", exposeInd='ind'),
      sim.DiscardIf(natural_death)
   ],
   matingScheme = sim.HeteroMating(
      [
         sim.CloneMating(subPops = [(0,0), (0,1), (0,2)], weight = -1),
         sim.RandomMating(
            ops = [
               sim.IdTagger(),
               sim.PedigreeTagger(),
               sim.InfoExec("smurf = 0.0"),
               sim.MendelianGenoTransmitter(),
               sim.PySelector(loci=[0], func=fitness_func),
               sim.PyQuanTrait(loci = sim.ALL_AVAIL, func = MaleEffect, infoFields = ['a', 'b']) 
            ],
            weight = 1,
            subPops = [(0,1)],
            numOffspring=(sim.UNIFORM_DISTRIBUTION, 1,5)
         )
      ],
      subPopSize = demo
   ),
   postOps = [
      sim.PyOperator(func=OutputStats, step=10)
   ],
   gen=args.G
)

pop = simu.extract(0)
sim.dump(pop, max=200)
