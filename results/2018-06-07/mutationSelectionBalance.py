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
parser.add_argument('-m', default=0.0030, type=float, help='Minimum value of the "a" parameter of the two phases model of aging. Default: 0.003.')
parser.add_argument('-M', default=0.0500, type=float, help='Maximum value of the "a" parameter of the two phases model of aging. Default: 0.050.')
parser.add_argument('-k', default=0.1911, type=float, help='rate of mortality of smurfs. Default: 0.1911')
parser.add_argument('-N', default=5000, type=int, help='Population size. Default: 5000.')
parser.add_argument('-G', default=200, type=int, help='Number of generations. Default: 200.')
parser.add_argument('-o', '--output', default='z1.txt', type=argparse.FileType('w'), help='Ouput file name. Default: z1.txt.')
parser.add_argument('-s', '--seed', default=115, type=int, help='Random number generator seed. Default: 115.')
parser.add_argument('-u', '--mutation', default=0.00000001, type=float, help='Mutation rate. Default 1.0E-08.')
args = parser.parse_args()

min_a = args.m
max_a = args.M
X_loci = 7

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

def AdditiveCodominant(geno, ind):
   if ind.totNumLoci() > 0:
      a = min_a + sum(geno) * (max_a - min_a) / ind.totNumLoci()
   else:
      a = args.a
   b = -10.0 * a
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
   b = -a * 10.0
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

def OutputStats(pop):
   outstring = str(pop.dvars().gen)
   for locus in range(pop.totNumLoci()):
      outstring += "\t%.3f" % pop.dvars().alleleFreq[locus][1]
   outstring += "\t%3d" % pop.dvars().maxOfInfo_males['age']
   outstring += "\t%3d" % pop.dvars().maxOfInfo_females['age']
   outstring += "\t%4f" % pop.dvars().meanOfInfo_males['a']
   outstring += "\t%4f\n" % pop.dvars().meanOfInfo_females['a']
   args.output.write(outstring)
   return True

###############################################################
#                         POPULATION                          #
###############################################################

# I simulate 28 loci, distributed along chromosomes X, 2 and 3
# more or less proportionally to their genetic lengths. They are
# all 10 cM apart from the next or previous one in the same chromosome.
pop = sim.Population(args.N, loci = [X_loci,10,11], ploidy = 2,
   chromTypes = [sim.CHROMOSOME_X, sim.AUTOSOME, sim.AUTOSOME],
   infoFields = ['age', 'a', 'b', 'smurf', 'ind_id',  'father_id', 'mother_id', 'luck'])
pop.dvars().seed = args.seed

# Below the population is split first by an age cutoff of 10 days
# (larvae and adults), and then by the aging stage (normal, smurf).
# Both splitters are combined in a product, so that they generate
# VSP 0 (normal larvae), 1 (smurf larvae), 2 (normal adults), and
# 3 (smurf adults). Then I create three additional and independent
# VSPs: males (VSP 4), females (VSP 5), and 0-day larvae (VSP 6).
# Finally, I join 1 and 3, although it does not seem necessary.
pop.setVirtualSplitter(
   sim.CombinedSplitter(
      splitters = [
         sim.ProductSplitter(
            splitters = [
               sim.InfoSplitter(field = 'age', cutoff = 10),
               sim.InfoSplitter(field = 'smurf', values = [0, 1])
            ]
         ),
         sim.SexSplitter(),
         sim.InfoSplitter(field = 'age', values = 0)
      ],
      vspMap = [(0), (2), (1,3), (4), (5), (6)],
      names = ['larvae', 'adults', 'smurfs', 'males', 'females', 'zero']
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
      sim.InitGenotype(freq = [1.0,0.0]),
      sim.InitInfo([0], infoFields = 'age'),
      sim.InitInfo([min_a], infoFields = 'a'),
      sim.InitInfo([-10 * min_a], infoFields = 'b'),
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
               sexSpecificRecombinator(rates=0.01, maleRates=0.0),
               sim.PyQuanTrait(loci = sim.ALL_AVAIL, func = AdditiveRecessive, infoFields = ['a', 'b']) 
            ],
            weight = 1,
            subPops = [(0,1)],
            numOffspring=(sim.UNIFORM_DISTRIBUTION, 10,50)
         )
      ],
      subPopSize = demo
   ),
   postOps = [
      sim.SNPMutator(u=args.mutation, subPops=[(0,5)]),
      sim.Stat(alleleFreq=sim.ALL_AVAIL, step=100),
      sim.Stat(maxOfInfo='age', meanOfInfo='a', subPops=[(0,3)], suffix='_males', step=100),
      sim.Stat(maxOfInfo='age', meanOfInfo='a', subPops=[(0,4)], suffix='_females', step=100),
      sim.PyOperator(func=OutputStats, step=100)
   ],
   gen=args.G
)

pop = simu.extract(0)
sim.dump(pop, max=200)
