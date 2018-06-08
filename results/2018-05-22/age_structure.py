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

parser = argparse.ArgumentParser(description = 'Simulates a population and prints the age structure')
parser.add_argument('-a', default=0.0039, type=float, help='additional fraction of smurfs per unit of time')
parser.add_argument('-b', default=-0.019, type=float, help='intercept of the probability of becoming smurf as a function of age.')
parser.add_argument('-k', default=0.1911, type=float, help='rate of mortality of smurfs.')
parser.add_argument('-N', default=5000, type=int, help='Population size.')
parser.add_argument('-G', default=200, type=int, help='Number of generations.')
parser.add_argument('-o', '--output', default='z1.txt', type=argparse.FileType('w'))
args = parser.parse_args()

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

def qtrait(geno):
   return (args.a, args.b)

def outputStructure(pop):
   numInd = {1: {0: {}, 1: {}}, 2: {0: {}, 1: {}}}
   sim.stat(pop, popSize=True, maxOfInfo='age')
   maxAge = int(pop.vars()['maxOfInfo']['age'])
   numInd[1][0] = {x: 0 for x in range(maxAge + 1)}
   numInd[1][1] = {x: 0 for x in range(maxAge + 1)}
   numInd[2][0] = {x: 0 for x in range(maxAge + 1)}
   numInd[2][1] = {x: 0 for x in range(maxAge + 1)}
#   with open(args.output, 'w') as out:
   args.output.write("# Generation: %d\n" % pop.dvars().gen)
   args.output.write("# Total size: %d; maximum age: %d\n" % (pop.vars()['popSize'], maxAge))
   args.output.write("# Age\tSmurfMales\tMales\tFemales\tSmurfFemales\n")
   args.output.write("# Parameters of mortality model: a=%.4f, b=%.4f, k=%.4f, t_0=%.4f\n" % (args.a, args.b, args.k, -args.b/args.a))
   for ind in pop.individuals():
      numInd[ind.sex()][int(ind.smurf)][int(ind.age)] += 1
   for age in range(maxAge + 1):
      args.output.write("%d\t%d\t%d\t%d\t%d\n" % (age, numInd[1][1][age], numInd[1][0][age], numInd[2][0][age], numInd[2][1][age]))
   args.output.close()
   return True

###############################################################
#                         POPULATION                          #
###############################################################

pop = sim.Population(args.N, loci = 1, ploidy = 2, infoFields = ['age', 'a', 'b', 'smurf', 'ind_id',  'father_id', 'mother_id', 'luck', 'fitness'])

pop.setVirtualSplitter(
   sim.CombinedSplitter(
      splitters = [
         sim.ProductSplitter(
            splitters = [
               sim.InfoSplitter(field = 'age', cutoff = 10),
               sim.InfoSplitter(field = 'smurf', values = [0, 1])
            ]
         )
      ],
      vspMap = [(0), (2), (1,3)],
      names = ['larvae', 'adults', 'smurfs']
   )
)

# This is to be able to call random from InfoExec:
exec('import random', pop.vars(), pop.vars())

###############################################################
#                         SIMULATION                          #
###############################################################

simu = sim.Simulator(pop, rep=1)

simu.evolve(
   initOps = [
      sim.InitSex(),
      sim.InitGenotype(freq = [0.9,0.1]),
      sim.InitInfo([0], infoFields = 'age'),
      sim.InitInfo([args.a], infoFields = 'a'),
      sim.InitInfo([args.b], infoFields = 'b'),
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
               sim.PyQuanTrait(loci = [0], func = qtrait, infoFields = ['a', 'b']),
               sim.InfoExec("smurf = 0.0"),
               sim.MendelianGenoTransmitter()
            ],
            weight = 1,
            subPops = [(0,1)],
            numOffspring=(sim.UNIFORM_DISTRIBUTION, 10,50)
         )
      ],
      subPopSize = demo
   ),
   gen=args.G
)

pop = simu.extract(0)
outputStructure(pop)
