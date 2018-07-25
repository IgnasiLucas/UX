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

parser = argparse.ArgumentParser(description = 'Simulates a cohort and prints the proportion alive every day.')
parser.add_argument('-m', '--model', default='two_phases', type=str, help='Aging model: "two_phases", "gompertz", or "weibull". Default: "two_phases".')
parser.add_argument('-a', default=0.0039, type=float, help='additional fraction of smurfs per unit of time')
parser.add_argument('-b', default=-0.019, type=float, help='intercept of the probability of becoming smurf as a function of age.')
parser.add_argument('-k', default=0.1911, type=float, help='rate of mortality of smurfs.')
parser.add_argument('-N', default=5000, type=int, help='Population size.')
parser.add_argument('-G', default=200, type=int, help='Number of days to simulate. Default: 200.')
parser.add_argument('-o', '--output', default='z1.txt', type=argparse.FileType('w'))
parser.add_argument('-r', '--replicates', default=1, type=int, help='Number of replicates. Default: 1.')
args = parser.parse_args()

###############################################################
#                    FUNCTIONS AND CLASSES                    #
###############################################################

def aging_model(model):
   def two_phases(smurf, pop):
      k = pop.dvars().k
      if smurf == 1.0 and random.random() < 1.0 - math.exp(k):
         return True
      else:
         return False
   def weibull(age, a, b):
      if random.random() < 1.0 - math.exp(-(a/(b+1)) * (age**(b+1) - (age - 1)**(b+1))):
         return True
      else:
         return False
   def gompertz(pop):
      return True
   if model == 'two_phases':
      return two_phases
   if model == 'gompertz':
      return gompertz
   if model == 'weibull':
      return weibull

def demo(gen, pop):
   return pop.popSize()

###############################################################
#                         POPULATION                          #
###############################################################

pop = sim.Population(args.N, loci = 0, ploidy = 2, infoFields = ['age', 'a', 'b', 'smurf', 'luck'])

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
pop.dvars().k = args.k
pop.dvars().model = args.model

###############################################################
#                         SIMULATION                          #
###############################################################

simu = sim.Simulator(pop, rep=args.replicates)

simu.evolve(
   initOps = [
      sim.InitInfo([0], infoFields = 'age'),
      sim.InitInfo([args.a], infoFields = 'a'),
      sim.InitInfo([args.b], infoFields = 'b'),
      sim.InitInfo(lambda: random.random(), infoFields = 'luck'),
      sim.InfoExec("smurf = 1.0 if model == 'two_phases' and ind.luck < ind.age * ind.a + ind.b else 0.0", exposeInd = 'ind'),
      sim.PyExec("Surviving = {'larvae': [], 'adults': [], 'smurfs': []}")
   ],
   preOps = [
      sim.InfoExec("age += 1"),
      sim.InfoExec("luck = random.random()"),
      sim.InfoExec("smurf = 1.0 if model == 'two_phases' and ind.luck <= (ind.age * ind.a + ind.b) else 0.0", exposeInd='ind'),
      sim.DiscardIf(aging_model(args.model))
   ],
   matingScheme = sim.CloneMating(subPops = sim.ALL_AVAIL, subPopSize = demo),
   postOps = [
      sim.Stat(popSize=True, subPops=[(0,0), (0,1), (0,2)]),
      sim.PyExec("Surviving['larvae'].append(subPopSize[0])"),
      sim.PyExec("Surviving['adults'].append(subPopSize[1])"),
      sim.PyExec("Surviving['smurfs'].append(subPopSize[2])"),
#      sim.PyEval(r'"{:d}\t{:d}\t{:d}\t{:d}\n".format(gen, subPopSize[0], subPopSize[1], subPopSize[2])', step=1),
      sim.TerminateIf('popSize == 0')
   ],
   gen=args.G
)

print("#Gen\t" + "\t".join(['Larvae\tAdults\tSmurfs' for a in range(args.replicates)]))
for day in range(args.G):
   line = "{:3d}".format(day)
   for rep in range(args.replicates):
      for kind in ['larvae', 'adults', 'smurfs']:
         try:
            line += "\t{:5d}".format(simu.vars(rep)['Surviving'][kind][day])
         except IndexError:
            line += "\t0"
   print(line)
