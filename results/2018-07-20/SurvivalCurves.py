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

parser = argparse.ArgumentParser(description = 'Simulates a cohort and prints the proportion alive every day.')
parser.add_argument('-m', '--model', default='two_phases', type=str, help='Aging model: "two_phases", "gompertz", or "weibull". Default: "two_phases".')
parser.add_argument('-a', default=0.0039, type=float, help='a parameter. Default: 0.0039.')
parser.add_argument('-b', default=-0.019, type=float, help='b parameter. Default: -0.019')
parser.add_argument('-k', default=0.1911, type=float, help='rate of mortality of smurfs, for two-phases model. Default: 0.1911')
parser.add_argument('-N', default=5000, type=int, help='Population size.')
parser.add_argument('-G', default=200, type=int, help='Number of days to simulate. Default: 200.')
parser.add_argument('-o', '--output', default='z1.txt', type=argparse.FileType('w'))
parser.add_argument('-r', '--replicates', default=1, type=int, help='Number of replicates. Default: 1.')
args = parser.parse_args()

if re.search("two|phases|smurf", args.model, re.I):
   args.model = "two_phases"
if re.search("gompert?z", args.model, re.I):
   args.model = "gompertz"
if re.search("weibull?", args.model, re.I):
   args.model = "weibull"

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
#      if random.random() < 1.0 - math.exp(-(a/b) * (age**b - (age - 1)**b)):
      if random.random() < 1.0 - math.exp(-(a/b) * (age + 1) ** b + (a/b) * age ** b):
         return True
      else:
         return False
   def gompertz(age, a, b):
#      if random.random() < 1.0 - math.exp( -(a * math.exp(b*age) / b) * (1.0 - math.exp(-b))):
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
   return pop.popSize()

###############################################################
#                         POPULATION                          #
###############################################################

pop = sim.Population(args.N, loci = 0, ploidy = 2, infoFields = ['age', 'a', 'b', 'smurf', 'luck', 't0'])

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
exec('import math', pop.vars(), pop.vars())
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
      sim.InfoExec("t0 = -ind.b / ind.a", exposeInd = 'ind'),
      sim.InfoExec("smurf = 1 if (model == 'two_phases' and ind.age > ind.t0 and ind.luck <= 1.0 - math.exp(-ind.a * ind.age + ind.a * ind.t0 - ind.a / 2.0)) else 0", exposeInd = 'ind'),
      sim.PyExec("Surviving = {'larvae': [], 'adults': [], 'smurfs': []}")
   ],
   preOps = [
      sim.InfoExec("luck = random.random()"),
      sim.InfoExec("smurf = 1 if ((ind.smurf == 1) or (model == 'two_phases' and ind.age > ind.t0 and ind.luck <= 1.0 - math.exp(-ind.a * ind.age + ind.a * ind.t0 - ind.a / 2.0))) else 0", exposeInd='ind'),
      sim.DiscardIf(aging_model(args.model)),
      sim.InfoExec("age += 1")
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

args.output.write("#Gen\t" + "\t".join(['Larvae\tAdults\tSmurfs' for a in range(args.replicates)]) + "\n")
for day in range(args.G):
   line = "{:3d}".format(day + 1)
   for rep in range(args.replicates):
      for kind in ['larvae', 'adults', 'smurfs']:
         try:
            line += "\t{:.4f}".format(simu.vars(rep)['Surviving'][kind][day] / args.N)
         except IndexError:
            line += "\t0.0000"
   args.output.write(line + "\n")
