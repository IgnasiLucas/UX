###############################################################
#                           MODULES                           #
###############################################################

import simuPOP as sim
import random
import math

###############################################################
#                          VARIABLES                          #
###############################################################

random.seed(115)
smurf_death_rate = 0.4

###############################################################
#                    FUNCTIONS AND CLASSES                    #
###############################################################

def natural_death(smurf):
   if smurf == 1.0 and random.random() < smurf_death_rate:
      return True
   else:
      return False

def qtrait(geno):
   return (0.0166667, -0.166667)

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


###############################################################
#                         POPULATION                          #
###############################################################

pop = sim.Population(2000, loci = 1, ploidy = 2, infoFields = ['age', 'a', 'b', 'smurf', 'ind_id',  'father_id', 'mother_id', 'luck', 'fitness'])

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

exec('import random', pop.vars(), pop.vars())

###############################################################
#                         SIMULATION                          #
###############################################################

simu = sim.Simulator(pop, rep=1)

simu.evolve(
   initOps = [
      sim.InitSex(),
      sim.InitGenotype(freq = [0.9,0.1]),
      sim.InitInfo(lambda: random.randint(0,75), infoFields = 'age'),
      sim.InitInfo([0.01667], infoFields = 'a'),
      sim.InitInfo([-0.16667], infoFields = 'b'),
#     sim.InitInfo(lambda : random.betavariate(10.05025, 2000.0), infoFields = 'a'),
#     sim.InitInfo(lambda : random.betavariate(83.33333, 2000.0), infoFields = 'b'),
      sim.InitInfo(lambda: random.random(), infoFields = 'luck'),
#     sim.InitInfo([0.0], infoFields = 'smurf'),
      sim.InfoExec("smurf = 1.0 if ind.luck < ind.age * ind.a + ind.b else 0.0", exposeInd = 'ind'),
      sim.IdTagger()
   ],
   preOps = [
      sim.InfoExec("age += 1"),
      sim.InfoExec("luck = random.random()"),
#      sim.InfoExec("smurf = round(random.random())"),
      sim.InfoExec("smurf = 1.0 if ind.luck < (ind.age * ind.a + ind.b) else 0.0", exposeInd='ind'),
      sim.DiscardIf(natural_death)#,
#      sim.PyEval(r'"gen.%03d\t%d\t%d\t%d\n" % (gen, larvae, adults, smurfs)',
#         stmts="larvae = pop.subPopSize((0,0))\n"
#               "adults = pop.subPopSize((0,1))\n"
#               "smurfs = pop.subPopSize((0,2))\n",
#         exposePop='pop')
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
#              sexSpecificRecombinator(intensity = 0.01, maleIntensity = 0)
               sim.MendelianGenoTransmitter()
            ],
#            polySex = sim.MALE,
#            polyNum = 2,
            weight = 1,
            subPops = [(0,1)],
            numOffspring=(sim.UNIFORM_DISTRIBUTION, 1,5)
         )
      ],
      subPopSize = 2000
   ),
   postOps = [
      sim.Stat(popSize = True, subPops = [(0,0), (0,1), (0,2)]),
      sim.Stat(maxOfInfo='age', subPops = [(0,0)], suffix='_larva'),
      sim.Stat(maxOfInfo='age', subPops = [(0,1)], suffix='_adult'),
      sim.Stat(maxOfInfo='age', subPops = [(0,2)], suffix='_smurf'),
      sim.PyEval(r'"gen.%03d\t%s\t%d\t%d\t%d\t%d\n" % (gen, subPopSize, popSize, maxOfInfo_larva["age"], maxOfInfo_adult["age"], maxOfInfo_smurf["age"])')
   ],
   gen=200
)
