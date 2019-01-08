###############################################################
#                           MODULES                           #
###############################################################

import simuOpt
# I need to set module's options before loading the simuPOP module. And I need to load the
# simuPOP module before defining the SexSpecificRecombinator class, which requires simuPOP.
# Thus, I cannot set the module's options from arguments passed in the command line to the
# main function, unfortunately.
simuOpt.setOptions(numThreads = 20, optimized = False, debug = 'DBG_ALL', alleleType = 'short', quiet = True)
# These are the available debugging codes:
#    'DBG_ALL', 'DBG_GENERAL', 'DBG_UTILITY', 'DBG_POPULATION', 'DBG_OPERATOR', 'DBG_SIMULATOR',
#    'DBG_INDIVIDUAL', 'DBG_MUTATOR', 'DBG_TRANSMITTER', 'DBG_INITIALIZER', 'DBG_STATOR', 'DBG_TAGGER',
#    'DBG_SELECTOR', 'DBG_MATING', 'DBG_MIGRATOR', 'DBG_PROFILE', 'DBG_BATCHTESTING',
#    'DBG_INTEROPERABILITY', 'DBG_COMPATIBILITY', 'DBG_DEVEL', 'DBG_WARNING'
import simuPOP as sim
import random
import math
import argparse
import numpy


###############################################################
#                    FUNCTIONS AND CLASSES                    #
###############################################################

def natural_death(aging_model):
   '''Returns a death function according to the aging model.'''
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
   if aging_model == 'two_phases':
      return two_phases
   if aging_model == 'weibull':
      return weibull
   if aging_model == 'gompertz':
      return gompertz

def fitness_func1(age):
   '''Age-specific probability of being chosen as a parent.'''
   # Age-specific reproductive value is taken from Lin et al. 2014 (Florida Entomologist 94(4):1434-43, Figure 4B).
   # It is defined as 'the contribution of individuals of age x and stage j to the future population'.
   Reproductive_Value = {
      0:   0,  1:  0,  2:  0,  3:  0,  4:  0,  5:  0,  6:  0,  7:  0,  8:  0,  9: 23,
      10: 27, 11: 14, 12: 17, 13: 19, 14: 21, 15: 21, 16: 23, 17: 24, 18: 26, 19: 26,
      20: 26, 21: 27, 22: 25, 23: 23, 24: 23, 25: 22, 26: 20, 27: 17, 28: 16, 29: 13,
      30: 12, 31:  9, 32:  6, 33:  7, 34:  7, 35:  4, 36:  8, 37:  8, 38:  3, 39:  4,
      40:  3, 41:  3, 42:  1, 43:  0, 44:  0, 45:  0, 46:  0, 47:  0, 48:  0, 49:  0,
   }
   if age < 50:
      return Reproductive_Value[age]
   else:
      return 0

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

def demo(pop):
   if pop.subPopSize([0,1]) == 0:
      return pop.popSize()
   else:
      return pop.dvars().N

def TweakAdditiveRecessive(aging_a1, aging_a2, aging_b, X_loci):
   def AdditiveRecessive(geno, ind):
      '''Assigns a and b values of the aging model to newborns, according to genotype.'''
      # males have a shorter 'geno' tuple, because there is only one
      # allele in each locus of the X chromosome. Assuming that X loci
      # are listed before the autosomal ones.
      a = aging_a1
      gene_effect = (aging_a2 - aging_a1) / ind.totNumLoci()
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
      b = aging_b
      return (a, b)
   return AdditiveRecessive

def MutationSelection(N=1000, generations=10000, X_loci=100, A_loci=0, AgingModel='two_phases', seed=2001, reps=1, InitMutFreq=0.001, aging_a1=0.003, aging_a2=0.05, aging_b=-0.019, aging_k=0.1911, MutRate=0.001, StatsStep=100, OutPopPrefix='z1', PrintFreqs=False, debug=False):
   '''Creates and evolves a population to reach mutation-selection balance.'''
   if debug:
      sim.turnOnDebug('DBG_ALL')
   else:
      sim.turnOffDebug('DBG_ALL')
   sim.setRNG('mt19937', seed)
   pop = sim.Population(N, loci = [X_loci, A_loci], ploidy = 2,
      chromTypes = [sim.CHROMOSOME_X, sim.AUTOSOME],
      infoFields = ['age', 'a', 'b', 'smurf', 'ind_id',  'father_id', 'mother_id', 'luck', 't0', 'fitness'])
   pop.setVirtualSplitter(
      sim.CombinedSplitter(
         splitters = [
            sim.ProductSplitter(
               splitters = [
                  sim.InfoSplitter(field = 'age', cutoff = 9),
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
   pop.dvars().k = aging_k
   pop.dvars().N = N
   pop.dvars().seed = seed
   pop.dvars().X_loci = X_loci
   pop.dvars().A_loci = A_loci
   pop.dvars().AgingModel = AgingModel
   exec("import random\nrandom.seed(seed)", pop.vars(), pop.vars())
   exec("import math", pop.vars(), pop.vars())
   simu = sim.Simulator(pop, rep = reps)
   simu.evolve(
      initOps = [
         sim.InitSex(),
         sim.InitGenotype(freq = [1 - InitMutFreq, InitMutFreq]),
         sim.InitInfo([0], infoFields = 'age'),
         sim.InitInfo([aging_a1], infoFields = 'a'),
         sim.InitInfo([aging_b], infoFields = 'b'),
         sim.InitInfo(lambda: random.random(), infoFields = 'luck'),
         sim.InfoExec('t0 = -ind.b / ind.a', exposeInd = 'ind'),
         sim.InfoExec('smurf = 1.0 if AgingModel == "two_phases" and (ind.smurf == 1 or (ind.age > ind.t0 and ind.luck < 1.0 - math.exp(-ind.a * ind.age + ind.a * ind.t0 - ind.a / 2.0))) else 0.0', exposeInd = 'ind'),
         sim.IdTagger(),
         sim.PyExec('XFreqChange={}'),
         sim.PyExec('AFreqChange={}')
      ],
      preOps = [
         sim.InfoExec('luck = random.random()'),
         sim.InfoExec('smurf = 1.0 if AgingModel == "two_phases" and (ind.smurf == 1 or (ind.age > ind.t0 and ind.luck < 1.0 - math.exp(-ind.a * ind.age + ind.a * ind.t0 - ind.a / 2.0))) else 0.0', exposeInd = 'ind'),
         sim.DiscardIf(natural_death(AgingModel)),
         sim.InfoExec('age += 1'),
         sim.PySelector(func = fitness_func1)
      ],
      matingScheme = sim.HeteroMating(
         [
            sim.CloneMating(subPops = [(0,0), (0,1), (0,2)], weight = -1),
            sim.RandomMating(
               ops = [
                  sim.IdTagger(),
                  sim.PedigreeTagger(),
                  sim.InfoExec('smurf = 0.0'),
                  sexSpecificRecombinator(rates = [ 0.75 / X_loci for x in range(X_loci) ] + [2.07 / A_loci for x in range(A_loci) ], maleRates = 0.0),
                  sim.PyQuanTrait(loci = sim.ALL_AVAIL, func = TweakAdditiveRecessive(aging_a1, aging_a2, aging_b, X_loci), infoFields = ['a', 'b'])
               ],
               weight = 1,
               subPops = [(0,1)],
               numOffspring = 1
            )
         ],
         subPopSize = demo
      ),
      postOps = [
         sim.SNPMutator(u=MutRate, subPops=[(0,5)]),
         sim.Stat(alleleFreq=sim.ALL_AVAIL, step=StatsStep),
         sim.IfElse('X_loci > 0',
            ifOps = [sim.PyExec('XFreqChange[gen] = [alleleFreq[x][1] for x in range(X_loci)]')],
            elseOps = [sim.PyExec('XFreqChange[gen] = []')],
            step = StatsStep
         ),
         sim.IfElse('A_loci > 0',
            ifOps = [sim.PyExec('AFreqChange[gen] = [alleleFreq[a][1] for a in range(X_loci, pop.totNumLoci())]', exposePop='pop')],
            elseOps = [sim.PyExec('AFreqChange[gen] = []')],
            step = StatsStep
         ),
         sim.IfElse(PrintFreqs,
            ifOps=[
               sim.PyEval(r"str(rep) + '\t' + str(gen) + '\t' + '\t'.join(map('{0:.4f}'.format, XFreqChange[gen])) + '\t\t' + '\t'.join(map('{0:.4f}'.format, AFreqChange[gen])) + '\n'")
            ],
            step = StatsStep
         ),
         sim.TerminateIf('sum([alleleFreq[x][0] * alleleFreq[x][1] for x in range(X_loci + A_loci)]) == 0')
      ],
      gen = generations
   )
   i = 0
   for pop in simu.populations():
      pop.save('{}_{}.pop'.format(OutPopPrefix, i))
      i += 1

def main():
   parser = argparse.ArgumentParser(description = 'Simulates a population with a specific aging model, where mutations in X-linked and/or autosomal loci affect parameter a of the model. The population is evolved to reach mutation-selection balance.')
   parser.add_argument('-a', '--aging_a1',   default=0.0030, type=float, help='Minimum value of the "a" parameter of the aging model. Default: 0.003, which is appropriate for the two-phases model. Recommended value for Weibul is 0.000485. And for Gompertz, 0.0053.')
   parser.add_argument('-2', '--aging_a2',   default=0.0500, type=float, help='Maximum value of the "a" parameter of the aging model. Default: 0.050, appropriate for the two-phases model.')
   parser.add_argument('-b', '--aging_b',    default=-0.019, type=float, help='Constant value of the "b" parameter of the aging model. Default: -0.019, which is appropriate for the two-phases model. Recommended value for the Weibul model is 2.4746. And for Gomperz, 0.0942.')
   parser.add_argument('-k', '--aging_k',    default=0.1911, type=float, help='Rate of mortality of smurfs in two-phases model of aging. Default: 0.1911.')
   parser.add_argument('-N', '--PopSize',    default=1000,   type=int,   help='Population size. Default: 1000.')
   parser.add_argument('-G', '--generations', default=10000, type=int,   help='Number of simulated generations, which actually represent days. Default: 10000.')
   parser.add_argument('-o', '--OutPopPrefix', default='z1', type=str,   help='Prefix of output file(s). Default: "z1".')
   parser.add_argument('-s', '--seed',       default=2001,   type=int,   help='Random number generator seed. Default: 2001.')
   parser.add_argument('-u', '--MutRate',    default=0.001,  type=float, help='Mutation rate. Default: 0.001.')
   parser.add_argument('-A', '--A_loci',     default=0,      type=int,   help='Number of autosomal loci. Default: 0.')
   parser.add_argument('-X', '--X_loci',     default=100,    type=int,   help='Number of X-linked loci. Default: 100.')
   parser.add_argument('-S', '--StatsStep',  default=100,    type=int,   help='Periodicity of statistics output. Default: 100 generations.')
   parser.add_argument('-q', '--InitMutFreq', default=0.001, type=float, help='Initial mutant frequency. Default: 0.001.')
   parser.add_argument('-m', '--AgingModel', default='two_phases', type=str, choices=['two_phases', 'weibull', 'gompertz'], help='Aging model. Choices: "two_phases", "weibull", "gompertz". Default: "two_phases".')
   parser.add_argument('-r', '--reps',       default=1,      type=int,   help='Number of population replicates to simulate. Default: 1.')
   parser.add_argument('-P', '--PrintFreqs', action='store_true',        help='Print mutant allele frequencies every StatsStep generations. Default: False.')
   parser.add_argument('-D', '--debug',      action='store_true',        help='Turn on debugging. Default: False.')
   args = parser.parse_args()
   MutationSelection(N = args.PopSize,
                     generations = args.generations,
                     X_loci = args.X_loci,
                     A_loci = args.A_loci,
                     AgingModel = args.AgingModel,
                     seed = args.seed,
                     reps = args.reps,
                     InitMutFreq = args.InitMutFreq,
                     aging_a1 = args.aging_a1,
                     aging_a2 = args.aging_a2,
                     aging_b = args.aging_b,
                     aging_k = args.aging_k,
                     MutRate = args.MutRate,
                     StatsStep = args.StatsStep,
                     OutPopPrefix = args.OutPopPrefix,
                     PrintFreqs = args.PrintFreqs,
                     debug = args.debug)

if __name__ == '__main__':
   main()
