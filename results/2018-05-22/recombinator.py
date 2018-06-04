import simuPOP as sim
from simuPOP.utils import viewVars

# This class is to recombine only in females. Note that it comes from an example
# in the simuPOP manual (sexSpecificRecombinator.py), but Bo Peng had to correct
# it after I told him it didn't work. The default loci and maleLoci were not well
# defined.
class sexSpecificRecombinator(PyOperator):
    def __init__(self, intensity=0, rates=0, loci=ALL_AVAIL, convMode=NO_CONVERSION,
            maleIntensity=0, maleRates=0, maleLoci=ALL_AVAIL, maleConvMode=NO_CONVERSION,
            *args, **kwargs):
        # This operator is used to recombine maternal chromosomes
        self.Recombinator = Recombinator(rates, intensity, loci, convMode)
        # This operator is used to recombine paternal chromosomes
        self.maleRecombinator = Recombinator(maleRates, maleIntensity,
            maleLoci, maleConvMode)
        #
        PyOperator.__init__(self, func=self.transmitGenotype, *args, **kwargs)
    #
    def transmitGenotype(self, pop, off, dad, mom):
        # Form the first homologous copy of offspring.
        self.Recombinator.transmitGenotype(mom, off, 0)
        # Form the second homologous copy of offspring.
        self.maleRecombinator.transmitGenotype(dad, off, 1)
        return True


#pop = Population(2, loci=1000, infoFields=['ind_id', 'father_idx', 'mother_idx'])
#simu = Simulator(pop, rep = 3)
simu = sim.Simulator(
   Population(2, loci=10000, infoFields=['ind_id', 'father_idx', 'mother_idx']),
   rep = 500
)
simu.evolve(
   initOps = [
      InitSex(maleProp = 0.5),
      InitGenotype(prop = [0.5, 0.5]),
      InitLineage()#mode = PER_ALLELE)
   ],
   matingScheme = RandomMating(
      numOffspring = 20,
      subPopSize = 2,
      sexMode = (NUM_OF_MALES, 1),
      ops = [
         sexSpecificRecombinator(rates = 0.0001, maleRates = 0),
         ParentsTagger()
      ]
   ),
   postOps = [
      Stat(inbreeding = ALL_AVAIL, step = 1),
      # This reports a line of 0s and 1s corresponding to the loci in a chromosome, where
      # 0 means the homologous alleles are identical by descent, and 1, the opposite. 
      PyEval(r'"Gen:%d %s %s %.4f\n" % ' +
         r'(gen, "".join([str(int(i - j != 0)) for i,j in zip(LIN[:10000],LIN[10000:20000])]), ' +
         r'"".join([str(int(i - j != 0)) for i,j in zip(LIN[20000:30000],LIN[30000:40000])]), sum(IBD_freq.values()) /len(IBD_freq))',
         stmts="LIN=pop.lineage()",
         exposePop='pop')
   ],
   gen = 20
)
