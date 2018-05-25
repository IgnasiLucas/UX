import simuPOP as sim
import random

random.seed(115)
loci_X = 7
loci_2 = 10
loci_3 = 11
N = 2000
# The chromosomes X, 2, and 3 are roughly 76, 108, and 111 cM long.
# Below, I distribute the number of loci indicated above at regular
# intervals along each chromosome. 
pop = sim.Population(size = N,
   loci = [loci_X, loci_2, loci_3],
   lociPos = list(range(0, 76, round(76.0/loci_X))) + list(range(0, 108, round(108.0/loci_2))) + list(range(0, 111, round(111.0/loci_3))),
   chromNames = ['X', '2', '3'],
   chromTypes = [sim.CHROMOSOME_X, sim.AUTOSOME, sim.AUTOSOME],
   ploidy = 2,
   infoFields = ['age', 't0', 'RIS', 'ind_id', 'father_id', 'mother_id'])

pop.setVirtualSplitter(sim.InfoSplitter(field='age', cutoff=10))


   
