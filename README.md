2018-01-29
==========
Simulations that show that the ratio between Cad and Va is actually sensitive
to the allele frequencies. This is a reproduction of one of Kelly (1999) figures
but with a variety of more realistic scenarios.

2018-04-19
==========
Here, I repeat the last simulations with a different parameterization. The
results are equivalent.

2018-05-22
==========
Here I wrote a python script, 'age_structure.py', that implements the two-
phases model of aging in Drosophila. I obtain some plots of the age structure
corresponding to different values of one of the basic parameters of the model,
controlling the rate at which adults become 'smurfs' (committed to die within
the next few days). The script may include some bugs, since I cannot obtain
a population with the values of the parameters estimated by Tricoire and
Rera (2015; PLoS ONE, 10(11): e0141920).

2018-06-07
==========
Starting from the basic python script, here I write one that implements a
simple mutation-selection balance. The idea is to let the 'a' parameter be
increased by deleterious, recessive mutations across the genome, and follow
the allele frequencies up.
