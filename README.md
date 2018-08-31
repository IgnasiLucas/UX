Introduction
============
Motivated by my collaboration with Pau Carazo and Zahida Sultanova, I am here
interested in determining if an experiment of the kind devised by Kelly (1999)
would be able to tell if any X-linked genetic variation that contributes to sex
differences in lifespan is due to either low-frequency deleterious alleles in
(or close to) mutation-selection equilibrium, or intermediate frequency alleles
with sexually antagonistic pleiotropic effects. After a shallow look at the
literature, I find myself unable to figure this out analytically, and I resort
to computer simulations, with the awesome Python package simuPOP.

2018-01-29
==========
Simulations that show that the ratio between Cad and Va is actually sensitive
to the allele frequencies. This is a reproduction of one of Kelly (1999, Genetics
Research 73(3):263-273) figures but with a variety of more realistic scenarios.

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
the next few days). See Tricoire and Rera (2015; PLoS ONE, 10(11): e0141920)
for details of the model.

2018-06-07
==========
I implement the mutation-selection balance in a python script with simuPOP.
The 'a' parameter of the two-phases model of aging is increased additively
by a number of X-linked or autosomal mutations. All mutations have the same
effect, which is smaller the larger the number of loci involved. I reproduce
the unguarded-X effect. However, the age difference between males and females
due to this process is at most of only 2 days. I recall in reality may be 5
or more.

2018-06-21
==========
Before simulating a selection balance with sexually antagonistic pleiotropy,
I need to figure out the fitness effects of mutations that increase the aging
rate. Otherwise, I would be testing values of gene effects until finding a
stable equilibrium. In the two-phases model, it is not easy to translate a
change in aging rate in terms of fitness. I use a numerical method appropriate
for the discrete-time simulations. 

2018-07-20
==========
Realizing that the two-phases of aging is just one among several models of aging,
I decide to implement also the two classic ones: Gompertz and Weibull. Here I
just write a script, SurvivalCurves.py, that simulates a cohort and prints out
the proportion alive every day under any of the three models. I compare the
simulated results with the theoretical survival curve to make sure the models
are well implemented.

2018-08-08
==========
I start writing the script that will simulate Kelly's experiment. Not done
yet.

2018-08-31
==========
Run the sexually antagonistic pleiotropy model to make sure that it produces a
stable equilibrium.
