#!/bin/bash
#
#				2018-07-20
#				----------
#
# Till now, I have used the two-phases model of aging, which is not a
# necessary part of the unguarded X hypothesis. I should implement other
# models, such as Weibull and Gompertz. Plus, I am not satisfied with the
# test run on 2018-05-22, because Tricoire and Rera's (2015) estimates of
# the parameters of the model produce older populations than observed in
# our experience (Pau Carazo).
#
# It turns out that there is a whole survival analysis theory that I was
# not aware about. It is applied both to demography and to the analysis
# of the lifetime of industrial products.
#
# There are four basic functions:
#
# 1. The probability density function of the time of the event (failure or
#    death), f(t), shows how the events are distributed along the lifetime.
# 2. The cumulative probability function, F(t) is the probability of an
#    event occurring at or before time t.
# 3. The survival or reliability function, S(t) is the complement of F(t),
#    and represents the probability of surviving at least until time t.
# 4. The hazard function is not a probability, but the instantaneous rate
#    of event occurrence per unit of time. h(t) = f(t)/S(t). Note that it
#    is conditional to having survived up to time t, and therefore it is
#    independent of the proportion of survivors. It informs about how components
#    fail: early in life due to manufacturing defects, randomly along the
#    lifetime, or increasingly with time because the product wears out.
#
# Having one function defined, the other ones can be found. Since models
# are usualy defined in these terms, the question is how to implement a
# specific survival model. Simulations are individual-based. Every cycle,
# individuals get one day older, and then they either die or not, before
# mating at that age. The probability with which they die at cycle t should
# represent the probability of dying at any point between age t-1 and age t.
# In principle, this can be obtained from the survival function as
# S(t-1) - S(t). However, note that this is the probability of dying at that
# interval before knowing if it survived that far or not; at birth, so to
# speak. That is, S(t-1)-S(t) represents the fraction of a cohort that will
# die between ages t-1 and t. When running a simulation, at each step we are
# only concerned with the individuals that survived so far.
#
# Ricklefs and Scheuerlein (2002, The Journals of Gerontology: Series A,
# 57(2):B69–B76) simulated some datasets under Gompertz and Weibull models:
#
#  "For each simulation, the time scale was divided into 0.1-year intervals,
#   and the probability of death during the interval was calculated from the
#   model as the difference in the expected survival to the beginning and end
#   of each time interval. Then, for each individual in the simulation a random
#   number from a uniform distribution between 0 and 1 was generated for each
#   age interval beginning at age 0. The age at death for each individual was
#   the midpoint of the first age interval in which its random number was less
#   than the probability of mortality during that interval."
#
# Apparently, they used S(t-1) - S(t) as the probability of death at each
# interval. However, after drawing as many random numbers as intervals, and
# comparing them with the interval's probability  of death, it is clear that
# any age of death is not picked with its alleged probability, because not only
# the random number has to be below the probability, but none of the previous
# intervals must have been picked first.
#
# The right probability of death between ages t-1 and t experienced by an
# individual who has survived up to age t-1 is:
#
#                                    S(t-1) - S(t)
#   P(t - 1 < T <= t | T > t - 1) = ---------------
#                                       S(t-1)
#
# An alternative implementation would be to assign a death time to every
# newborn, according to the f(t). Note that when S(t) = exp(-kt), as in the
# case of the smurfs in the two-phases model, the above equation results
# in 1-exp(-k), just as I was using in previous simulations. 
#
#
# Weibull
# =======
#
# In the Weibull model, S(t) = exp(-(kt)^p). Thus, the probability of death
# between ages t-1 and t, having survived up tot age t-1 must be:
#
#   1 - exp{-k^p · (t^p - (t-1)^p}
#
# Where k and p are the two parameters of the Weibull distribution. However,
# Tricoire and Rera (2015) report estimated parameters a=0.000485 and b=2.4746
# for the Weibull model, without specifying the parameterization. It is likely
# that they used the definition of the hazard function of the Weibull model
# as this:
#
#    h(t) = a·t^b
#
# In terms of k and p, as in S(t)=exp(-(kt)^p), the hazard function is:
#
#    h(t) = p·k·(k·t)^(p-1)
#
# I believe that the estimates by Tricoire and Rera are equivalent to
# p = k+1 = 3.4746, and k = exp((log(a)-log(b+1))/(b+1)) = -2.5548
