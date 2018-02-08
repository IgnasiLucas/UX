# The purpose of this script is to show the relationship between the ratio
# of two functions of the genetic variation affecting a fitness-related trait
# in a population and the average frequency of the deleterious alleles. The
# two functions are: the covariance between average effect of alleles and the
# homozygous dominance effect; and the additive genetic variance.
#
# The number of loci may vary among simulations. All loci are assumed diallelic,
# with one allele causing an increase in the character, and one allele causing
# a decrease. The genotypic values, in increasing order, are '-a', 'd', and 'a',
# where 'a' and 'd' depend on the locus. Most loci may have small effects, and
# only some will have large effects. Thus, 'a' may be taken from an exponential
# distribution. Because a similar total variance should be expected, it seems
# reasonable to make the mean of that exponential distribution inversely proportional
# to the number of loci involved.
#
# Because the trait is assumed to be related to fitness, I will
# assume that low alleles are the deleterious ones. However, no further relationship
# is assumed between the allelic frequency and the size of the effect, 'a'. This
# way, I represent the allele frequency may be affected by factors other than the
# mutation-selection balance.
#
# The low alleles will also be partially recessive. To do so, 'd' will be sampled
# from a uniform distribution between 0 and 'a'.
#
# In a fraction of the simulations there will be a certain portion of loci where the
# low allele's frequency is sampled from a Beta distribution with intermediate mean.
# In those simulations, the portion of loci with intermediate allele frequencies
# will be determined by a uniform distribution.

library(RColorBrewer)
NumSim <- 2000
FracSimWithIntermediateFreq <- 0.7
MinPortionIntermediate <- 0.005
MaxPortionIntermediate <- 0.50
PopulationSize <- 200
MinLoci <- 50
MaxLoci <- 200
FreqClassThreshold <- 0.005

# Parameters of the high-mean Beta distribution. The mean is a / (a + b), where
# 'a' is the first shape parameter of the Beta distribution, and 'b' is the second.
QShape1High <- 7
QShape2High <- 10

# Parameters of the low-mean Beta distribution of allele frequencies.
QShape1Low <- 0.1
QShape2Low <- 20

Ratio <- numeric(length=NumSim)
MeanQ <- numeric(length=NumSim)
PropLociAboveThreshold <- numeric(length=NumSim)
NumLociList <- numeric(length=NumSim)
NumVarLociList <- numeric(length=NumSim)
PortionIntermediateList <- numeric(length=NumSim)

for (i in 1:NumSim) {
   NumLoci <- ceiling(runif(1, min=MinLoci, max=MaxLoci))
   NumLociList[i] <- NumLoci
   NumVarLoci <- NumLoci
   PortionIntermediate <- 0
   if (runif(1, min=0, max=1) <= FracSimWithIntermediateFreq) {
      PortionIntermediate <- runif(1, min=MinPortionIntermediate, max=MaxPortionIntermediate)
   }
   PortionIntermediateList[i] <- PortionIntermediate
   Cad <- 0
   Va <- 0
   Qlist <- numeric(length=NumLoci)
   for (j in 1:NumLoci) {
      a <- rexp(1, rate=1/NumLoci)
      if (runif(1, min=0, max=1) <= PortionIntermediate) {
         q <- rbeta(1, QShape1High, QShape2High)
      } else {
         q <- rbeta(1, QShape1Low, QShape2Low)
      }
      if (PopulationSize > 0) {
         SampledQ <- rbinom(1, 2 * PopulationSize, q) / (2 * PopulationSize)
         q <- SampledQ
      }
      p <- 1 - q
      if (q == 0) {
         NumVarLoci <- NumVarLoci - 1
      }
      Qlist[j] <- q
      d <- runif(1, min=0, max=a)
      Cad <- Cad + 2 * p * q * (p - q) * d * (a + d * (q - p))
      Va <- Va + 2 * p * q * (a + d * (q - p)) ^ 2
   }
   NumVarLociList[i] <- NumVarLoci
   PropLociAboveThreshold[i] <- sum(Qlist >= FreqClassThreshold) / NumLoci
   Ratio[i] <- Cad / Va
   MeanQ[i] <- mean(Qlist[Qlist > 0])
}
pal <- colorRampPalette(c('red', 'blue'))
colourNumVarLoci <- pal(50)[as.numeric(cut(NumVarLociList, breaks=50))]
colourPropQAboveThreshold <- pal(50)[as.numeric(cut(PropLociAboveThreshold, breaks=50))]
colourPortionIntermediate <- pal(50)[as.numeric(cut(PortionIntermediateList, breaks=50))]
png(filename='Kelly1999.png', width=1000, height=1000)
par(mex=2)
plot(c(0.005, 1.0), c(-0.5, 2), type='n', log='x', xlab='Mean q', ylab='Cad/Va', cex.axis=2, cex.lab=2)
points(MeanQ, Ratio, pch='.', cex=5, col=colourPortionIntermediate)
abline(h=0)
dev.off()
rm(a, q, p, i, j, d, Va, Cad, pal, Qlist, SampledQ, NumVarLoci)
