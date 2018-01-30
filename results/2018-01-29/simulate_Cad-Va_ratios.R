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
# way, I allow for the allele frequency to be affected by factors other than the
# mutation-selection bias.
#
# The low alleles will also be partially recessive. To do so, 'd' will be sampled
# from a uniform distribution between 0 and 'a'. 


NumSim <- 100
MinLoci <- 20
MaxLoci <- 400
PropHigh <- seq(from=0.0, to=0.25, by=0.005)
QShape1High <- 7
QShape2High <- 10
QShape1Low <- 0.1
QShape2Low <- 20
plot(c(0.001, 1.0), c(-0.5, 2), type='n', log='x', xlab='Mean q', ylab='Cad/Va')
for (i in 1:NumSim) {
   ratio <- numeric(length=length(PropHigh))
   MedianQ <- numeric(length=length(PropHigh))
   for (k in 1:length(PropHigh)) {
      NumLoci <- ceiling(runif(1, min=MinLoci, max=MaxLoci))
      Cad <- 0
      Va <- 0
      Qlist <- numeric(length=NumLoci)
      for (j in 1:NumLoci) {
         a <- rexp(1, rate=200/NumLoci)
         toss <- runif(1, min=0, max=1)
         if (toss <= PropHigh[k]) {
            q <- rbeta(1, QShape1High, QShape2High)
         } else {
            q <- rbeta(1, QShape1Low, QShape2Low)
         }
         p <- 1 - q
         Qlist[j] <- q
         d <- runif(1, min=0, max=a)
         Cad <- Cad + 2 * p * q * (p - q) * d * (a + d * (q - p))
         Va <- Va + 2 * p * q * (a + d * (q - p)) ^ 2
      }
      ratio[k] <- Cad / Va
      MedianQ[k] <- mean(Qlist)
   }
   points(MedianQ, ratio, pch='.')
}
abline(h=0, col='red')
