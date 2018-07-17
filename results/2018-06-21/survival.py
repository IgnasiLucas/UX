import argparse
import numpy

parser = argparse.ArgumentParser(description='Survival function of the smurf model of aging.')
parser.add_argument('-a', '--min_a', default=0.001, type=float, help='')
parser.add_argument('-A', '--max_a', default=0.050, type=float, help='')
parser.add_argument('-b', '--min_b', default=-0.500, type=float, help='')
parser.add_argument('-B', '--max_b', default=-0.010, type=float, help='')
parser.add_argument('-d', '--death_rate', default=0.1911, type=float, help='')
parser.add_argument('-n', '--num_values', default=50, type=int, help='')
parser.add_argument('-x', '--oldest', default=100, type=int, help='')
args=parser.parse_args()

assert (args.min_a <= args.max_a)
assert (args.min_b <= args.max_b)

def survival(a, b, d, x):
   t0 = -b / a
   t1 = (1.0 - b) / a
   S = 1.0
   if x <= t0:
      return S
   else:
      i = 1
      while i <= min(x, t1):
         S *= 1 - a*i - b
         i += 1
      i = 1
      while i <= min(x, t1):
         P = 1.0
         j = 1
         while j <= i - 1:
            P *= 1 - a * j - b
            j += 1
         S += (a * i + b) * numpy.exp((i-x-1)*d) * P
         i += 1
   return S

Survival = {}
if args.min_a < args.max_a:
   for a in numpy.linspace(args.min_a, args.max_a, num=args.num_values):
      if args.min_b < args.max_b:
         for b in numpy.linspace(args.min_b, args.max_b, num=args.num_values):
            for x in range(1,args.oldest + 1):
               try:
                  Survival[(a,b)][x] = survival(a, b, args.death_rate, x)
               except KeyError:
                  Survival[(a,b)] = {}
                  Survival[(a,b)][x] = survival(a, b, args.death_rate, x)
      else:
         b = args.min_b
         for x in range(1, args.oldest + 1):
            try:
               Survival[(a,b)][x] = survival(a, b, args.death_rate, x)
            except KeyError:
               Survival[(a,b)] = {}
               Survival[(a,b)][x] = survival(a, b, args.death_rate, x)
else:
   a = args.min_a
   if args.b < args.max_b:
      for b in numpy.linspace(args.min_b, args.max_b, num=args.num_values):
         for x in range(1, args.oldest + 1):
            try:
               Survival[(a,b)][x] = survival(a, b, args.death_rate, x)
            except KeyError:
               Survival[(a,b)] = {}
               Survival[(a,b)][x] = survival(a, b, args.death_rate, x)
   else:
      b = args.min_b
      for x in range(1, args.oldest + 1):
         try:
            Survival[(a,b)][x] = survival(a, b, args.death_rate, x)
         except KeyError:
            Survival[(a,b)] = {}
            Survival[(a,b)][x] = survival(a, b, args.death_rate, x)

for x in range(1, args.oldest + 1):
   line = "{}".format(x)
   for ab in sorted(Survival):
      line += "\t{:.6f}".format(Survival[ab][x])
   print(line)
