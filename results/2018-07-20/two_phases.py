import argparse
import numpy
import math

parser = argparse.ArgumentParser(description='Survival function of the smurf model of aging.')
parser.add_argument('-a', '--min_a', default=0.001, type=float, help='')
parser.add_argument('-t', '--t0', default=10.0, type=float, help='')
parser.add_argument('-A', '--max_a', type=float, help='')
parser.add_argument('-b', '--min_b', type=float, help='')
parser.add_argument('-B', '--max_b', type=float, help='')
parser.add_argument('-d', '--death_rate', default=0.1911, type=float, help='')
parser.add_argument('-n', '--num_values', default=5, type=int, help='')
parser.add_argument('-x', '--oldest', default=100, type=int, help='')
parser.add_argument('-o', '--output', default='z1.txt', type=argparse.FileType('w'), help='')
args=parser.parse_args()

def survival(a, b, d, x):
   t0 = -b / a
   t1 = (1.0 - b) / a
   S = 1.0
   if x <= t0:
      return S
   else:
      i = math.ceil(t0)
      while i <= min(x, t1):
         S *= 1.0 - a*i - b
         i += 1
      i = math.ceil(t0)
      while i <= min(x, t1):
         P = 1.0
         j = math.ceil(t0)
         while j <= i - 1:
            P *= 1.0 - a * j - b
            j += 1
         S += (a * i + b) * numpy.exp((i-x-1)*d) * P
         i += 1
   return S

if not(args.max_a is None):
   if args.max_a > args.min_a:
      a_values = list(numpy.linspace(args.min_a, args.max_a, num=args.num_values))
   else:
      a_values = [args.min_a for i in range(args.num_values)]
else:
   a_values = [args.min_a for i in range(args.num_values)]

if not(args.min_b is None):
   if not(args.max_b is None):
      if args.max_b > args.min_b:
         b_values = list(numpy.linspace(args.min_b, args.max_b, num=args.num_values))
      else:
         b_values = [args.min_b for i in range(args.num_values)]
   else:
      b_values = [args.min_b for i in range(args.num_values)]
else:
   b_values = []
   for a in a_values:
      b_values.append(-a * args.t0)

Survival = {}

for i in range(args.num_values):
   for x in range(1, args.oldest + 1):
      try:
         Survival[(a_values[i],b_values[i])][x] = survival(a_values[i], b_values[i], args.death_rate, x)
      except KeyError:
         Survival[(a_values[i],b_values[i])] = {}
         Survival[(a_values[i],b_values[i])][x] = survival(a_values[i], b_values[i], args.death_rate, x)

header1 = "#a  "
header2 = "#b  "
header3 = "#t0 "
header4 = "#t1 "
for ab in sorted(Survival):
   header1 += "\t{:8.4f}".format(ab[0])
   header2 += "\t{:8.4f}".format(ab[1])
   header3 += "\t{:8.2f}".format(-ab[1]/ab[0])
   header4 += "\t{:8.2f}".format((1.0 - ab[1]) / ab[0])
header1 += "\n"
header2 += "\n"
header3 += "\n"
header4 += "\n"
args.output.write(header1)
args.output.write(header2)
args.output.write(header3)
args.output.write(header4)
for x in range(1, args.oldest + 1):
   line = "{}".format(x)
   for ab in sorted(Survival):
      line += "\t{:.8f}".format(Survival[ab][x])
   line += "\n"
   args.output.write(line)
