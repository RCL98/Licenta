import numpy as np
from scipy import integrate, inf
from math import exp, log2, log10, sqrt, ceil
from random import sample, seed

def check_ineq(n, Ns, f):
  diff = np.empty([len(Ns)])
  vectf = np.vectorize(f)
  inter = n * integrate.quad(f, 0, inf, limit = 1000, epsabs=0, epsrel=1.49e-14)[0]
  div = np.vectorize(lambda x: x / n)
  for ind, N in enumerate(Ns):
    xs = np.linspace(0, N, num = N)
    xs = div(xs)
    # print(np.sum(vectf(xs)))
    diff[ind] = np.sum(vectf(xs)) - inter
  return diff

def hash(val: int, k: int):
  return (val & ((1 << k) - 1))

def simulate_points(N: int, n: int, c: int):
  D = []
  k = -ceil(log2(n/sqrt(N)))
  for x in range(1, N):
    if hash(x, k) == 0:
      D.append(x)
  print(f"B: {len(D)/ceil(sqrt(N)):.4f}, logN: {len(D)/(c * ceil(sqrt(N))):.4f}")

if __name__ == "__main__":
  # Ns = [10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
  # cs = [1 + 1/10, 1 + 1/8, 1 + 1/6, 1 + 1/4, 1 + 1/2, 1 + 1]
  # for c in cs:
  #   print(c)
  #   for N in Ns:
  #     simulate_points(N, c, 0.75)
  #   print()
  # f = lambda x: (x**2)*exp(-(x**2) / 2)
  # Ns = [10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
  # diffs = check_ineq(10, Ns, f)
  # print(diffs)
  for i in range(5):
    if i % 2:
      i = i - 1
      continue
