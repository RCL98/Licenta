from generate_primes import logarithm_test_numbers
from random import seed, sample, randint
from Silver_Pohlig_Hellman_algoritm import silver_pohlig_hellman_alg, shoup_alg, pohlig_hellman_alg
from timeit import default_timer as timer
from Shanks_algoritm import Shanks
from Pollard_algoritmi import pollard_rho
from sympy.ntheory import factorint
from functools import partial
from math import ceil, sqrt, log
from aritmetica_modulara import countBits, avg

def test_pohlig_primes(bitlist: list, algQe_type = pohlig_hellman_alg, solver = 'shanks',
                       shanks_type = 'classic', pollard_cycle = 'brent', pollard_func = 'teske', rand_seed = 0):
  seed(rand_seed)
  if solver == 'shanks':
    Solver = partial(Shanks, type = shanks_type, r = 2.25)
  else:
    Solver = partial(pollard_rho, cycle=pollard_cycle, func=pollard_func, rand_seed=rand_seed)
  for ind, bits in enumerate(bitlist):
    p, g, exp, h = logarithm_test_numbers(bits, silvPohHell = True)
    print(list(factorint(p - 1).keys())[-1], countBits(p))
    t = timer()
    x = silver_pohlig_hellman_alg(g, h, p, fact = factorint, qeAlg = algQe_type, solver = Solver)
    if x != exp:
      print("WRONG!")
    print(f"Bits: {bits}, time: {(timer() - t) * 1000: .4f}")

def num_of_tests(bits, numOfTests):
  if bits <= 5:
    tests = 3
  elif bits > 5 and bits <= 10:
    tests = numOfTests * 100
  if bits > 10 and bits < 15:
    tests = numOfTests * 10
  elif bits >= 15 and bits < 25:
    tests = numOfTests * 5
  elif bits >= 25 and bits < 40:
    tests = numOfTests
  elif bits >= 40:
    tests = numOfTests // 100
  return tests

def number_of_distinct_primes(bits):
  nr_primes = {10: 75, 15: 1612, 20: 38635, 25: 985818}
  try:
    primes = nr_primes[bits]
  except:
    primes = ceil((2 ** bits - 1)/log(2 ** bits - 1))
  return primes

def test_qe_solver(listOfBits: list, numOfTests: int, contests: list, rand_seeds = [0]):
  n = len(contests)
  times = [[[0 for _ in rand_seeds] for _ in range(n)] for _ in listOfBits]
  # iterations = [[[0 for _ in rand_seeds] for _ in range(n)] for _ in listOfBits]
  # iterations_avg = [[0 for _ in range(n)] for _ in listOfBits]
  times_avg = [[0 for _ in range(n)] for _ in listOfBits]

  for ind_b, bits in enumerate(listOfBits):
    tests = num_of_tests(bits, numOfTests)
    p_set = set()
    print("\nStarted bits: ", bits)
    print("Number of tests: ", tests)
    print("Seeds:", end=' ')
    for ind_s, s in enumerate(rand_seeds):
      print(s, end=', ')
      seed(s)
      i = 0
      distinct_primes = number_of_distinct_primes(bits)
      last_sum = 0
      while i < tests:
        p, g, exp, h = logarithm_test_numbers(bits, safe=False)
        p_iter = 0

        while (p in p_set or ceil(sqrt(p)) > 30000000) and len(p_set) < distinct_primes and p_iter < 10 ** 2:
          p, g, exp, h = logarithm_test_numbers(bits, safe=False)
          p_iter += 1
        p_set.add(p)

        bit_seed = randint(1, s)
        for ind_c, c in enumerate(contests):

          if c[0] == 'shoup':
            qeAlg = partial(shoup_alg)
          else:
            qeAlg = partial(pohlig_hellman_alg)

          if c[1] == 'shanks':
            Solver = partial(Shanks, type=c[2], r=2.25)
          else:
            Solver = partial(pollard_rho, cycle=c[2], func=c[3], rand_seed=bit_seed)

          t = timer()
          x = silver_pohlig_hellman_alg(g, h, p, fact = factorint, qeAlg = qeAlg, solver = Solver)
          if x != exp:
            print("WRONG " + c[0] + " " + c[1] + ' ' + c[2] + ' ' + c[3])
          times[ind_b][ind_c][ind_s] += (timer() - t) * 1000
          # iterations[ind_b][ind_c][ind_s] += iter
        i += 1

    for ind_c, c in enumerate(contests):
      # iterations[ind_b][ind_c][ind_s] /= tests * ceil(sqrt(p))
      times[ind_b][ind_c][ind_s] /= tests
    print()
    for ind_c, c in enumerate(contests):
      # iterations_avg[ind_b][ind_c] = avg([it for it in iterations[ind_b][ind_c]]) {iterations_avg[ind_b][ind_c]:.4f},
      times_avg[ind_b][ind_c] = avg([it for it in times[ind_b][ind_c]])
      print(
        f"Qe_alg: {c[0]},  solver: {c[1]}: type: {c[2]}, function: {c[3]}, time: {times_avg[ind_b][ind_c]:.3f}")


if __name__ == "__main__":
  numOfTests = 1000
  bitlist = [15, 20, 25, 30, 40, 50, 55, 50]
  seed(42)
  rand_seeds = sample(range(1, 4000), 3)
  contests = [('shoup', 'pollard', 'teske', 'pollard_mod'), ('pohlig-hellman', 'pollard', 'teske', 'pollard_mod')]
  test_qe_solver(bitlist, numOfTests, contests, rand_seeds)
  # bitlist = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
  # import decimal
  # # test_pohlig_primes(bitlist, algQe_type=shoup_alg, shanks_type='general')
  # for bits in bitlist:
  #   p, g, exp, h = logarithm_test_numbers(bits, silvPohHell=True)
  #   print(list(factorint(p - 1).keys())[-1], countBits(p))
  #   p = decimal.Decimal(p)
  #   print(f"P: {p:.4e}")

