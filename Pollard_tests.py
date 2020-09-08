from timeit import default_timer as timer
from generate_primes import logarithm_test_numbers
from aritmetica_modulara import sqrt, avg
import random as rand
from math import ceil, log
from Pollard_algoritmi import pollard_rho

def num_of_tests(bits, numOfTests):
  if bits <= 5:
    tests = 3
  elif bits > 5 and bits <= 10:
    tests = numOfTests * 100
  if bits > 10 and bits < 15:
    tests = numOfTests * 10
  elif bits >= 15 and bits < 20:
    tests = numOfTests * 5
  elif bits >= 20 and bits < 25:
    tests = numOfTests * 3
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

def test_pollard(numOfTests: int, listOfBits: list, contests: dict, img_path: str, rand_seeds = [0]):
  import matplotlib.pyplot as plt
  from matplotlib.lines import Line2D
  from itertools import cycle
  marker = cycle(rand.sample(Line2D.markers.keys(), len(contests)))
  plt.tight_layout(pad=0.05)

  times = [[[0 for _ in rand_seeds] for _ in contests] for _ in listOfBits]
  times_avg = [[0 for _ in contests] for _ in listOfBits]
  iterations = [[[0 for _ in rand_seeds] for _ in contests] for _ in listOfBits]
  iterations_avg = [[0 for _ in contests] for _ in listOfBits]

  for ind_b, bits in enumerate(listOfBits):
    print("Bits: ", bits)
    tests = num_of_tests(bits, numOfTests)
    print("Nr of tests: ", tests)
    p_set = set()
    print("Seed: ", end = '')
    for ind_s, seed in enumerate(rand_seeds):
      rand.seed(seed)
      print(seed, end=', ')
      distinct_primes = number_of_distinct_primes(bits)
      for _ in range(tests):
        p, g, exp, h = logarithm_test_numbers(bits, safe=False)
        p_iter = 0
        while p in p_set and len(p_set) < distinct_primes - 1 and p_iter < 10**2:
          p, g, exp, h = logarithm_test_numbers(bits, safe=False)
          p_iter += 1
        p_set.add(p)

        bit_seed = rand.randint(1, seed)
        for ind_c, c in enumerate(contests):
          t = timer()
          x, iter = pollard_rho(g, h, p, cycle= c[0], func = c[1], rand_seed=bit_seed)
          if x != exp:
            print("WRONG " + c[0] + " " + c[1])
          times[ind_b][ind_c][ind_s] += (timer() - t) * 1000
          iterations[ind_b][ind_c][ind_s] += iter

      for ind_c, c in enumerate(contests):
        iterations[ind_b][ind_c][ind_s] /= tests
        times[ind_b][ind_c][ind_s] /= tests
    print()
    for ind_c, c in enumerate(contests):
      iterations_avg[ind_b][ind_c] = avg(iterations[ind_b][ind_c])/ceil(sqrt(p))
      times_avg[ind_b][ind_c] = avg(times[ind_b][ind_c])
      print(f"Cycle: {c[0]}, function: {c[1]}, itertions: {iterations_avg[ind_b][ind_c]:.4f}, time: {times_avg[ind_b][ind_c]:.3f}")
    print()

  plt.figure(0)
  for ind_c, c in enumerate(contests):
    the_times = [t[ind_c] for t in times_avg]
    plt.plot(listOfBits, the_times, marker=next(marker), label=c[0] + '_' + c[1])
  plt.legend()
  plt.ylabel("Milisecunde")
  plt.xlabel("biti")
  # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  plt.yscale('log')
  plt.savefig(img_path + '/pollard_time.png', bbox_inches='tight')

  plt.figure(1)
  print("The avg number of iterations: ", end=''); avgits = []
  for ind_c, c in enumerate(contests):
    the_iterations = [it[ind_c] for it in iterations_avg]
    avg_c = avg(the_iterations)
    print(c[0] + '-' + c[1] + f": {avg_c:.4f}", end=' ')
    avgits.append(avg_c)
    plt.plot(listOfBits, the_iterations, marker=next(marker), label= c[0] + '_' + c[1])
  plt.legend()
  plt.ylabel("Iteratii")
  plt.xlabel("biti")
  # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  plt.yscale('log')
  plt.savefig(img_path + '/pollard_iter.png', bbox_inches='tight')

  # plt.figure(2)
  # print("The avg grouop operations: ", end=''); avgos = []
  # for ind_c, c in enumerate(contests):
  #   the_gos = [go[ind_c] for go in group_operations]
  #   avg_c = avg(the_gos)/ceil(sqrt(p))
  #   print(c[0] + '-' + c[1] + f": {avg_c:.4f}", end=' ')
  #   avgos.append(avg_c)
  #   plt.plot(listOfBits, the_gos, marker=next(marker), label= c[0] + '_' + c[1])
  # plt.legend()
  # plt.ylabel("Operații de grup")
  # plt.xlabel("biti")
  # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  # plt.yscale('log')
  # plt.savefig(img_path + '/pollard_iter.png', bbox_inches='tight')
  plt.show()
  return avgits

def test_pollard_dellays(numOfTests: int, listOfBits: list, contests: dict, img_path: str, rand_seeds = [0]):
  iterations = [[[0 for _ in rand_seeds] for _ in contests] for _ in listOfBits]
  iterations_avg = [[0 for _ in contests] for _ in listOfBits]

  lambs = [[[0 for _ in rand_seeds] for _ in contests] for _ in listOfBits]
  lambs_avg = [[0 for _ in contests] for _ in listOfBits]

  mus = [[[0 for _ in rand_seeds] for _ in contests] for _ in listOfBits]
  mus_avg = [[0 for _ in contests] for _ in listOfBits]

  for ind_b, bits in enumerate(listOfBits):
    print("Bits: ", bits)
    tests = num_of_tests(bits, numOfTests)
    print("Nr of tests: ", tests)
    p_set = set()
    print("Seed: ", end = '')
    for ind_s, seed in enumerate(rand_seeds):
      rand.seed(seed)
      print(seed, end=', ')
      distinct_primes = number_of_distinct_primes(bits)
      for _ in range(tests):
        p, g, exp, h = logarithm_test_numbers(bits, safe=False)
        p_iter = 0
        while p in p_set and len(p_set) < distinct_primes - 1 and p_iter < 10**2:
          p, g, exp, h = logarithm_test_numbers(bits, safe=False)
          p_iter += 1
        p_set.add(p)

        bit_seed = rand.randint(1, seed)
        for ind_c, c in enumerate(contests):
          x, iter, lamb, mu = pollard_rho(g, h, p, cycle= c[0], func = c[1], rand_seed=bit_seed)
          if x != exp:
            print("WRONG " + c[0] + " " + c[1])
          iterations[ind_b][ind_c][ind_s] += iter
          lambs[ind_b][ind_c][ind_s] += lamb
          mus[ind_b][ind_c][ind_s] += mu

      for ind_c, c in enumerate(contests):
        iterations[ind_b][ind_c][ind_s] /= tests
        lambs[ind_b][ind_c][ind_s] /= tests
        mus[ind_b][ind_c][ind_s] /= tests
    print()
    for ind_c, c in enumerate(contests):
      iterations_avg[ind_b][ind_c] = avg(iterations[ind_b][ind_c])/ceil(sqrt(p))
      lambs_avg[ind_b][ind_c] = avg(lambs[ind_b][ind_c]) / ceil(sqrt(p))
      mus_avg[ind_b][ind_c] = avg(mus[ind_b][ind_c]) / ceil(sqrt(p))
      print(f"Cycle: {c[0]}, function: {c[1]}, itertions: {iterations_avg[ind_b][ind_c]:.4f}, lambs: {lambs_avg[ind_b][ind_c]:.4f} mus: {mus_avg[ind_b][ind_c]:.4f}, "
            f"rho: {lambs_avg[ind_b][ind_c] + mus_avg[ind_b][ind_c]:.4f} delay: {iterations_avg[ind_b][ind_c]/(lambs_avg[ind_b][ind_c] + mus_avg[ind_b][ind_c]):.4f}")
    print()

  for ind_c, c in enumerate(contests):
    the_iterations = avg([it[ind_c] for it in iterations_avg])
    the_lambs = [lamb[ind_c] for lamb in lambs_avg]
    the_mus = [mu[ind_c] for mu in mus_avg]
    the_rhos = avg([x + y for x, y in zip(the_lambs, the_mus)])
    print(
      f"Cycle: {c[0]}, function: {c[1]}, itertions_avg: {the_iterations:.4f}, rho: {the_rhos:.4f}, delays {the_iterations/the_rhos:.4f}")

if __name__ == "__main__":
  numOfTests = 100
  bits = [15, 20, 25, 30]
  rand.seed(42)
  rand_seeds = rand.sample(range(0, 5000), 10)
  contests = [('teske_spec', 'pollard'), ('teske_spec', 'pollard_mod'), ('teske_spec', 'teske'), ('teske_spec', 'combined')]
  # contests = [('brent', 'pollard'), ('brent', 'pollard_mod'), ('brent', 'teske'), ('brent', 'combined')]
  img_path = 'imagini/'
  iter = test_pollard_dellays(numOfTests, bits, contests, img_path, rand_seeds)

