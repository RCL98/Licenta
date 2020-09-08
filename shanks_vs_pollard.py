from math import ceil, sqrt, log
import Pollard_algoritmi as poll
import Shanks_algoritm as shanks
from generate_primes import logarithm_test_numbers
from random import seed, randint, sample
from aritmetica_modulara import avg, countBits
from timeit import default_timer as timer

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


def test_shanks_vs_pollard(numOfTests: int, listOfBits: list, contests: list, rand_seeds=[0]):
  n = len(contests) + 2
  times = [[[0 for _ in rand_seeds] for _ in range(n)] for _ in listOfBits]
  iterations = [[[0 for _ in rand_seeds] for _ in range(n)] for _ in listOfBits]
  iterations_avg = [[0 for _ in range(n)] for _ in listOfBits]
  times_avg = [[0 for _ in range(n)] for _ in listOfBits]

  for ind_b, bits in enumerate(listOfBits):
    tests = num_of_tests(bits, numOfTests)
    p_set = set()
    print("\nStarted bits: ", bits)
    print("Number of tests: ", tests)
    print("Seeds:", end=' ')
    for ind_s, s in enumerate(rand_seeds):
      print(s, end=',\n')
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
        print(f"sqrt(p): {ceil(sqrt(p))}")

        t = timer()
        x, iter = poll.pollard_simple(g, h, p, rand_seed=bit_seed)
        if x != exp:
          if x == None:
            continue
          print("WRONG Pollard simplu")
        times[ind_b][n - 2][ind_s] += (timer() - t) * 1000
        iterations[ind_b][n - 2][ind_s] += iter
        # print("Poll simplu Done!:", (timer() - t) * 1000)

        t = timer()
        x, sm, gs, _ = shanks.shanks_classic(g, h, p, p - 1)
        if x != exp:
          if x == None:
            continue
          print("WRONG Shanks clasic")
        times[ind_b][n - 1][ind_s] += (timer() - t) * 1000
        iterations[ind_b][n - 1][ind_s] += sm + gs
        # print("Shanks Done!", (timer() - t) * 1000)

        for ind_c, c in enumerate(contests):
          t = timer()
          x, iter = poll.pollard_rho(g, h, p, cycle=c[0], func=c[1], rand_seed=bit_seed)
          if x != exp:
            print("WRONG " + c[0] + " " + c[1])
          times[ind_b][ind_c][ind_s] += (timer() - t) * 1000
          iterations[ind_b][ind_c][ind_s] += iter
          new_sum = sum([t[0] for t in times[ind_b][:]])
        print(f"Test {i} done, time: {new_sum - last_sum:.4f}")
        last_sum = new_sum
        i += 1

        # t = timer()
        # x, op, _ = shanks.shanks_two_grumpys_one_baby(g, h, p, p - 1)
        # if x != exp:
        #   print("WRONG Shanks clasic")
        # times[ind_b][n - 1][ind_s] += (timer() - t) * 1000
        # iterations[ind_b][n - 1][ind_s] += op

      for ind_c, c in enumerate(contests):
        iterations[ind_b][ind_c][ind_s] /= tests * ceil(sqrt(p))
        times[ind_b][ind_c][ind_s] /= tests
      iterations[ind_b][n - 2][ind_s] /= tests * ceil(sqrt(p));
      times[ind_b][n - 2][ind_s] /= tests
      iterations[ind_b][n - 1][ind_s] /= tests * ceil(sqrt(p));
      times[ind_b][n - 1][ind_s] /= tests
      # iterations[ind_b][n - 1][ind_s] /= tests * ceil(sqrt(p)); times[ind_b][n - 1][ind_s] /= tests
    print()
    for ind_c, c in enumerate(contests):
      iterations_avg[ind_b][ind_c] = avg([it for it in iterations[ind_b][ind_c]])
      times_avg[ind_b][ind_c] = avg([it for it in times[ind_b][ind_c]])
      print(
        f"Cycle: {c[0]},  function: {c[1]}: itertions: {iterations_avg[ind_b][ind_c]:.4f}, time: {times_avg[ind_b][ind_c]:.3f}")
    iterations_avg[ind_b][n - 2] = avg([it for it in iterations[ind_b][n - 2]])
    times_avg[ind_b][n - 2] = avg([it for it in times[ind_b][n - 2]])
    print(f"Pollard simple: itertions: {iterations_avg[ind_b][n - 2] :.4f}, time: {times_avg[ind_b][n - 2]:.3f}")
    iterations_avg[ind_b][n - 1] = avg([it for it in iterations[ind_b][n - 1]])
    times_avg[ind_b][n - 1] = avg([it for it in times[ind_b][n - 1]])
    print(f"Shanks clasic: itertions: {iterations_avg[ind_b][n - 1] :.4f}, time: {times_avg[ind_b][n - 1]:.3f}")
    # iterations_avg[ind_b][n - 1] = avg([it for it in iterations[ind_b][n - 1]])
    # times_avg[ind_b][n - 1] = avg([it for it in times[ind_b][n - 1]])
    # print(f"Two giants: itertions: {iterations_avg[ind_b][n - 1]:.4f}, time: {times_avg[ind_b][n - 1]:.3f}")

    print("Finished bits: ", bits)

  print("\nThe avg number of iterations: ", end='');
  avgits = []
  for ind_c, c in enumerate(contests):
    avg_c = avg([it[ind_c] for it in iterations_avg])
    print(c[0] + '-' + c[1] + f": {avg_c:.4f}", end=' ')
    avgits.append(avg_c)
  avg_c = avg([it[n - 2] for it in iterations_avg])
  print(f"pollard simple: {avg_c:.4f}", end=', ')
  avg_c = avg([it[n - 1] for it in iterations_avg])
  print(f"shanks clasic: {avg_c:.4f}", end=', ')
  # avg_c = avg([it[n - 1] for it in iterations_avg])
  # print(f"two giants: {avg_c:.4f}", end=' .')

  return iterations_avg

# def test_shanks_vs_pollard(numOfTests: int, listOfBits: list, contests: dict, rand_seeds = [0]):
#   n = len(contests) + 2
#   times = [[[0 for _ in rand_seeds] for _ in range(n)] for _ in listOfBits]
#   iterations = [[[0 for _ in rand_seeds] for _ in range(n)] for _ in listOfBits]
#   iterations_avg = [[0 for _ in range(n)] for _ in listOfBits]
#   times_avg = [[0 for _ in range(n)] for _ in listOfBits]
#   for ind_b, bits in enumerate(listOfBits):
#     tests = num_of_tests(bits, numOfTests)
#     p_set = set()
#     print("\nStarted bits: ", bits)
#     print("Seeds:", end = ' ')
#     for ind_s, s in enumerate(rand_seeds):
#       print(s, end=', ')
#       seed(s)
#       i = 0
#       distinct_primes = number_of_distinct_primes(bits)
#       while i < tests:
#         p, g, exp, h = logarithm_test_numbers(bits, safe=False)
#         p_iter = 0
#         while (p in p_set or ceil(sqrt(p)) > 30000000) and len(p_set) < distinct_primes and p_iter < 10**2:
#           p, g, exp, h = logarithm_test_numbers(bits, safe=False)
#           p_iter += 1
#         p_set.add(p)
#
#         bit_seed = randint(1, s)
#         # print(f"sqrt(p): {ceil(sqrt(p))}")
#
#         t = timer()
#         x, iter = poll.pollard_simple(g, h, p, rand_seed=bit_seed)
#         if x != exp:
#           if x == None:
#               continue
#           print("WRONG Pollard simplu")
#         times[ind_b][n - 2][ind_s] += (timer() - t) * 1000
#         iterations[ind_b][n - 2][ind_s] += iter
#         # print("Poll simplu Done!:", (timer() - t) * 1000)
#
#         t = timer()
#         x, sm, gs, _ = shanks.shanks_classic(g, h, p, p - 1)
#         if x != exp:
#           if x == None:
#               continue
#           print("WRONG Shanks clasic")
#         times[ind_b][n - 1][ind_s] += (timer() - t) * 1000
#         iterations[ind_b][n - 1][ind_s] += sm + gs
#         # print("Shanks Done!", (timer() - t) * 1000)
#
#
#         for ind_c, c in enumerate(contests):
#           t = timer()
#           x, iter = poll.pollard_rho(g, h, p, cycle=c[0], func=c[1], rand_seed=bit_seed)
#           if x != exp:
#             print("WRONG " + c[0] + " " + c[1])
#           times[ind_b][ind_c][ind_s] += (timer() - t) * 1000
#           iterations[ind_b][ind_c][ind_s] += iter
#         # print("Polls done")
#         i += 1
#
#         # t = timer()
#         # x, op, _ = shanks.shanks_two_grumpys_one_baby(g, h, p, p - 1)
#         # if x != exp:
#         #   print("WRONG Shanks clasic")
#         # times[ind_b][n - 1][ind_s] += (timer() - t) * 1000
#         # iterations[ind_b][n - 1][ind_s] += op
#
#       for ind_c, c in enumerate(contests):
#         iterations[ind_b][ind_c][ind_s] /= tests * ceil(sqrt(p))
#         times[ind_b][ind_c][ind_s] /= tests
#       iterations[ind_b][n - 2][ind_s] /= tests * ceil(sqrt(p)); times[ind_b][n - 2][ind_s] /= tests
#       iterations[ind_b][n - 1][ind_s] /= tests * ceil(sqrt(p)); times[ind_b][n - 1][ind_s] /= tests
#       # iterations[ind_b][n - 1][ind_s] /= tests * ceil(sqrt(p)); times[ind_b][n - 1][ind_s] /= tests
#     print()
#     for ind_c, c in enumerate(contests):
#       iterations_avg[ind_b][ind_c] = avg([it for it in iterations[ind_b][ind_c]])
#       times_avg[ind_b][ind_c] = avg([it for it in times[ind_b][ind_c]])
#       print(f"Cycle: {c[0]},  function: {c[1]}: itertions: {iterations_avg[ind_b][ind_c]:.4f}, time: {times_avg[ind_b][ind_c]:.3f}")
#     iterations_avg[ind_b][n - 2] = avg([it for it in iterations[ind_b][n - 2]])
#     times_avg[ind_b][n - 2] = avg([it for it in times[ind_b][n - 2]])
#     print(f"Pollard simple: itertions: {iterations_avg[ind_b][n - 2] :.4f}, time: {times_avg[ind_b][n - 2]:.3f}")
#     iterations_avg[ind_b][n - 1] = avg([it for it in iterations[ind_b][n - 1]])
#     times_avg[ind_b][n - 1] = avg([it for it in times[ind_b][n - 1]])
#     print(f"Shanks clasic: itertions: {iterations_avg[ind_b][n - 1] :.4f}, time: {times_avg[ind_b][n - 1]:.3f}")
#     # iterations_avg[ind_b][n - 1] = avg([it for it in iterations[ind_b][n - 1]])
#     # times_avg[ind_b][n - 1] = avg([it for it in times[ind_b][n - 1]])
#     # print(f"Two giants: itertions: {iterations_avg[ind_b][n - 1]:.4f}, time: {times_avg[ind_b][n - 1]:.3f}")
#
#     print("Finished bits: ", bits)
#
#   print("\nThe avg number of iterations: ", end='');
#   avgits = []
#   for ind_c, c in enumerate(contests):
#     avg_c = avg([it[ind_c] for it in iterations_avg])
#     print(c[0] + '-' + c[1] + f": {avg_c:.4f}", end=' ')
#     avgits.append(avg_c)
#   avg_c = avg([it[n - 2] for it in iterations_avg])
#   print(f"pollard simple: {avg_c:.4f}", end=', ')
#   avg_c = avg([it[n - 1] for it in iterations_avg])
#   print(f"shanks clasic: {avg_c:.4f}", end=', ')
#   # avg_c = avg([it[n - 1] for it in iterations_avg])
#   # print(f"two giants: {avg_c:.4f}", end=' .')
#
#   return iterations_avg

# if __name__ == "__main__":
#   numOfTests = 100
#   bits = [15]
#   seed(42)
#   rand_seeds = sample(range(0, 5000), 1)
#   contests = [('brent', 'pollard'), ('floyd', 'pollard'), ('teske', 'pollard')]
#   iters = test_shanks_vs_pollard(numOfTests, bits, contests, rand_seeds)