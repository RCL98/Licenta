from aritmetica_modulara import sqrt, inverse, extended_euclid, avg
from math import gcd, ceil
import random as rand
from sympy.ntheory import pollard_pm1, factorint
from generate_primes import logarithm_test_numbers
from timeit import default_timer as timer
from functools import partial

def pollard_p1(numar, limita = 15, reincearcari = 1, step = 1, seed = 42):
  a = 2
  rand.seed(seed)
  for _ in range(reincearcari):
    for n in range(2, limita + 1):
      a = pow(a, n, numar)
      if n % step == 0:
        div = gcd(a - 1, numar)
        if div > 1 and div < numar:
          return div
    a = rand.randint(2, numar - 2)
  return None

def pollard_rho_fact(numar, val_init_x = 2, val_init_b = 1, reincearcari = 1, seed = 42, max_iter = None, G = None):
  x = val_init_x
  b = val_init_b
  for _ in range(reincearcari):
    y = x
    if not G:
      G = lambda a: ((pow(a, 2, numar) + b)) % numar
    i = 0
    while 1:
      if max_iter and i > max_iter:
        break
      x = G(x)
      y = G(G(y))
      d = gcd(x - y, numar)
      if d == 1:
        continue
      if d == numar:
        break
      return d
      x = rand.randint(1, numar - 1)
      b = rand.randint(1, numar // 2)
      G = None
  return None

def found_match(g, h, p, a_x, a_y, b_x, b_y):
  u = (a_x - a_y) % (p - 1)
  v = (b_y - b_x) % (p - 1)
  d, s, t = extended_euclid(v, p - 1)
  if d == 1:
    return (u * s) % (p - 1)
  else:
    # print("D: ",d)
    w = ((u * s) % (p - 1)) // d
    sol_pos = [w + k * ((p - 1) // d) for k in range(d)]
    for sol in sol_pos:
      if pow(g, sol, p) == h:
        return sol

def func_pollard(x, a, b, g, h, p):
  if x < p // 3:
    return ((g * x) % p, (a + 1) % (p - 1), b)
  elif x >= p // 3 and x < 2 * p // 3:
    return ((x * x) % p, (2 * a) % (p - 1), (2 * b) % (p - 1))
  else:
    return ((h * x) % p, a, (b + 1) % (p - 1))

def func_pollard_mod(x, a, b, M, N, m, n, p):
  if x < p // 3:
    return ((M * x) % p, (m + 1) % (p - 1), b)
  elif x >= p // 3 and x < 2 * p // 3:
    return ((x * x) % p, (2 * a) % (p - 1), (2 * b) % (p - 1))
  else:
    return ((N * x) % p, a, (b + n) % (p - 1))

def func_teske(x, a, b, unit, p, Ms, ms, ns):
  s = min(x // unit, 19)
  return ((x * Ms[s]) % p, (a + ms[s]) % (p - 1), (b + ns[s]) % (p - 1))

def func_combined(x, a, b, unit, p, Ms, ms, ns):
  s = min(x // unit, 19)
  if s in Ms.keys():
    return ((x * Ms[s]) % p, (a + ms[s]) % (p - 1), (b + ns[s]) % (p - 1))
  return ((x * x) % p, (2 * a) % (p - 1), (2 * b) % (p - 1))


def floyd_cycle(x, a_x, b_x, y, a_y, b_y, g, h, p, function, max_iter):
  iterations = 0
  for l in range(max_iter):
    x, a_x, b_x = function(x, a_x, b_x)
    for _ in range(2):
      y, a_y, b_y = function(y, a_y, b_y)
    iterations += 3
    if x == y:
      sol = found_match(g, h, p, a_x, a_y, b_x, b_y)
      return sol, iterations

def brent_cycle(x, a_x, b_x, y, a_y, b_y, g, h, p, function, max_iter):
  iterations, solution, r = 0, 1, 1
  while solution == 1 and iterations < max_iter:
    x, a_x, b_x = y, a_y, b_y
    k = 0
    while k < r:
      y, a_y, b_y = function(y, a_y, b_y)
      iterations += 1
      if x == y:
        solution = found_match(g, h, p, a_x, a_y, b_x, b_y)
        return solution, iterations
      k += 1
    r *= 2
  return None

def pollard_rho_pld(g, h, p, max_iter = 10**9, cycle = 'floyd', f = "pollard", rand_seed = 0):
  rand.seed(rand_seed)
  a_0 = rand.randint(0, p - 1)
  x = pow(g, a_0, p)
  if f == 'pol':
    func = partial(func_pollard, g = g, h = h, p = p)
  else:
    m = rand.randint(0, p - 1); n = rand.randint(0, p - 1)
    M = pow(g, m, p); N = pow(g, n, p)
    func = partial(func_pollard_mod, M = M, N = N, m = m, n = n, p = p)
  if cycle == 'floyd':
    return floyd_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
  else:
    return brent_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)

def pollard_teske_rho_pld(g, h, p,  max_iter = 10**9, cycle = 'floyd', f = 'teske', rand_seed = 0):
  rand.seed(rand_seed)
  a_0 = rand.randint(0, p - 1)
  ms = rand.sample(range(1, p), 20)
  ns = rand.sample(range(1, p), 20)
  Ms = {}
  for ind, (n, m) in enumerate(zip(ns, ms)):
    Ms[ind] = (pow(g, m, p) * pow(h, n, p)) % p
  unit = p // 20
  x = pow(g, a_0, p)

  if f == 'teske':
    func = partial(func_teske, p = p, unit = unit, Ms = Ms, ms = ms, ns = ns)
  else:
    for u in rand.sample(range(1, 20), 4):
      del Ms[u]
    func = partial(func_combined, p = p, unit = unit, Ms = Ms, ms = ms, ns = ns)
  if cycle == 'floyd':
    return floyd_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
  else:
    return brent_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)

def pollard_rho(g, h, p,  max_iter = 10**10, cycle = 'floyd', func = 'pollard', rand_seed = 0):
  if f == 'pollard' or  f == 'pollard_mod':
    return pollard_rho_pld(g, h, p, max_iter, cycle, func, rand_seed)
  elif f == 'teske' or f == 'combined':
    return pollard_teske_rho_pld(g, h, p, max_iter, cycle, func, rand_seed)
  else:
    print("The" + f + "doesn't exist!")

def test_pollard(numOfTests: int, listOfBits: list, contests: dict, rand_seed = 0):
  import matplotlib.pyplot as plt
  from matplotlib.lines import Line2D
  from itertools import cycle
  rand.seed(rand_seed)
  marker = cycle(sample(Line2D.markers.keys(), len(rs) + 1))
  plt.tight_layout(pad=0.05)

  times = [[0 for _ in contests] for _ in listOfBits]
  iterations = [[0 for _ in contests] for _ in listOfBits]
  for ind_b, bits in enumerate(listOfBits):
    print("Bits: ", bits)

    for _ in range(numOfTests):
      p, g, exp, h = logarithm_test_numbers(bits)
      for ind_c, c in enumerate(contests):
        t = timer()
        x, iter = pollard_rho(g, h, p, cycle= c[0], func = c[1], rand_seed=rand_seed)
        if x != exp:
          print("WRONG " + c[0] + " " + c[1])
        times[ind_b][ind_c] += (timer() - t) * 1000
        iterations[ind_b][ind_c] += iter
    for ind_c, c in enumerate(contests):
      iterations[ind_b][ind_c] //= numOfTests
      times[ind_b][ind_c] /= numOfTests
      print(f"Cycle: {c[0]}, function: {c[1]}, itertions: {iterations[ind_b][ind_c]}, time: {imes[ind_b][ind_c]:.3f}")
    print()

  plt.figure(0)
  for ind_c, c in enumerate(contests):
    the_times = [t[ind_c] for t in times]
    plt.plot(listOfBits, the_times, marker=next(marker), label=c[0] + '_' + c[1])
  plt.legend()
  plt.ylabel("Milisecunde")
  plt.xlabel("biti")
  # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  plt.yscale('log')
  plt.savefig(f'imagini/pollard_iter.png', bbox_inches='tight')

  plt.figure(0)
  for ind_c, c in enumerate(contests):
    the_iterations = [it[ind_c] for it in iterations]
    plt.plot(listOfBits, the_iterations, marker=next(marker), label= c[0] + '_' + c[1])
  plt.legend()
  plt.ylabel("Milisecunde")
  plt.xlabel("biti")
  # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  plt.yscale('log')
  plt.savefig(f'imagini/pollard_iter.png', bbox_inches='tight')

  plt.show()

def test_pollard_rho(numOfTests: int, listOfBits: list, rand_seed = 0):
  import matplotlib.pyplot as plt
  from matplotlib.lines import Line2D
  from itertools import cycle
  plt.tight_layout(pad=0.05)
  marker = cycle(['*', 'o', 's', 'p'])

  rand.seed(rand_seed)
  floyd_times_pollard, brent_times_pollard = [], []
  floyd_times_tekse, brent_times_teske = [], []
  floyd_iterations_pollard, brent_iterations_pollard = [], []
  floyd_iterations_tekse, brent_iterations_teske = [], []
  for bits in listOfBits:
    fl_time_pol, brt_time_pol = 0, 0; fl_iter_pol, brt_iter_pol = 0, 0
    fl_time_tes, brt_time_tes = 0, 0; fl_iter_tes, brt_iter_tes = 0, 0
    print("Bits: ", bits)
    for tests in range(numOfTests):
      p, g, exp, h = logarithm_test_numbers(bits)

      t = timer()
      x, iter = pollard_rho_pld(g, h, p, rand_seed = rand_seed)
      if x != exp:
        print("WRONG FLOYD")
      fl_time_pol += (timer() - t) * 1000; fl_iter_pol += iter

      t = timer()
      x, iter = pollard_teske_rho_pld(g, h, p, rand_seed=rand_seed)
      if x != exp:
        print("WRONG FLOYD TESKE")
      fl_time_tes += (timer() - t) * 1000; fl_iter_tes += iter

      t = timer()
      x, iter = pollard_rho_pld(g, h, p, cycle = 'b', rand_seed = rand_seed)
      if x != exp:
        print("WRONG BRENT")
      brt_time_pol += (timer() - t) * 1000; brt_iter_pol += iter

      t = timer()
      x, iter = pollard_teske_rho_pld(g, h, p, cycle='b', rand_seed=rand_seed)
      if x != exp:
        print("WRONG BRENT TESKE")
      brt_time_tes += (timer() - t) * 1000; brt_iter_tes += iter

    print(f"FL iter: {fl_iter_pol //numOfTests}, time: {fl_time_pol /numOfTests:.3f}")
    print(f"FL-Teske iter: {fl_iter_tes // numOfTests}, time: {fl_time_tes / numOfTests:.3f}")
    print(f"BRT iter: {brt_iter_pol // numOfTests}, time: {brt_time_pol / numOfTests:.3f}")
    print(f"BRT-Tekse iter: {brt_iter_tes // numOfTests}, time: {brt_time_tes / numOfTests:.3f}\n")
    floyd_times_pollard.append(fl_time_pol/numOfTests); floyd_times_tekse.append(fl_time_tes/numOfTests);
    brent_times_pollard.append(brt_time_pol/numOfTests); brent_times_teske.append(brt_time_tes/numOfTests)
    floyd_iterations_pollard.append(fl_iter_pol//numOfTests); floyd_iterations_tekse.append(fl_iter_tes//numOfTests)
    brent_iterations_pollard.append(brt_iter_pol//numOfTests); brent_iterations_teske.append(brt_iter_tes//numOfTests)


  print(f"PF timp: {avg(floyd_times_pollard):.3f}, iter: {avg(floyd_iterations_pollard):.2e}")
  print(f"PB timp: {avg(brent_times_pollard):.3f}, iter: {avg(brent_iterations_teske):.2e}")
  print(f"TF timp: {avg(floyd_times_tekse):.3f}, iter: {avg(floyd_iterations_tekse):.2e}")
  print(f"TB timp: {avg(brent_times_teske):.3f}, iter: {avg(brent_iterations_teske):.2e}")


  plt.figure(0)
  plt.plot(listOfBits, floyd_times_pollard, marker=next(marker), label="PF")
  plt.plot(listOfBits, brent_times_pollard, marker=next(marker), label="PB")
  plt.plot(listOfBits, floyd_times_tekse, marker=next(marker), label="TF")
  plt.plot(listOfBits, brent_times_teske, marker=next(marker), label="TB")
  plt.legend()
  plt.ylabel("Milisecunde")
  plt.xlabel("biti")
  plt.yscale('log')
  # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  plt.savefig(f'imagini/pollard_timp.png', bbox_inches='tight')

  plt.figure(1)
  plt.plot(listOfBits, floyd_iterations_pollard, marker=next(marker), label="PF")
  plt.plot(listOfBits, brent_iterations_pollard, marker=next(marker), label="PB")
  plt.plot(listOfBits, floyd_iterations_tekse, marker=next(marker), label="TF")
  plt.plot(listOfBits, brent_iterations_teske, marker=next(marker), label="TB")
  plt.legend()
  plt.ylabel("Milisecunde")
  plt.xlabel("biti")
  # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  plt.yscale('log')
  plt.savefig(f'imagini/pollard_iter.png', bbox_inches='tight')
  plt.show()


# if __name__ == "__main__":
#   test_pollard_rho(10, [10, 15, 20, 25, 30, 35, 40, 45, 50, 55], 312)


