from aritmetica_modulara import sqrt
from sympy.ntheory import factorint
from math import ceil, log2
from functools import reduce
import operator
import generate_primes as gnp
from timeit import default_timer as timer
from random import randint, seed
import numpy as np
from sys import getsizeof
import bisect

def shanks_ordin_necunoscut(g, h, p):
  rad_p = ceil(sqrt(p - 1))
  n = 100
  lista_n = []
  while n < rad_p:
    lista_n.append(n)
    n *= 10
  lista_n.append(rad_p)
  last_n = 1
  e = 1
  lista = {e: 0}
  for n in lista_n:
    for j in range(last_n, n + 1):
      e = (e * g) % p
      if e == h:
        return j
      lista[e] = j

    u = pow(g, p - n - 1, p)
    l = h
    for i in range(0, n + 1):
      try:
        j = lista[l]
        return i * n + j
      except KeyError:
        l = (l * u) % p
    last_n = n + 1
  return None

def shanks_classic(g, h, p, ordin):
  n = ceil(sqrt(ordin))
  small_steps, giant_steps = 1, 0
  e = 1
  lista = {e: 0}
  for j in range(1, n):
    e = (e * g) % p
    if e == h:
      return j, small_steps, giant_steps, getsizeof(lista)
    else:
      lista[e] = j
      small_steps += 1
  u = pow(g, p - n - 1, p)
  s = pow(g, n, p)
  ep = em = h
  u = pow(g, p - n - 1, p)
  e = h
  giant_steps = 1
  for i in range(0, n):
    try:
      j = lista[e]
      return i * n  + j, small_steps, giant_steps, getsizeof(lista)
    except KeyError:
      e = (e * u) % p
      giant_steps += 1
  return None

def find_index(elements, value):
  left, right = 0, len(elements) - 1

  while left <= right:
    middle = (left + right) // 2
    middle_element = elements[middle][0]

    if middle_element == value:
      return elements[middle][1]

    if middle_element < value:
      left = middle + 1
    elif middle_element > value:
      right = middle - 1
  return None

def shanks_classic_binary_search(g, h, p, ordin):
  n = ceil(sqrt(ordin))
  small_steps, giant_steps = 1, 0
  e = 1
  lista = [(e,0)]
  for j in range(1, n):
    e = (e * g) % p
    if e == h:
      return j, small_steps, giant_steps, getsizeof(lista)
    else:
      lista.append((e,j))
      small_steps += 1
  lista.sort(key=lambda x: x[0])
  u = pow(g, p - n - 1, p)
  s = pow(g, n, p)
  ep = em = h
  u = pow(g, p - n - 1, p)
  e = h
  giant_steps = 1
  for i in range(0, n):
    j = find_index(lista, e)
    if j != None:
      return i * n  + j, small_steps, giant_steps, getsizeof(lista)
    else:
      e = (e * u) % p
      giant_steps += 1
  return None

def shanks_classic_with_memory(g, h_values, p, ordin):
  n = ceil(sqrt(ordin))
  small_steps = 1
  e = 1
  lista = {e: 0}
  for j in range(1, n):
    e = (e * g) % p
    lista[e] = j
    small_steps += 1

  u = pow(g, p - n - 1, p)
  exponents, giant_steps_list = [], []
  for h in h_values:
    e = h
    giant_steps = 1
    for i in range(0, n):
      try:
        j = lista[e]
        exponents.append(i * n + j)
        giant_steps_list.append(giant_steps)
        break
      except KeyError:
        e = (e * u) % p
        giant_steps += 1
  return exponents, small_steps, giant_steps_list

def shanks_general(g, h, p, ordin, r):
  n = ceil(ordin ** (1/r))
  m = ceil(ordin/n)
  small_steps, giant_steps = 1, 0
  e = 1
  lista = {e: 0}
  for j in range(1, m):
    e = (e * g) % p
    if e == h:
      return j, small_steps, giant_steps
    else:
      lista[e] = j
      small_steps += 1

  u = pow(g, p - m - 1, p)
  e = h
  giant_steps = 1
  for i in range(0, n):
    try:
      j = lista[e]
      return i * m + j, small_steps, giant_steps
    except KeyError:
      e = (e * u) % p
      giant_steps += 1
  return None

def shanks_general_with_memory(g, h_values, p, ordin, r):
  n = ceil(ordin ** (1/r))
  m = ceil(ordin/n)
  small_steps = 1
  e = 1
  lista = {e: 0}
  for j in range(1, m):
    e = (e * g) % p
    lista[e] = j
    small_steps += 1
  u = pow(g, p - m - 1, p)
  exponents, giant_steps_list = [], []
  for h in h_values:
    e = h
    giant_steps = 1
    for i in range(0, n):
      try:
        j = lista[e]
        exponents.append(i * m + j)
        giant_steps_list.append(giant_steps)
        break
      except KeyError:
        e = (e * u) % p
        giant_steps += 1
  return exponents, small_steps, giant_steps_list

def optim_factors(ordin):
  factori = [pow(x, e) for x, e in factorint(ordin).items()]
  factori.sort()
  if len(factori) == 1:
    return Shanks_ordin_cunoscut(g, h, p, ordin)
  if len(factori) == 2:
    l = ceil(sqrt(factori[1]))
    m = ceil(sqrt(factori[0]))
  else:
    fact_1 = reduce(operator.mul, factori[:-1])
    if fact_1 > factori[-1]:
      l = ceil(sqrt(fact_1))
      m = ceil(sqrt(factori[-1]))
    else:
      l = ceil(sqrt(factori[-1]))
      m = ceil(sqrt(fact_1))
  return l, m

def shanks_factor(g, h, p, l, m):
  small_steps, giant_steps = 1, 0
  e = 1
  lista = {e: 0}
  for j in range(1, l ** 2):
    e = (e * g) % p
    if e == h:
      return j, small_steps, giant_steps
    else:
      lista[e] = j
      small_steps += 1

  u = pow(g, p - l ** 2 - 1, p)
  e = h
  giant_steps = 1
  for i in range(0, m ** 2):
    try:
      j = lista[e]
      return i * (l ** 2) + j, small_steps, giant_steps
    except KeyError:
      e = (e * u) % p
      giant_steps += 1
  return None

def shanks_factor_with_memory(g, h_values, p, l, m):
  small_steps = 1
  e = 1
  lista = {e: 0}
  for j in range(1, l ** 2):
    e = (e * g) % p
    lista[e] = j
    small_steps += 1
  u = pow(g, p - l ** 2 - 1, p)
  exponents, giant_steps_list = [], []
  for h in h_values:
    e = h
    giant_steps = 1
    for i in range(0, m ** 2):
      try:
        j = lista[e]
        exponents.append(i * (l ** 2) + j)
        giant_steps_list.append(giant_steps)
        break
      except KeyError:
        e = (e * u) % p
        giant_steps += 1
  return exponents, small_steps, giant_steps_list

def shanks_middle_to_edges(g, h, p, ordin):
  n = ceil(sqrt(ordin))
  K = ordin // 2
  small_steps, giant_steps = 1, 0
  e = pow(g, K, p)
  lista = {e: 0}
  for j in range(1, n):
    e = (e * g) % p
    if e == h:
      return j, small_steps, giant_steps
    else:
      lista[e] = j
      small_steps += 1
  u = pow(g, p - n - 1, p)
  s = pow(g, n, p)
  ep = em = h
  giant_steps = 2
  for i in range(0, n):
    try:
      j = lista[em]
      return j + K - i * n, small_steps, giant_steps
    except KeyError:
      try:
        j = lista[ep]
        return K + j + i * n, small_steps, giant_steps
      except KeyError:
        ep = (ep * u) % p
        em = (em * s) % p
        giant_steps += 2
  return None

def Shanks(g, h, p, ordin = None):
  if ordin:
    return Shanks_classic(g, h, p, ordin)
  return Shanks_ordin_necunoscut(g, h, p)

def avg(lst):
    return reduce(lambda a, b: a + b, lst) / len(lst)

def extrapolate(data, numOfPoints):
  ratio = 0
  for i in range(1, len(data)):
    ratio += data[i] // data[i - 1]
  ratio = ratio//(len(data) - 1)
  newPoints = [ratio * data[-1]]
  for _ in range(1, numOfPoints):
    newPoints.append(ratio * newPoints[-1])
  return newPoints

def shakns_benchmark(numOfTests: int, maxBits: int, rand_seed = 0):
    import matplotlib.pyplot as plt
    seed(rand_seed)
    times, small_steps_list, giant_steps_list, total_steps, list_sizes, bits, time_bin, list_sizes_bin = [], [], [], [], [], [], [], []
    for numOfBits in range(10, maxBits + 5, 5):
      time, small, giant, total, list_sz, tb, lb = 0, 0, 0, 0, 0, 0, 0
      bits.append(numOfBits)
      for _ in range(numOfTests):
        p, g, e, h = gnp.logarithm_test_numbers(numOfBits)
        t_start = timer()
        x, small_steps, giant_steps, list_size = shanks_classic(g, h, p, p - 1)
        time += (timer() - t_start)*1000
        small += small_steps; giant += giant_steps; total += small_steps + giant_steps; list_sz += list_size
        t_start = timer()
        x, small_steps, giant_steps, list_size = shanks_classic_binary_search(g, h, p, p - 1)
        tb += (timer() - t_start) * 1000; lb += list_size
      times.append(time/numOfTests)
      time_bin.append(tb/numOfTests)
      small_steps_list.append(small/numOfTests)
      giant_steps_list.append(giant/numOfTests)
      total_steps.append(total/numOfTests)
      list_sizes.append(list_sz/numOfTests)
      list_sizes_bin.append(lb/numOfTests)

    next_xs = list(range(numOfBits + 5, numOfBits + 25, 5))
    plt.figure(1)
    plt.xlim(5, numOfBits + 25)
    next_ys = extrapolate(times, len(next_xs))
    next_ys_bin = extrapolate(time_bin, len(next_xs))
    plt.plot(np.append(bits, next_xs), np.append(times, next_ys), marker = 'o', color = 'b', linestyle = 'dashed', label = 'dispersie')
    plt.plot(bits, times, marker = 'o', c='b')
    plt.plot(np.append(bits, next_xs), np.append(time_bin, next_ys_bin), marker = 's', c ='r', linestyle = 'dashed', label ='cautare binara')
    plt.plot(bits, time_bin, marker = 's', c='r')
    plt.xlabel('biti')
    plt.ylabel('milisecunde')
    plt.yscale('log')
    plt.legend(loc="upper left")
    plt.grid(True)
    plt.savefig('imagini/shanks_classic_timp.png')

    plt.figure(2)
    plt.xlim(5, numOfBits + 25)
    next_ys = extrapolate(small_steps_list, len(next_xs))
    plt.plot(np.append(bits, next_xs), np.append(small_steps_list, next_ys), 'b--')
    plt.plot(bits, small_steps_list, c='b')
    plt.xlabel('biti')
    plt.ylabel('pasi mici')
    plt.yscale('log')
    plt.grid(True)
    plt.savefig('imagini/shanks_classic_pasi_mici.png')

    plt.figure(3)
    plt.xlim(5, numOfBits + 25)
    next_ys = extrapolate(giant_steps_list, len(next_xs))
    plt.plot(np.append(bits, next_xs), np.append(giant_steps_list, next_ys), 'b--')
    plt.plot(bits, giant_steps_list, c='b')
    plt.xlabel('biti')
    plt.ylabel('pasi giant')
    plt.yscale('log')
    plt.grid(True)
    plt.savefig('imagini/shanks_classic_pasi_mari.png')

    plt.figure(4)
    plt.xlim(5, numOfBits + 25)
    next_ys = extrapolate(total_steps, len(next_xs))
    plt.plot(np.append(bits, next_xs), np.append(total_steps, next_ys), 'b--')
    plt.plot(bits, total_steps, c='b')
    plt.xlabel('biti')
    plt.ylabel('pasi totali')
    plt.yscale('log')
    plt.grid(True)
    plt.savefig('imagini/shanks_classic_pasi_total.png')

    plt.figure(5)
    plt.xlim(5, numOfBits + 25)
    next_ys = extrapolate(list_sizes, len(next_xs))
    next_ys_bin = extrapolate(list_sizes_bin, len(next_xs))
    plt.plot(np.append(bits, next_xs), np.append(list_sizes, next_ys),  marker = 'o', color = 'b', linestyle = 'dashed', label = 'dispersie')
    plt.plot(bits, list_sizes, marker = 'o', c='b')
    plt.plot(np.append(bits, next_xs), np.append(list_sizes_bin, next_ys_bin), marker = 's', c ='r', linestyle = 'dashed', label ='cautare binara')
    plt.plot(bits, list_sizes_bin, marker = 's', c='r')
    plt.xlabel('biti')
    plt.ylabel('bytes')
    plt.yscale('log')
    plt.legend(loc="upper left")
    plt.grid(True)
    plt.savefig('imagini/shanks_classic_memory.png')

    plt.show()

def test_shanks(numOfTests: int, numOfBits: int, type: str, r = 2):
  times = []
  for _ in range(numOfTests):
    p, g, e, h = gnp.logarithm_test_numbers(1, numOfBits)
    if type == 'c':
      t_start = time.time()
      x, small_steps, giant_steps = Shanks_ordin_cunoscut(g[0], h[0], p[0], p[0] - 1)
    elif type == 'g':
      t_start = time.time()
      x, small_steps, giant_steps = Shanks_general(g[0], h[0], p[0], p[0] - 1, r)
    else:
      l, m = optim_factors(p[0] - 1)
      t_start = time.time()
      x, small_steps, giant_steps = Shanks_factor(g[0], h[0], p[0], l, m)
    t_stop = time.time()
    if e[0] != x:
      print(f"e:{e[0]}  x:{x}")
    else:
      t = (t_stop - t_start)/1000
      times.append(t)
      print(f"Test passed! p:{p[0]}, time:{t:.7f} ms, small_steps:{small_steps} giant_steps:{giant_steps}")
  print(f"Avg time:{avg(times)}")

def test_shanks_same_p(numOfTests: int, numOfBits: int, r = 2, radn_seed = 0):
  seed(radn_seed)
  times_c, times_g, times_f = [], [], []
  p, g, exps, hs = gnp.logarithm_test_numbers_same_p(numOfTests, numOfBits)
  l, m = optim_factors(p - 1)

  t_start_c = timer()
  x_cs, small_steps_c, giant_steps_c = shanks_classic_with_memory(g, hs, p, p - 1)
  t_stop_c = timer()


  t_start_g = timer()
  x_gs, small_steps_g, giant_steps_g = shanks_general_with_memory(g, hs, p, p - 1, r)
  t_stop_g = timer()


  t_start_f = timer()
  x_fs, small_steps_f, giant_steps_f = shanks_factor_with_memory(g, hs, p, l, m)
  t_stop_f = timer()

  for e, x_c, x_g, x_f in zip(exps, x_cs, x_gs, x_fs):
    if x_c != e:
      print(f"e:{e}  x_c:{x_c}")
    if x_g != e:
      print(f"e:{e}  x_g:{x_g}")
    if x_f != e:
      print(f"e:{e}  x_f:{x_f}")

  print(f"t_c: {(t_stop_c - t_start_c):.7f}, t_g: {(t_stop_g - t_start_g)}, t_f:{(t_stop_f - t_start_f)}")
  print(f"sm_c: {small_steps_c}, sm_g: {small_steps_g}, sm_f:{small_steps_f}")
  print()
  for (e, g_c, g_g, g_f) in zip(exps, giant_steps_c, giant_steps_g, giant_steps_f):
    print(f"e:{e}, giant_steps_c:{g_c}, giant_steps_g:{g_g} giant_steps_f:{g_f}")
  print()
  print(f"giant_steps_c_avg:{avg(giant_steps_c)}, giant_steps_g:{avg(giant_steps_g)}, giant_steps_f:{avg(giant_steps_f)}")

def test_shanks_middle(numOfTests: int, numOfBits: int, randSeed = 0):
  seed(randSeed)
  giant_steps_c_list, small_steps_c_list = [], []
  giant_steps_m_list, small_steps_m_list = [], []
  for iter in range(numOfTests):
    p, g, exp, h = gnp.logarithm_test_numbers(numOfBits, True)
    print("iter:", iter)
    x_c, small_steps_c, giant_steps_c = shanks_classic(g, h, p, p - 1)
    if x_c != exp:
      print(f"Classic x:{x_c}, exp:{exp}")
    else:
      giant_steps_c_list.append(giant_steps_c)
      small_steps_c_list.append(small_steps_c)
      print(f"C_passed x:{x_c}, exp:{exp}, sm:{small_steps_c}, gs:{giant_steps_c}")
    x_m, small_steps_m, giant_steps_m = shanks_middle_to_edges(g, h, p, p - 1)
    if x_m != exp:
      print(f"Middle x:{x_m}, exp:{exp}\n")
    else:
      giant_steps_m_list.append(giant_steps_m)
      small_steps_m_list.append(small_steps_m)
      print(f"M_passed x:{x_c}, exp:{exp}, sm:{small_steps_m}, gs:{giant_steps_m}\n")
  print(f"Classic: smavg:{avg(small_steps_c_list)}, gsavg:{avg(giant_steps_c_list)}")
  print(f"Middle: smavg:{avg(small_steps_m_list)}, gsavg:{avg(giant_steps_m_list)}")

# if __name__ == "__main__":
# #   from memory_profiler import memory_usage
#     shakns_benchmark(30, 30, 802)
# #     for _ in range(100):
# #       p, g, e, h = gnp.logarithm_test_numbers(25)
# #       x, sm, gm, ls = shanks_classic_binary_search(g, h, p, p-1)
# #       if x != e:
# #         print(x,e)
# #   # mem_usage = memory_usage((shanks_classic, (g, h, p, p - 1)), max_iterations=1)
# #   # print('Memory usage (in chunks of .1 seconds): %s' % mem_usage)
# #   # print('Maximum memory usage: %s' % sum(mem_usage))
