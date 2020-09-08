from aritmetica_modulara import sqrt, extended_euclid, avg
from math import gcd, ceil, floor
import random as rand
from functools import partial
from  sys import getsizeof

A = (sqrt(5) - 1)/2

def found_match(g, h, p, a_x, a_y, b_x, b_y):
  u = (a_x - a_y) % (p - 1)
  v = (b_y - b_x) % (p - 1)
  if v == 0:
    return None
  d, s, t = extended_euclid(v, p - 1)
  if d == 1:
    return (u * s) % (p - 1)
  else:
    # print("D: ",d)
    w = ((u * s) % (p - 1)) // d
    sol_pos = [w] + [(w + k * ((p - 1) // d)) % (p - 1) for k in range(1, d)]
    for sol in sol_pos:
      if pow(g, sol, p) == h:
        return sol
  return None

def hash_function(num, r):
  return floor(((A * num) % 1) * r) + 1

def func_pollard(x, a, b, g, h, p):
  s = hash_function(x, 3)
  if s == 1: #x < p // 3:
    return ((g * x) % p, (a + 1) % (p - 1), b)
  elif s == 2: #x >= p // 3 and x < 2 * p // 3:
    return ((x * x) % p, (2 * a) % (p - 1), (2 * b) % (p - 1))
  else:
    return ((h * x) % p, a, (b + 1) % (p - 1))

def func_pollard_mod(x, a, b, M, N, m, n, p):
  s = hash_function(x, 3)
  if s == 1: #x < p // 3:
    return ((M * x) % p, (a + m) % (p - 1), b)
  elif s == 2: #x >= p // 3 and x < 2 * p // 3:
    return ((x * x) % p, (2 * a) % (p - 1), (2 * b) % (p - 1))
  else:
    return ((N * x) % p, a, (b + n) % (p - 1))

def func_teske(x, a, b, unit, p):
  s = hash_function(x, 20) - 1 #min(x // unit, 19)
  return ((x * Ms[s]) % p, (a + ms[s]) % (p - 1), (b + ns[s]) % (p - 1))

def func_combined(x, a, b, unit, p):
  s = hash_function(x, 20) - 1 #min(x // unit, 19)
  if s in Ms.keys():
    return ((x * Ms[s]) % p, (a + ms[s]) % (p - 1), (b + ns[s]) % (p - 1))
  return ((x * x) % p, (2 * a) % (p - 1), (2 * b) % (p - 1))

def floyd_cycle(x, a_x, b_x, y, a_y, b_y, g, h, p, function, max_iter):
  for iter in range(max_iter):
    x, a_x, b_x = function(x, a_x, b_x)
    for _ in range(2):
      y, a_y, b_y = function(y, a_y, b_y)
    if x == y:
      sol = found_match(g, h, p, a_x, a_y, b_x, b_y)
      return sol, iter * 3

def brent_cycle(x, a_x, b_x, y, a_y, b_y, g, h, p, function, max_iter):
  iterations, solution, r = 0, 1, 1
  init, init_a, init_b = x, a_x, a_y
  while solution == 1 and iterations < max_iter:
    x, a_x, b_x = y, a_y, b_y
    k = 0
    while k < r:
      y, a_y, b_y = function(y, a_y, b_y)
      iterations += 1
      if x == y:
        lamb = 1
        tortoise, (hare, a_hare, b_hare) = x, function(x, a_x, a_y)
        while tortoise != hare:
          hare, a_hare, b_hare = function(hare, a_hare, b_hare)
          lamb += 1

        tortoise, a_tort, b_tort = init, init_a, init_b
        (hare, a_hare, b_hare) = init, init_a, init_b
        for i in range(lamb):
          (hare, a_hare, b_hare) = function(hare, a_hare, b_hare)
        mu = 0
        while tortoise != hare:
          tortoise, a_tort, b_tort = function(tortoise, a_tort, b_tort)
          hare, a_hare, b_hare = function(hare, a_hare, b_hare)
          mu += 1

        solution = found_match(g, h, p, a_x, a_y, b_x, b_y)
        return solution, iterations, lamb, mu
      k += 1
    r *= 2
  return None

def teske_cycle(x, a_x, b_x, y, a_y, b_y, g, h, p, function, max_iter):
  L = [(x, a_x, b_x)] * 8
  it, iterations = 0, 0
  d = 0
  while iterations < max_iter:
    while 1:
      x, a_x, b_x = function(x, a_x, b_x)
      match = next(iter([i for i in L if i[0] == x]), None)
      if match:
        solution = found_match(g, h, p, a_x, match[1], b_x, match[2])
        return solution, iterations
      it += 1; iterations += 1
      if it >= 3 * d:
        break
    for i in range(7):
      L[i] = L[i + 1]
    L[7] = (x, a_x, b_x); d = it; it = 0
  return None

def teske_cycle_special(x, a_x, b_x, y, a_y, b_y, g, h, p, function, max_iter):
  L = [(x, a_x, b_x)] * 8
  it, iterations = 0, 0
  d = 0
  init, init_a, init_b = x, a_x, a_y
  while iterations < max_iter:
    while 1:
      x, a_x, b_x = function(x, a_x, b_x)
      match = next(iter([i for i in L if i[0] == x]), None)
      if match:
        lamb = 1
        tortoise, (hare, a_hare, b_hare) = match[0], function(match[0], match[1], match[2])
        while tortoise != hare:
          hare, _, _ = function(hare, 0, 0)
          lamb += 1

        tortoise = init
        hare = init
        for i in range(lamb):
          (hare, a_hare, b_hare) = function(hare, 0, 0)
        mu = 0
        while tortoise != hare:
          tortoise, a_tort, b_tort = function(tortoise, 0, 0)
          hare, a_hare, b_hare = function(hare, 0, 0)
          mu += 1

        solution = found_match(g, h, p, a_x, match[1], b_x, match[2])
        return solution, iterations, lamb, mu
      it += 1; iterations += 1
      if it >= 3 * d:
        break
    for i in range(7):
      L[i] = L[i + 1]
    L[7] = (x, a_x, b_x); d = it;
  return None

def pollard_rho_pld(g, h, p, max_iter = 10**12, cycle = 'floyd', f = "pollard", rand_seed = 0):
  rand.seed(rand_seed)
  a_0 = rand.randint(0, p - 1)
  x = pow(g, a_0, p)
  if f == 'pollard':
    func = partial(func_pollard, g = g, h = h, p = p)
  else:
    m = rand.randint(0, p - 1); n = rand.randint(0, p - 1)
    M = pow(g, m, p); N = pow(h, n, p)
    func = partial(func_pollard_mod, M = M, N = N, m = m, n = n, p = p)
  if cycle == 'floyd':
    rez = floyd_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
    while rez[0] == None:
      new_seed = rand.randint(0, max(p - 1, rand_seed))
      return pollard_rho_pld(g, h, p, max_iter, cycle, f, new_seed)
    return rez
  elif cycle == 'brent':
    rez = brent_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
    while rez[0] == None:
      new_seed = rand.randint(0, max(p - 1, rand_seed))
      return pollard_rho_pld(g, h, p, max_iter, cycle, f, new_seed)
    return rez
  elif cycle == 'teske':
    rez = teske_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
    while rez[0] == None:
      new_seed = rand.randint(0, max(p - 1, rand_seed))
      return pollard_rho_pld(g, h, p, max_iter, cycle, f, new_seed)
    return rez
  elif cycle == 'teske_spec':
    rez = teske_cycle_special(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
    while rez[0] == None:
      new_seed = rand.randint(0, max(p - 1, rand_seed))
      return pollard_rho_pld(g, h, p, max_iter, cycle, f, new_seed)
    return rez

def pollard_teske_rho_pld(g, h, p,  max_iter = 10**12, cycle = 'floyd', f = 'teske', rand_seed = 0):
  rand.seed(rand_seed)
  a_0 = rand.randint(0, p - 1)
  global ms, ns, Ms
  ms = rand.sample(range(1, p), 20)
  ns = rand.sample(range(1, p), 20)
  Ms = {}
  for ind, (n, m) in enumerate(zip(ns, ms)):
    Ms[ind] = (pow(g, m, p) * pow(h, n, p)) % p
  unit = p // 20
  x = pow(g, a_0, p)

  if f == 'teske':
    func = partial(func_teske, p = p, unit = unit)
  else:
    for u in rand.sample(range(1, 20), 4):
      del Ms[u]
    func = partial(func_combined, p = p, unit = unit)
  if cycle == 'floyd':
    rez = floyd_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
    while rez[0] == None:
      new_seed = rand.randint(0, max(p - 1, rand_seed))
      return pollard_teske_rho_pld(g, h, p,  max_iter, cycle, f, new_seed)
    return rez
  elif cycle == 'brent':
    rez = brent_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
    while rez[0] == None:
      new_seed = rand.randint(0, max(p - 1, rand_seed))
      return pollard_teske_rho_pld(g, h, p, max_iter, cycle, f, new_seed)
    return rez
  elif cycle == 'teske':
    rez = teske_cycle(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
    while rez[0] == None:
      new_seed = rand.randint(0, max(p - 1, rand_seed))
      return pollard_teske_rho_pld(g, h, p, max_iter, cycle, f, new_seed)
    return rez
  elif cycle == 'teske_spec':
    rez = teske_cycle_special(x, a_0, 0, x, a_0, 0, g, h, p, func, max_iter)
    while rez[0] == None:
      new_seed = rand.randint(0, max(p - 1, rand_seed))
      return pollard_teske_rho_pld(g, h, p, max_iter, cycle, f, new_seed)
    return rez

def pollard_rho(g, h, p,  max_iter = 10**12, cycle = 'floyd', func = 'pollard', rand_seed = 0):
  if func == 'pollard' or  func == 'pollard_mod':
    return pollard_rho_pld(g, h, p, max_iter, cycle, func, rand_seed)
  elif func == 'teske' or func == 'combined':
    return pollard_teske_rho_pld(g, h, p, max_iter, cycle, func, rand_seed)
  else:
    print("The" + f + "doesn't exist!")



def pollard_simple(g, h, p,  max_iter = 10**10, rand_seed = 0):
  rand.seed(rand_seed)
  a_0 = rand.randint(0, p - 1)
  x = pow(g, a_0, p)
  a_x, b_x = a_0, 0
  L = {x: (a_x, b_x)}
  dict_limit = 1.15 * (10 ** 10)
  for iter in range(max_iter):
    if getsizeof(L) > dict_limit:
      return None, None
    x, a_x, b_x = func_pollard(x, a_x, b_x, g, h, p)
    y = L.get(x, None)
    if y:
      a_y, b_y = y
      u = (a_x - a_y) % (p - 1)
      v = (b_y - b_x) % (p - 1)
      if v == 0:
        new_seed = rand.randint(a_0, p - 1)
        return pollard_floyd(g, h, p, max_iter, new_seed)
      d, s, t = extended_euclid(v, p - 1)
      if d == 1:
        return (u * s) % (p - 1), len(L)
      else:
        w = ((u * s) % (p - 1)) // d
        sol_pos = [w] + [(w + k * ((p - 1) // d)) % (p - 1) for k in range(1, d)]
        for sol in sol_pos:
          if pow(g, sol, p) == h:
            return sol, len(L)
    L[x] = (a_x, b_x)
  return None

def pollard_floyd(g, h, p,  max_iter = 10**10, rand_seed = 0):
  rand.seed(rand_seed)
  a_0 = rand.randint(0, p - 1)
  x = pow(g, a_0, p)
  y = x
  a_x, a_y = a_0, a_0
  b_x, b_y = 0, 0
  for iter in range(max_iter):
    x, a_x, b_x = func_pollard(x, a_x, b_x, g, h, p)
    for _ in range(2):
      y, a_y, b_y = func_pollard(y, a_y, b_y, g, h, p)
    if x == y:
      u = (a_x - a_y) % (p - 1)
      v = (b_y - b_x) % (p - 1)
      if v == 0:
        new_seed = rand.randint(a_0, p - 1)
        return pollard_floyd(g, h, p, max_iter, new_seed)
      d, s, t = extended_euclid(v, p - 1)
      if d == 1:
        return (u * s) % (p - 1), 3 * iter
      else:
        w = ((u * s) % (p - 1)) // d
        sol_pos = [w] + [(w + k * ((p - 1) // d)) % (p - 1) for k in range(1, d)]
        for sol in sol_pos:
          if pow(g, sol, p) == h:
            return sol,  3 * iter
  return None

def pollard_brent(g, h, p,  max_iter = 10**10, rand_seed = 0):
  rand.seed(rand_seed)
  a_0 = rand.randint(0, p - 1)
  x = pow(g, a_0, p)
  y = x
  iterations, solution, r = 0, 1, 1
  a_x, a_y = a_0, a_0
  b_x, b_y = 0, 0
  while iterations < max_iter:
    x, a_x, b_x = y, a_y, b_y
    k = 0
    while k < r:
      y, a_y, b_y = func_pollard(y, a_y, b_y, g, h, p)
      iterations += 1
      if x == y:
        u = (a_x - a_y) % (p - 1)
        v = (b_y - b_x) % (p - 1)
        if v == 0:
          new_seed = rand.randint(a_0, p - 1)
          return pollard_brent(g, h, p, max_iter, new_seed)
        d, s, t = extended_euclid(v, p - 1)
        if d == 1:
          return (u * s) % (p - 1), iterations
        else:
          w = ((u * s) % (p - 1)) // d
          sol_pos = [w] + [(w + k * ((p - 1) // d)) % (p - 1) for k in range(1, d)]
          for sol in sol_pos:
            if pow(g, sol, p) == h:
              return sol, iterations
      k += 1
    r *= 2
  return None
