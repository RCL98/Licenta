from aritmetica_modulara import sqrt
from sympy.ntheory import factorint
from math import ceil
from functools import reduce
from operator import mul
from sys import getsizeof

def shanks_order_unkown(g, h, p):
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
  e = 1
  lista = {e: 0}
  dict_limit = 1.15 * (10 ** 10)
  for j in range(1, n):
    if getsizeof(lista) > dict_limit:
      return None, None, None, None
    e = (e * g) % p
    if e == h:
      return j, j, 0, getsizeof(lista)
    else:
      lista[e] = j
  u = pow(g, p - n - 1, p)
  e = h
  for i in range(0, n):
    try:
      j = lista[e]
      return i * n  + j, n, i + 1, getsizeof(lista)
    except KeyError:
      e = (e * u) % p
  return None

def shanks_classic_better_avg(g, h, p, ordin):
  n = ceil(sqrt(ordin/2))
  e = 1
  lista = {e: 0}
  for j in range(1, n):
    e = (e * g) % p
    if e == h:
      return j, j, 0, getsizeof(lista)
    else:
      lista[e] = j
  u = pow(g, p - n - 1, p)
  e = h
  for i in range(0, 2 * n):
    try:
      j = lista[e]
      return i * n  + j, n, i + 1, getsizeof(lista)
    except KeyError:
      e = (e * u) % p
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
  small_steps = 1
  e = 1
  lista = [(e,0)]
  for j in range(1, n):
    e = (e * g) % p
    if e == h:
      return j, small_steps, 0, getsizeof(lista)
    else:
      lista.append((e,j))
      small_steps += 1
  lista.sort(key=lambda x: x[0])
  u = pow(g, p - n - 1, p)
  e = h
  for i in range(0, n):
    j = find_index(lista, e)
    if j != None:
      return i * n  + j, small_steps, i + 1, getsizeof(lista)
    else:
      e = (e * u) % p
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
    for i in range(0, n):
      try:
        j = lista[e]
        exponents.append(i * n + j)
        giant_steps_list.append(i + 1)
        break
      except KeyError:
        e = (e * u) % p
  return exponents, small_steps, giant_steps_list, getsizeof(lista)

def shanks_general(g, h, p, ordin, r):
  n = ceil(ordin ** (1/r))
  m = ceil(ordin/n)
  small_steps = 1
  e = 1
  lista = {e: 0}
  for j in range(1, m):
    e = (e * g) % p
    if e == h:
      return j, small_steps, 0
    else:
      lista[e] = j
      small_steps += 1

  u = pow(g, p - m - 1, p)
  e = h
  for i in range(0, n):
    try:
      j = lista[e]
      return i * m + j, small_steps, i + 1
    except KeyError:
      e = (e * u) % p
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
    for i in range(0, n):
      try:
        j = lista[e]
        exponents.append(i * m + j)
        giant_steps_list.append(i + 1)
        break
      except KeyError:
        e = (e * u) % p
  return exponents, small_steps, giant_steps_list, getsizeof(lista)

def optim_factors(ordin):
  factori = [pow(x, e) for x, e in factorint(ordin).items()]
  factori.sort()
  if len(factori) == 1:
    return Shanks_ordin_cunoscut(g, h, p, ordin)
  if len(factori) == 2:
    l = ceil(sqrt(factori[1]))
    m = ceil(sqrt(factori[0]))
  else:
    fact_1 = reduce(mul, factori[:-1])
    if fact_1 > factori[-1]:
      l = ceil(sqrt(fact_1))
      m = ceil(sqrt(factori[-1]))
    else:
      l = ceil(sqrt(factori[-1]))
      m = ceil(sqrt(fact_1))
  return l, m

def shanks_factor(g, h, p, l, m):
  small_steps = 1
  e = 1
  lista = {e: 0}
  for j in range(1, l ** 2):
    e = (e * g) % p
    if e == h:
      return j, small_steps, 0
    else:
      lista[e] = j
      small_steps += 1

  u = pow(g, p - l ** 2 - 1, p)
  e = h
  for i in range(0, m ** 2):
    try:
      j = lista[e]
      return i * (l ** 2) + j, small_steps, i + 1
    except KeyError:
      e = (e * u) % p
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
    for i in range(0, m ** 2):
      try:
        j = lista[e]
        exponents.append(i * (l ** 2) + j)
        giant_steps_list.append(i + 1)
        break
      except KeyError:
        e = (e * u) % p
  return exponents, small_steps, giant_steps_list, getsizeof(lista)

def shanks_centered(g, h, p, ordin):
  n = ceil(sqrt(ordin))
  K = ordin // 2
  small_steps = 1
  e = pow(g, K, p)
  lista = {e: 0}
  for j in range(1, n):
    e = (e * g) % p
    if e == h:
      return K + j, small_steps, 0
    else:
      lista[e] = j
      small_steps += 1
  u = pow(g, p - n - 1, p)
  s = pow(g, n, p)
  ep = em = h
  for i in range(0, n):
    try:
      j = lista[em]
      return j + K - i * n, small_steps, 2 * (i + 1)
    except KeyError:
      try:
        j = lista[ep]
        return K + j + i * n, small_steps, 2 * (i + 1)
      except KeyError:
        ep = (ep * u) % p
        em = (em * s) % p
  return None

def shanks_interleved(g, h, p, ordin):
  n = ceil(sqrt(ordin))
  e = 1
  u = pow(g, p - n - 1, p)
  s = h
  baby, giant = {e:0}, {s:0}
  # lista = {e:0, s:0}
  for j in range(1, n + 1):
    e = (e * g) % p
    if e == h:
      return j, 2 * j + 1, getsizeof(baby) + getsizeof(giant) #getsizeof(lista)
    else:
      i = giant.get(e, None)
      # i = lista.get(e, None)
      if i:
        return (j + i*n) % ordin, 2 * j + 1, getsizeof(baby) + getsizeof(giant) #getsizeof(lista)
      else:
        baby[e] = j
        # lista[e] = j
      s = (s * u) % p
      i = baby.get(s, None)
      # i = lista.get(s, None)
      if i:
        return (j*n + i) % ordin, 2 * (j + 1), getsizeof(baby) + getsizeof(giant)#getsizeof(lista)
      else:
        giant[s] = j
        # lista[s] = j
  return None

def shanks_two_grumpys_one_baby(g, h, p, ordin):
  n = ceil(sqrt(ordin))
  m = ceil(n/2)
  u = pow(g, m, p)
  s = pow(g, p - (m + 1) - 1, p)
  b = 1
  giant1 = h
  giant2 = pow(h, 2, p)
  babys, giants1, giants2 = {b: 0}, {giant1: 0}, {giant2: 0}
  for i in range(1, n + 1):
    b = (b * g) % p
    if b == h:
      return i, 3 * (i - 1) + 1, getsizeof(babys) + getsizeof(giants1) + getsizeof(giants2)
    else:
      j = giants1.get(b, None)
      if j:
        return (i - j * m) % ordin, 3 * (i - 1) + 1, getsizeof(babys) + getsizeof(giants1) + getsizeof(giants2)
      else:
        j = giants2.get(b, None)
        if j:
          r = (i + j * (m + 1)) // 2
          l = ordin // 2
          for _ in range(3):
            if pow(g, r, p) == h:
              return r, 3 * (i - 1) + 1, getsizeof(babys) + getsizeof(giants1) + getsizeof(giants2)
            r += l
        else:
          babys[b] = i

    giant1 = (giant1 * u) % p
    j = babys.get(giant1, None)
    if j:
      return (j - i * m) % ordin, 3 * (i - 1) + 2, getsizeof(babys) + getsizeof(giants1) + getsizeof(giants2)
    else:
      j = giants2.get(giant1, None)
      if j:
        return (i * m + j * (m + 1)) % ordin, 3 * (i - 1) + 2, getsizeof(babys) + getsizeof(giants1) + getsizeof(giants2)
      else:
        giants1[giant1] = i

    giant2 = (giant2 * s) % p
    j = babys.get(giant2, None)
    if j:
      r = (j + i * (m + 1)) // 2
      l = ordin // 2
      for _ in range(3):
        if pow(g, r, p) == h:
          return r, 3 * i , getsizeof(babys) + getsizeof(giants1) + getsizeof(giants2)
        r += l
    else:
      j = giants1.get(giant2, None)
      if j:
        return (j * m + i * (m + 1)) % ordin, 3 * i, getsizeof(babys) + getsizeof(giants1) + getsizeof(giants2)
      else:
        giants2[giant2] = i
  return None

def Shanks(g: int, h: int, p: int, ordin: int, type = 'classic', r = 2):
  if ordin:
    if type == 'classic':
      return shanks_classic(g, h, p, ordin)
    elif type == 'inter':
      return shanks_interleved(g, h, p, ordin)
    elif type == '2giants':
      return shanks_two_grumpys_one_baby(g, h, p, ordin)
    elif type == 'factor':
      l, m = optim_factors(ordin)
      return shanks_factor(g, h, p, l, m)
    elif type == 'centered':
      return shanks_centered(g, h, p, ordin)
    elif type == 'general':
      return shanks_general(g, h, p, ordin, r)
    elif type == 'binary':
      return shanks_classic_binary_search(g, h, p, ordin)
    elif type == 'betteravg':
      return shanks_classic_better_avg(g, h, p, ordin)
  return Shanks_order_unkown(g, h, p)





