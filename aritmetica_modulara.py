import random
from functools import reduce
from sympy.ntheory import factorint
import math
from copy import deepcopy
import numpy as np

def avg(lst):
  return reduce(lambda a, b: a + b, lst) / len(lst)

def countBits(number):
  return int((math.log(number) /
              math.log(2)) + 1);

def chinese_lemma(module, resturi):
  N = reduce(lambda a, b: a * b, module)
  rezultat = 0
  for (n_i, x_i) in zip(module, resturi):
    N_i = N // n_i
    rezultat += x_i * inverse(N_i, n_i) * N_i
  return rezultat % N

def sqrt(value, modulus=None, prime = True):
  if modulus is None:
    if value < 0:
      raise ValueError("Square root of negative value")
    # http://stackoverflow.com/questions/15390807/integer-square-root-in-python

    x = value
    y = (x + 1) / 2
    while y < x:
      x = y
      y = (x + value / x) / 2
    result = x
  else:
    if modulus <= 0:
      raise ValueError("Modulus must be positive")
    if prime:
      result = tonelli_shanks(value % modulus, modulus)
    else:
      fact = factorint(modulus)
      moduls, rests = [], []
      for mod, exp in fact.items():
        rest = tonelli_shanks(value % mod, mod)
        moduls += [mod] * exp
        rests += [rest] * exp
      result = chinese_lemma(moduls, rests)
  return result

def inverse(value, modulus):
  if modulus == 0:
    raise ZeroDivisionError("Modulus cannot be zero")
  if modulus < 0:
    raise ValueError("Modulus cannot be negative")
  value = value % modulus
  if value == 0:
    return 0
  r_p, r_n = value, modulus
  s_p, s_n = 1, 0
  while r_n > 0:
    q = r_p // r_n
    r_p, r_n = r_n, r_p - q * r_n
    s_p, s_n = s_n, s_p - q * s_n
  if r_p != 1:
    raise ValueError("No inverse value can be computed " + str(r_p))
  while s_p < 0:
    s_p += modulus
  return s_p

def ordinul(numar, modul):
  factorizare = factorint(modul - 1)
  N = modul - 1
  for (factor, mult) in factorizare.items():
    N = N // (factor ** mult)
    a = pow(numar, N, modul)
    while a != 1:
      a = pow(a, factor, modul)
      N = N * factor
  return N

def generator_Zp(p):
  factorizare = factorint(p - 1)
  N = p - 1
  gasit = False
  while not gasit:
    a = random.randint(1, N)
    for factor in factorizare.keys():
      b = pow(a,  N // factor, p)
      if b == 1:
        break
    else:
      gasit = True
  return a

def extended_euclid(a, b):
  xprec, x = 1, 0
  yprec, y = 0, 1
  while b:
    q = a // b
    x, xprec = xprec - q * x, x
    y, yprec = yprec - q * y, y
    a, b = b, a % b
  return a, xprec, yprec

def gauss_mod(A: np.ndarray, b: np.ndarray, modulus: int) -> list:
  n = len(A)
  # for row, el in zip(A, b):
  #   row.append(el)
  A = np.c_[A, b]
  m = len(A[0])
  nr_sol = m - 1

  A = A.tolist()
  b = b.tolist()
  for i in range(0, nr_sol):
    # Search for maximum in this column
    maxEl = abs(A[i][i])
    maxRow = i #if maxEl < modulus - 1:
    for k in range(i + 1, n):
      if abs(A[k][i]) > maxEl:
        maxEl = abs(A[k][i])
        maxRow = k
        if maxEl == modulus - 1:
          break

    # Swap maximum row with current row (column by column)
    if maxRow != i:
      for k in range(i, m):
        # tmp = A[maxRow][k]
        A[maxRow][k], A[i][k]  = A[i][k], A[maxRow][k]
        # A[i][k] = tmp

    # Make all rows below this one 0 in current column
    inv_pivot = inverse(A[i][i], modulus)
    for k in range(i + 1, n):
      c = (-A[k][i] * inv_pivot) % modulus
      if c == 0:
        continue
      for j in range(i, m):
        if i == j:
          A[k][j] = 0
        else:
          A[k][j] = (A[k][j] + c * A[i][j]) % modulus

  # Solve equation Ax=b for an upper triangular matrix A

  x = [0 for i in range(nr_sol)]
  for i in range(nr_sol - 1, -1, -1):
    if A[i][i] != 0:
      x[i] = (A[i][nr_sol] * inverse(A[i][i], modulus)) % modulus
    else:
      return None
    for k in range(i - 1, -1, -1):
      A[k][nr_sol] = (A[k][nr_sol] - A[k][i] * x[i]) % modulus
  return x

def hensel_lifting(solutions, M, b, power, modulus):
  new_solutios = [0] * len(solutions)
  A = np.array(M)
  y = np.array(b)
  x = np.array(solutions)
  for i in range(1, power):
    if modulus == 3:
      z = (y - A @ x)/ modulus ** i
    y_2 = ((y - A @ x)//modulus ** i) % modulus
    A_mod = A % modulus
    x += np.array(gauss_mod(A_mod, y_2, modulus)) * modulus ** i
  return x.tolist()


