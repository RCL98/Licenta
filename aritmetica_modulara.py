import random
from functools import reduce
from sympy.ntheory import factorint
import math

def avg(lst):
  return reduce(lambda a, b: a + b, lst) / len(lst)

def tonelli_shanks(n, p):
  """Tonelli-shanks algorithm for computing the square root
  of n modulo a prime p.

  n must be in the range [0..p-1].
  p must be at least even.

  The return value r is the square root of modulo p. If non-zero,
  another solution will also exist (p-r).

  Note we cannot assume that p is really a prime: if it's not,
  we can either raise an exception or return the correct value.
  """

  # See https://rosettacode.org/wiki/Tonelli-Shanks_algorithm

  if n in (0, 1):
    return n

  if p % 4 == 3:
    root = pow(n, (p + 1) // 4, p)
    if pow(root, 2, p) != n:
      raise ValueError("Cannot compute square root")
    return root

  s = 1
  q = (p - 1) // 2
  while not (q & 1):
    s += 1
    q >>= 1

  z = n.__class__(2)
  while True:
    euler = pow(z, (p - 1) // 2, p)
    if euler == 1:
      z += 1
      continue
    if euler == p - 1:
      break
    # Most probably p is not a prime
    raise ValueError("Cannot compute square root")

  m = s
  c = pow(z, q, p)
  t = pow(n, q, p)
  r = pow(n, (q + 1) // 2, p)

  while t != 1:
    for i in range(0, m):
      if pow(t, 2 ** i, p) == 1:
        break
    if i == m:
      raise ValueError("Cannot compute square root of %d mod %d" % (n, p))
    b = pow(c, 2 ** (m - i - 1), p)
    m = i
    c = b ** 2 % p
    t = (t * b ** 2) % p
    r = (r * b) % p

  if pow(r, 2, p) != n:
    raise ValueError("Cannot compute square root")

  return r

def sqrt(value, modulus=None):
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
    result = tonelli_shanks(value % modulus, modulus)

  return result


def inverse(value, modulus):
  if modulus == 0:
    raise ZeroDivisionError("Modulus cannot be zero")
  if modulus < 0:
    raise ValueError("Modulus cannot be negative")
  r_p, r_n = value, modulus
  s_p, s_n = 1, 0
  while r_n > 0:
    q = r_p // r_n
    r_p, r_n = r_n, r_p - q * r_n
    s_p, s_n = s_n, s_p - q * s_n
  if r_p != 1:
    raise ValueError("No inverse value can be computed" + str(r_p))
  while s_p < 0:
    s_p += modulus
  return s_p

def lema_chineza_rez(module, resturi):
  N = reduce(lambda a, b: a * b, module)
  rezultat = 0
  for (n_i, x_i) in zip(module, resturi):
    N_i = N // n_i
    rezultat += x_i * inverse(N_i, n_i) * N_i
  return rezultat % N

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

# random.seed(33)

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

def countBits(number):
  # log function in base 2
  # take only integer part
  return int((math.log(number) /
              math.log(2)) + 1);

# # print(extended_euclid(3233, 560))
# print(sqrt(29))