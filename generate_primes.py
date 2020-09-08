import random as rand
import math
import aritmetica_modulara as artmod
from sympy.ntheory import factorint

def obtine_parametrii(numar: int) -> (int, int):
  numar -= 1
  q = 0
  while (numar & 1) == 0:
    q += 1
    numar >>= 1
  return q, numar

def alege_martori(numar: int, cardinal: int) -> list:
  martori = set()
  if numar//cardinal < 2:
    cardinal = numar // 2
  while len(martori) < cardinal:
    martori.add(rand.randint(2, numar - 2))
  return martori

def miller_rabin_test(numar: int, iteratii: int) -> bool:
  if numar <= 3:
    return True
  elif (numar & 1) == 0:
    return False

  q, d = obtine_parametrii(numar)
  # martori = rand.sample(range(2, numar - 2), iteratii)
  martori = alege_martori(numar, iteratii)
  for a in martori:
    g = math.gcd(a, numar)
    if g > 1 and g < numar:
      return False

    a = pow(a, d, numar)
    if a == 1:
      continue
    for __ in range(q):
      if a == numar - 1:
        break
      a = pow(a, 2, numar)
    else:
      return False
  return True

def get_primes(numBits: int, cardinal: int) -> list:
  mr_ranges = ((220, 30), (280, 20), (390, 15), (512, 10),
               (620, 7), (740, 6), (890, 5), (1200, 4),
               (1700, 3), (3700, 2))
  try:
    mr_iteratii = list(filter(lambda x: numBits < x[0],
                                mr_ranges))[0][1]
  except IndexError:
    mr_iteratii = 2

  numere_prime = []
  while len(numere_prime) < cardinal:
    numar = rand.randint(2 ** (numBits - 1), 2 ** numBits - 1) #rand.getrandbits(numBits)
    if miller_rabin_test(numar, mr_iteratii) == True:
      numere_prime.append(numar)
  return numere_prime

def get_prime(numOfBits: int, safe = False):
  if safe == False:
    return get_primes(numOfBits, 1)[0]
  else:
    iter = 0
    while 1 or iter < 10**10:
      p = get_primes(numOfBits, 1)[0]
      t = factorint(p - 1)
      if len(t) == 2:
        if 2 in t.keys() and t[2] == 1:
          return p
      iter += 1
    print("No safe prime found in max iterations!")

def primes_upto(limit):
  is_prime = [False] * 2 + [True] * (limit - 1)
  for n in range(int(limit ** 0.5 + 1.5)):  # stop at ``sqrt(limit)``
    if is_prime[n]:
      for i in range(n * n, limit + 1, n):
        is_prime[i] = False
  return len([i for i, prime in enumerate(is_prime) if prime])

def sieveOfEratosthenes(n):
  # Create a boolean array "prime[0..n]" and initialize
  #  all entries it as true. A value in prime[i] will
  # finally be false if i is Not a prime, else true.
  prime = [True for i in range(n)]
  p = 2
  while (p * p <= n):
    # If prime[p] is not changed, then it is a prime
    if (prime[p] == True):
      # Update all multiples of p
      for i in range(p * p, n, p):
        prime[i] = False
    p += 1
  return sum(prime)

def pohlig_primes(magnitude):
  m = 2
  N = 1
  while math.log2(N) < magnitude or not miller_rabin_test(N + 1, 15):
    N *= m
    m = m + 1
  p = N + 1
  # print(p)
  # print("p =", m - 1, "! + 1 is prime")
  return p

def logarithm_test_numbers(numOfBits: int, centered = False, silvPohHell = False, safe = False):
  if not silvPohHell:
    prime = get_prime(numOfBits, safe)
    while prime < 10:
      prime = get_prime(numOfBits, safe)
  else:
    prime = pohlig_primes(numOfBits)
    while prime < 10:
      prime = get_prime(numOfBits, safe)

  generator = artmod.generator_Zp(prime)
  if centered:
    exp = rand.randint((prime - 2) // 2 - (prime - 2) // 5,
                       (prime - 2) // 2 + (prime - 2) // 5 + 1)
  else:
    exp = rand.randint(2, prime - 2)
  h = pow(generator, exp, prime)
  return prime, generator, exp, h

def logarithm_test_numbers_same_p(numOfValues: int, numOfBts: int):
  prime = get_primes(numOfBts, 1)[0]
  generator = artmod.generator_Zp(prime)
  exps = rand.sample(range(2, prime - 2), numOfValues)
  hs = [pow(generator, exp, prime) for exp in exps]
  return prime, generator, exps, hs

def worst_shanks_numbers(numOfBits: int):
  prime = get_primes(numOfBits, 1)[0]
  while prime < 10:
    prime = get_primes(numOfBits, 1)[0]
  generator = artmod.generator_Zp(prime)
  n = math.ceil(artmod.sqrt(prime - 1))
  exp =  n * (n - 1)
  if exp >= prime - 1:
    exp = n * (n - 2)
  h = pow(generator, exp, prime)
  return prime, generator, exp, h
