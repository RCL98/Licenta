# import numpy as np
import pandas as pd
import random as rand
import math
import aritmetica_modulara as artmod
from sympy.ntheory import factorint

def genereaza_numere(numBits: int, cardinal: int) -> list:
  numere = []
  for ind in range(cardinal):
    numere.append(rand.getrandbits(numBits))
    #binNumbers.append(np.random.binomial(size = numOfBits, n = 1, p = 0.5))
  return numere

def obtine_parametrii(numar: int) -> (int, int):
  numar -= 1
  q = 0
  while (numar & 1) == 0:
    q += 1
    numar >>= 1
  return q, numar

def alege_martori(numar: int, cardinal: int) -> list:
  martori = set()
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
    numar = rand.getrandbits(numBits)
    if miller_rabin_test(numar, mr_iteratii) == True:
      numere_prime.append(numar)
  return numere_prime

def pohlig_primes(magnitude):
  m = 2
  N = 1
  while math.log10(N) < magnitude or not miller_rabin_test(N + 1, 15):
    N *= m
    m = m + 1
  p = N + 1
  print(p)
  print("p =", m - 1, "! + 1 is prime")
  return p

# def logarithm_test_numbers(size: int, numOfBts: int, to_csv = False) ->  pd.DataFrame:
#   primes = get_primes(numOfBts, size)
#   generators = [artmod.generator_Zp(p) for p in primes]
#   exponents = [rand.randint(2, p - 1) for p in primes]
#   h_values = [pow(g, e, p) for (g, e, p) in zip(generators, exponents, primes)]
#   if to_csv:
#     plds = pd.DataFrame(list(zip(primes, generators, exponents, h_values)), columns = ['Primes', 'Generators', 'Exponents', 'h_values'])
#     plds.to_csv('plds.csv', index=False, header=True, sep='\t', encoding='utf-8')
#   return primes, generators, exponents, h_values

def logarithm_test_numbers(numOfBts: int, centered = False):
  prime = get_primes(numOfBts, 1)[0]
  generator = artmod.generator_Zp(prime)
  if centered:
    exp = rand.randint((prime - 1) // 2 - (prime - 1) // 5,
                       (prime - 1) // 2 + (prime - 1) // 5 + 1)
  else:
    exp = rand.randint(2, prime - 1)
  h = pow(generator, exp, prime)
  return prime, generator, exp, h

def logarithm_test_numbers_same_p(numOfValues: int, numOfBts: int):
  prime = get_primes(numOfBts, 1)[0]
  generator = artmod.generator_Zp(prime)
  exps = rand.sample(range(2, prime - 1), numOfValues)
  hs = [pow(generator, exp, prime) for exp in exps]
  return prime, generator, exps, hs

# plds = logarithm_test_numbers(35, 40)



# rand.seed(45)
# # q = pohlig_primes(320)
# prime = numere_prime_posibile(24, 5)
# # print(prime[0] * prime[1], prime[0:2])
# prime_gen = {}
# for prim in prime:
#   prime_gen[prim] = artmod.generator_Zp(prim)
#
# for (prim, gen) in prime_gen.items():
#   print(f"Nr prim: {prim} - Generator: {gen}")


