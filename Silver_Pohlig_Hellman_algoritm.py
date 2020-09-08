from Shanks_algoritm import Shanks
from Pollard_algoritmi import pollard_rho
from aritmetica_modulara import sqrt, inverse, chinese_lemma, countBits
from sympy.ntheory import factorint

def pohlig_hellman_alg(g, h, p, q, e, rads, solver = Shanks):
  if g == h:
    return 1
  if h == 1:
    return 0
  if e == 1:
    if countBits(p) <= 10:
      return rads[h]
    return solver(g, h, p, q)[0]

  if countBits(p) <= 10:
    f = lambda h, g: rads[h]
  else:
    f = lambda h_i, g_i: solver(g_i, h_i, p, q)[0]

  z = h; u = inverse(g, p); t = pow(q, e - 1); r = pow(g, q ** (e - 1), p)
  b = [0] * e
  for i in range(e):
    w = pow(h, t, p)
    b[i] = f(w, r)
    h = (h * pow(u, b[i], p)) % p
    u = pow(u, q, p)
    t = t//q
  f = lambda i: b[i] * (q ** i)
  return sum(f(i) for i in range(e)) % p

def shoup_alg(g, h, p, q, e, rads, solver = Shanks):
  if g == h:
    return 1
  if e == 1:
    if countBits(p) <= 10:
      return rads[h]
    return solver(g, h, p, q)[0]

  f = e // 2;
  u = shoup_alg(pow(g, q ** (e - f), p), pow(h, q ** (e - f), p), p, q, f, rads, solver)
  v = shoup_alg(pow(g, q ** f, p), (h * pow(g, (q ** e) - u, p)) % p, p, q, e - f, rads, solver)
  return ((q ** f)*v + u) % (p - 1)

def pld_rezolv_qe(g, h, p, q, e, solver = Shanks):
  if g == h:
    return 1
  if e == 1:
    return solver(g, h, p, q)[0]

  x = 0
  g_e = pow(g, q ** (e - 1), p)
  u = inverse(g, p)
  for e_i in reversed(range(e)):
    u = pow(g, p - 1 - x, p)
    h_e = pow(h * u, q ** e_i, p)
    x += solver(g_e, h_e, p, q)[0] * (q ** (e - 1 - e_i))
  return x

def silver_pohlig_hellman_alg(g, h, p, fact = factorint, qeAlg = pld_rezolv_qe, solver = Shanks):
  N = p - 1
  factorizare_prima = fact(N)
  y = [0] * len(factorizare_prima)
  for ind, (factor, mult) in enumerate(factorizare_prima.items()):
    l = N // (factor ** mult)
    g_i = pow(g, l, p)
    h_i = pow(h, l, p)
    T = {1: 0}
    if countBits(factor) <= 10:
      for i in range(1, factor):
        # print(pow(g_i, i * (factor ** (mult - 1)), p))
        T[pow(g_i, i * (factor ** (mult - 1)), p)] = i
    y[ind] = qeAlg(g_i, h_i, p, factor, mult, T, solver = solver)
    j = 1
  module = [pow(factor, mult, p) for (factor, mult) in factorizare_prima.items()]
  return chinese_lemma(module, y)


