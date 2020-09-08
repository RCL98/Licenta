
from random import randint, seed, sample
from ast import literal_eval as make_tuple
from aritmetica_modulara import sqrt, inverse, chinese_lemma, matmul, gauss_mod, avg, hensel_lifting
from sympy.ntheory import factorint
from math import exp, sqrt, log, ceil, log2
from generate_primes import logarithm_test_numbers, sieveOfEratosthenes, primes_upto
from timeit import default_timer as timer
import numpy as np
# from primefac import factorint
# from multiprocessing import Queue, cpu_count, Pool
from subprocess import Popen, PIPE
# from copy import deepcopy
import asyncio
from aiostream import stream

# def worker(congruences, s):
#     seed(s + 10)
#     j = 0
#     while True:
#         k = randint(2, p)
#         res = is_Bsmooth(B, pow(g, k, p))
#         if res[0]:
#             congruences.append((res[1], k))
#         if len(congruences) >= max_equations: break
#         # if j % 50 == 0:
#         #     print(len(congruences))
#         # j += 1

# L complexity notation
def L(x):
    x = 2 ** ceil(log2(x))
    return exp(sqrt(2)*sqrt(log(x) *log(log(x))))

def shanks_complex(x):
    return 2 ** ceil(log2(x)/2)

def primes_up_to_B(B: int):
    a = [1] * B
    for i in range(2, ceil(sqrt(B))):
        if a[i]:
            for j in range(i * i, B, i): a[j] = 0
    return [i for i in range(len(a))if a[i]==1][2:]

# determine if n is B-smooth (uses fast factoring)
def is_Bsmooth(B: int, primes_B: list, n: int):
    factorization = {}
    for prime in primes_B:
        e = 0
        while n % prime == 0:
            e += 1
            n = n // prime
        factorization[prime] = e
        if n == 1:
            break
    if n != 1:
        return False, factorization
    return True, factorization

async def worker(g, inv_h, p, B, primes_B, rand_seed):
    seed(rand_seed)
    while 1:
        k = randint(1, p)
        factorization = {}
        n = (inv_h * pow(g, k, p)) % p
        for prime in primes_B:
            e = 0
            while n % prime == 0:
                e += 1
                n = n // prime
            factorization[prime] = e
            if n == 1:
                break
        if n == 1:
            yield factorization, k

async def worker_M(g, inv_h, p, B, primes_B, q, rand_seed):
    ok = True
    seed(rand_seed)
    while ok:
        k = randint(1, p)
        factorization = {}
        n = (inv_h * pow(g, k, p)) % p
        for prime in primes_B:
            e = 0
            while n % prime == 0:
                e += 1
                n = n // prime
            factorization[prime] = e
            if n == 1:
                break
        if n == 1:
            try:
                q.put_nowait((factorization, k))
            except asyncio.QueueFull:
                ok = False
            #


# find congruences in accordance with Hoffstein, Phipher, Silverman (3.32)
async def find_congruences_lect11_multyQueue(g, h, p, B: int, primes: list, max_equations: int):
    congruences, bases = [], []
    unique = lambda l: list(set(l))
    inv_h = inverse(h, p)
    new_seeds = sample(range(1, 5000), max_equations)
    q = asyncio.Queue(max_equations)
    tasks = []
    for s in new_seeds:
        task = asyncio.create_task(worker_M(g, inv_h, p, B, primes, q, s))
        tasks.append(task)
    await asyncio.gather(*tasks)
    # await q.join()
    while True:
        try:
            item = q.get_nowait()
            congruences.append((item[0], item[1]))
        except asyncio.QueueEmpty:
            break
    bases = unique([base for c in [c[0].keys() for c in congruences] for base in c])
    return bases, congruences

# find congruences in accordance with Hoffstein, Phipher, Silverman (3.32)
async def find_congruences_lect11_paralel(g, h, p, B: int, primes: list, max_equations: int):
    congruences, bases = [], []
    unique = lambda l: list(set(l))
    inv_h = inverse(h, p)
    # while True:
    new_seeds = sample(range(1, 5000), 4)
    #[worker(g, h, p, B, primesB, s) for s in new_seeds]
    combined_stream = stream.merge(worker(g, inv_h, p, B, primes, new_seeds[0]), worker(g, inv_h, p, B, primes, new_seeds[1]),
                                   worker(g, inv_h, p, B, primes, new_seeds[2]), worker(g, inv_h, p, B, primes, new_seeds[3]))
    async with combined_stream.stream() as streamer:
        async for item in streamer:
            congruences.append((item[0], item[1]))
            bases = unique([base for c in [c[0].keys() for c in congruences] for base in c])
            if len(congruences) >= max_equations: break
    return bases, congruences

# find congruences in accordance with Hoffstein, Phipher, Silverman (3.32)
def find_congruences_lect11(g, h, p, B: int, primes: list, max_equations: int):
    congruences, bases = [], []
    unique = lambda l: list(set(l))
    inv_h = inverse(h, p)
    while True:
        k = randint(2, p)
        _ = is_Bsmooth(B, primes, (inv_h * pow(g, k, p)) % p )
        if _[0]:
            congruences.append((_[1], k))
            if len(congruences) >= max_equations: break
    bases = unique([base for c in [c[0].keys() for c in congruences] for base in c])
    return bases, congruences

# find congruences in accordance with Hoffstein, Phipher, Silverman (3.32)
def find_congruences(g, h, p, B: int, primes: list, max_equations: int):
    congruences, bases = [], []
    unique = lambda l: list(set(l))
    while True:
        k = randint(2, p)
        _ = is_Bsmooth(B, primes, (h * pow(g, k, p)) % p )
        if _[0]:
            congruences.append((_[1], k))
            if len(congruences) >= max_equations: break
    bases = unique([base for c in [c[0].keys() for c in congruences] for base in c])
    return bases, congruences

# convert the linear system  to dense matrices 
def to_matrices(bases, congruences):
    M = [[c[0][base] if base in c[0].keys() else 0 \
            for base in bases] for c in congruences]
    b = [c[1] for c in congruences]
    return M, b

# convert the linear system  to dense matrices
def to_matrices_lect11(bases, congruences):
    M = [[c[0][base] if base in c[0].keys() else 0 \
            for base in bases]  + [1] for c in congruences]
    b = [c[1] for c in congruences]
    return M, b

# use sage to solve (potentially) big systems of equations:
def msolve(M, b, p):
    # sage_cmd = 'L1 = {};L2 = {};R = IntegerModRing({});M = Matrix(R, L1); b = vector(R, L2);print(M.solve_right(b))'.format(M,b,p-1)
    # with open('E:\\Sage\\run.sage', 'w') as output_file: output_file.write(sage_cmd)
    # P = Popen('test.bat', stdout=PIPE, stderr=PIPE)
    # stdout, stderr = P.communicate()
    # stdout = str(stdout).split('\\n')
    # cmd_result = make_tuple(stdout[-2])
    # # if cmd_result[0] == 0: print('sage failed with error {}'.format(cmd_result)); exit()
    # return cmd_result
    factors = factorint(p - 1)
    modules = list(factors.keys())
    mod_solutions, logs = [], [0] * len(M[0])
    for mod in modules:
        M_mod = np.array(M) % mod
        b_mod = np.array(b) % mod
        results = gauss_mod(M_mod, b_mod, mod)
        if results:
            if factors[mod] >= 2:
                results = hensel_lifting(results, M, b, factors[mod], mod) # hensel lifting
            # M_m = np.array(M) % mod
            # b_b = np.array(b) % mod
            # x = np.array(results)
            # if np.array_equal((M_m @ x) % mod, b_b):
            #     print("DA")
            # else:
            #     C = (M_m @ x) % mod
            mod_solutions.append(results)
        else: return None
    modules = [mod ** factors[mod] for mod in modules]
    for i in range(len(M[0])):
        rests = [int(xs[i]) for xs in mod_solutions]
        logs[i] = chinese_lemma(modules, rests)
    return logs

# solves a linear equation 
def evaluate(eq, dlogs, k, p):
    return (sum([dlogs[term] * exp for term, exp in eq.items()]) + k) % (p-1)

def check_congruences(congruences, dlogs):
    # print('checking congruences:');
    passed = True
    for c in congruences:
        if evaluate(c[0], dlogs) != c[1]: passed = False
    if passed: o = 0 #print('Passed!\n')
    else:
        #print('Failed, try running again?');
        exit()
    return passed

def check_dlogs(exponents, bases):
    #print('checking dlog exponents:');
    passed = True
    for exponent, base in zip(exponents, bases):
        if pow(g, exponent, p) != base: passed = False
        # else: print('{}^{} = {} (mod {})'.format(g,exponent, base, p))
    if passed: o = 0 #print('Passed!\n')
    else:
        #print('Failed, try running again.');
        exit()
    return passed

async def search_for_K(g, h, p, B, primes_B):
    for i in range(10**9):
        k = randint(1, p)
        factorization = {}
        n = (h * pow(g, p - k - 1, p)) % p
        for prime in primes_B:
            e = 0
            while n % prime == 0:
                e += 1
                n = n // prime
            factorization[prime] = e
            if n == 1:
                break
        if n == 1:
            return k, factorization
    return None

def SELD(g, h, p, B, max_equations, primesB, rand_seed):
    # print('searching for congruences.')
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    k, factorization = loop.run_until_complete(search_for_K(g, h, p, B, primesB))
    exponents = None
    while exponents == None:
        seed(rand_seed)
        bases, congruences = find_congruences(g, h, p, B, primesB, max_equations)
        M, b = to_matrices(bases, congruences)
        exponents = msolve(M, b, p)
        M_m = np.array(M)
        b_b = np.array(b)
        x = np.array(exponents)
        if np.array_equal(M_m @ x, b_b):
            print("DA")
        else:
            C = M_m @ x
        rand_seed = randint(1, 5000)
    # dictionary of bases and exponents
    dlogs = {b: int(exp) for (b,exp) in zip(bases, exponents)}


    soln = evaluate(factorization, dlogs, k, p)
    if pow(g, soln, p) == h:
        return soln
    # # print('searching for k such that h*g^-k is B-smooth.')
    # for i in range(10**9):
    #     k = randint(2, p)
    #     c = is_Bsmooth(B, primesB, (h * pow(g, p - k - 1, p)) % p)
    #     if c[0]:
    #         soln = evaluate(c[1], dlogs, k, p)
    #         if pow(g, soln, p) == h:
    #             return soln
    #         else:
    #             print('Failed exp.', p)
    return SELD(g, h, p, B, max_equation, primesB, randint(1, 500))

def SELD_lect11(g, h, p, B, max_equations, primesB, rand_seed):
    exponents = None
    it = 0
    t_cong = 0
    while exponents == None:
        seed(rand_seed)
        # t = timer()
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)
        bases, congruences = loop.run_until_complete(find_congruences_lect11_paralel(g, h, p, B, primesB, max_equations))
        # bases, congruences = find_congruences_lect11(g, h, p, B, primesB, max_equations)
        # bases, congruences = find_congruences_lect11_multyprocess(g, h, p, B, primesB, max_equations)
        # t_cong = (timer() - t) * 1000
        # print(f"Find congurences time: {(timer() - t) * 1000:.3f}")
        # t = timer()
        M, b = to_matrices_lect11(bases, congruences)
        exponents = msolve(M, b, p)
        # print(f"Solve linear system: {(timer() - t) * 1000:.3f}")
        rand_seed = randint(1, 5000)
        # print(f"Try: {it}\n")
        it += 1
        if exponents != None and pow(g, exponents[-1], p) != h:
            # M = np.array(M)
            # b = np.array(b)
            # exponents = np.array(exponents)
            # if not np.array_equal(M @ exponents, b):
            #     print("Failed System")
            exponents = None
    # dictionary of bases and exponents
    return exponents[-1] % (p - 1), t_cong, it

if __name__ == '__main__':
    i = 0
    ts, tts = [], []
    its = []
    seed(42)
    for q in range(5):
        for _ in range(500):
            ok = 1
            while ok:
                ok = 0
                p, g, exponent, h = logarithm_test_numbers(15, safe = False)
            x = 2 ** ceil(log2(p - 1))
            B = ceil(exp(1/sqrt(2) * sqrt(log(x) * log(log(x))))/2)
            # print(L(p), shanks_complex(p))
            primesB = primes_up_to_B(B)
            max_equations = primes_upto(B) + 5
            # print(B, max_equations, p)
            t = timer()
            x, ta, it = SELD_lect11(g, h, p, B, max_equations, primesB, 42)
            ts.append((timer() - t) * 1000)
            its.append(it)
            if x == exponent:
                # print("DA")
                i += 1
            else:
                if pow(g, x, p) != h:
                    print("Nu")
        avg_ts = avg(ts)
        tts.append(avg_ts)
        print("\n", i, avg_ts)
        print("\n", i, sum(its))
    print(avg(tts))