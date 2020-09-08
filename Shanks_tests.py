import generate_primes as gnp
import numpy as np
from timeit import default_timer as timer
from  math import ceil, log2, log
from aritmetica_modulara import sqrt, avg, generator_Zp
from random import randint, seed, sample
from Shanks_algoritm import Shanks

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

def r_optim(n, m):
  return 2*log(n)/log(n/m)

def test_shanks_same_p(numOfTests: int, numOfBits: int, rs = [2], radn_seed = 0):
  seed(radn_seed)
  times_g, smalls_g, giants_g, memsizes_g = [], [], [], []
  p, g, exps, hs = gnp.logarithm_test_numbers_same_p(numOfTests, numOfBits)
  l, m = optim_factors(p - 1)
  r_opt = r_optim(p-1,numOfTests)
  print(f"r optim:{r_opt:.2f}")

  print("General started!")
  t_start_c = timer()
  x_cs, small_steps_c, giant_steps_c, memsize_c = shanks_classic_with_memory(g, hs, p, p - 1)
  time_c = (timer() - t_start_c) * 1000
  print("General finished!")

  print("Factor started!")
  t_start_f = timer()
  x_fs, small_steps_f, giant_steps_f, memsize_f = shanks_factor_with_memory(g, hs, p, l, m)
  time_f = (timer() - t_start_f) * 1000
  print("Factor finished!")

  rs += [r_opt]
  for r in rs:
    print(f"Gen {r:.3f} started!")
    t_start_g = timer()
    x_gs, small_steps_g, giant_steps_g, memsize_g = shanks_general_with_memory(g, hs, p, p - 1, r)
    times_g.append((timer() - t_start_g)*1000)
    smalls_g.append(small_steps_g)
    giants_g.append(avg(giant_steps_g))
    memsizes_g.append(memsize_g)
    print(f"Gen {r:.3f} finished!")

  print(f"t_c: {time_c/numOfTests:.3f} ms, t_f: {time_f/numOfTests:.3f} ms, t_g:",end=' ')
  for r, t in zip(rs, times_g):
    print(f"r= {r}: {t/numOfTests:0.3f} ms",end=" ")
  print(f"\n\nsm_c: {small_steps_c/ceil(sqrt(p)):.3f}, sm_f: {small_steps_f/ceil(sqrt(p)):.3f}, sm_g:",end=' ')
  for r, s in zip(rs, smalls_g):
    print(f"r= {r}: {s/ceil(sqrt(p)):.3f}",end=" ")
  print(f"\n\ngc_avg: {avg(giant_steps_c)/ceil(sqrt(p)):.3f}, gf_avg: {avg(giant_steps_f)/ceil(sqrt(p)):.3f}, gg_avg:",end=' ')
  for r, g in zip(rs, giants_g):
    print(f"r= {r}: {g/ceil(sqrt(p)):.3f}",end=" ")
  print(f"\n\nsize_c: {memsize_c/ceil(sqrt(p)):.3f} bytes, size_f: {memsize_f/ceil(sqrt(p)):.3f} bytes, size_g:",end=' ')
  for r, m in zip(rs, memsizes_g):
    print(f"r= {r}: {m/ceil(sqrt(p)):.3f} bytes",end=" ")

def test_shanks_general(testStart=100, testEnd=900, testIter=100, bitStart=10, bitEnd=30, rs=[2], rand_seed=42):
  import matplotlib.pyplot as plt
  from matplotlib.lines import Line2D
  from itertools import cycle
  # seed(412)
  # marker = cycle(sample(Line2D.markers.keys(), len(rs) + 1))
  plt.tight_layout(pad=0.05)
  marker = cycle(['*', 'o', 's', 'p', '+'])

  seed(rand_seed)
  times = [[[0 for k in range(testStart, testEnd + testIter, testIter)] for j in range(len(rs) + 1)] for i in
           range(bitStart, bitEnd + 5, 5)]
  total_tests = list(range(testStart, testEnd + testIter, testIter));fig_ind = 0
  for bit_ind, bits in enumerate(range(bitStart, bitEnd + 10, 10)):
    print(f"Started: {bits}")
    p = gnp.get_primes(bits, 1)[0]
    g = generator_Zp(p)
    fig = plt.figure(fig_ind);fig_ind +=1
    # fig.suptitle(f'p = {bits} biti')
    r_optims = []
    for test_ind, tests in enumerate(total_tests):
      print(f"Curent Test: {tests}", end=' ')
      exps = sample(range(2, p - 1), tests)
      hs = [pow(g, exp, p) for exp in exps]
      r_opt = r_optim(p - 1, tests)
      r_optims.append(r_opt)
      for r_ind, r in enumerate(rs + [r_opt]):
        t_start = timer()
        res = shanks_general_with_memory(g, hs, p, p - 1, r)
        times[bit_ind][r_ind][test_ind] = (timer() - t_start) * 1000
      print(f"Time: {sum([sum(t) for t in times[bit_ind]]):.2e}")
    for ind in range(len(rs)):
      plt.plot(total_tests, times[bit_ind][ind], marker=next(marker), label=f"r:{rs[ind]}")
    plt.plot(total_tests, times[bit_ind][len(rs)], marker=next(marker), label="r optim")
    plt.legend()
    plt.ylabel("Milisecunde")
    plt.xlabel("Nr. de logaritmi")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.savefig(f'imagini/shanks_general_p_{bits}.png', bbox_inches='tight')
    fig = plt.figure(fig_ind);fig_ind += 1
    plt.plot(total_tests, r_optims, marker='X')
    plt.ylabel("r optim")
    plt.xlabel("Nr. de logaritmi")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
    plt.savefig(f'imagini/shanks_general_p_{bits}_r_optim.png', bbox_inches='tight')
    print(f"Finished: {bits}")
  plt.show()

def test_shanks_middle(numOfTests: int, numOfBits: list, randSeed = 0, centered = True):
  import matplotlib.pyplot as plt
  seed(randSeed)
  giant_steps_class_list, times_clasic = [], []
  giant_steps_centered_list, times_centered = [], []
  for bits in numOfBits:
    time_clas, time_cent = [], []
    giant_clas, giant_cent = [], []
    for iter in range(numOfTests):
      p, g, exp, h = gnp.logarithm_test_numbers(bits, centered)
      # print("iter:", iter)
      t = timer()
      x_c, _, giant_steps_c, _ = shanks_classic(g, h, p, p - 1)
      giant_clas.append(giant_steps_c)
      passed_time = (timer() - t)*1000
      time_clas.append(passed_time)
      # print(f"C_passed x:{x_c}, exp:{exp}, sm:{small_steps_c}, gs:{giant_steps_c} time:{passed_time:.3f}")

      t = timer()
      x_m, _, giant_steps_m = shanks_centered(g, h, p, p - 1)
      giant_cent.append(giant_steps_m)
      passed_time = (timer() - t) * 1000
      time_cent.append(passed_time)
      # print(f"M_passed x:{x_c}, exp:{exp}, sm:{small_steps_m}, gs:{giant_steps_m} time:{passed_time:.3f} \n")
    print()
    giant_steps_class_list.append(avg(giant_clas))
    giant_steps_centered_list.append(avg(giant_cent))
    times_clasic.append(avg(time_clas))
    times_centered.append(avg(time_cent))

  # fig = plt.figure(0)
  # plt.plot(numOfBits, times_clasic, color='b', marker = 'x', label = 'clasic')
  # plt.plot(numOfBits, times_centered, color='r', marker = 'o', label = 'centrat')
  # plt.ylabel("milisecunde")
  # plt.xlabel("biti")
  # plt.legend()
  # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  # plt.savefig(f'imagini/shanks_centered_timp.png', bbox_inches='tight')
  #
  # fig_1 = plt.figure(1)
  # plt.plot(numOfBits, giant_steps_class_list, color='b', marker = 'x', label = 'clasic')
  # plt.plot(numOfBits, giant_steps_centered_list, color='r', marker = 'o', label = 'centrat')
  # plt.ylabel("OM")
  # plt.xlabel("biti")
  # plt.legend()
  # # plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
  # plt.savefig(f'imagini/shanks_centered_oms.png', bbox_inches='tight')
  #
  # plt.show()
  #
  # # print(f"Classic: gsavg:{avg(giant_steps_c_list)} claisc_time:{avg(times_clasic):.3f}")
  # # print(f"Centered: gsavg:{avg(giant_steps_m_list)} centered_time:{avg(times_centered):.3f}")

def test_shanks_clasic_inter_grumpy(numOfTests: int, numOfBits: list, randSeeds = [0]):
  total_steps_class_list, times_clasic_list, total_memory_clasic = [], [], []
  total_steps_inter_list, times_inter_list, total_memory_inter = [], [], []
  total_steps_grumpy_list, times_grumpy_list, total_memory_grumpy = [], [], []
  for bits in numOfBits:
    seed_steps_class_list, seed_times_clasic_list, seed_total_memory_clasic = [], [], []
    seed_steps_inter_list, seed_times_inter_list, seed_total_memory_inter = [], [], []
    seed_steps_grumpy_list, seed_times_grumpy_list, seed_total_memory_grumpy = [], [], []
    print("Started: ", bits)
    if bits >= 50:
      numOfTests = 50
    elif bits >= 40:
      numOfTests = 100
    elif bits >= 30:
      numOfTests = 500
    for s in randSeeds:
      seed(s)
      time_clas, time_inter, time_grumpy = [], [], []
      steps_clas, steps_inter, steps_grumpy = [], [], []
      mem_clas, mem_inter, mem_grumpy = [], [], []
      for iter in range(numOfTests):
        p, g, exp, h = gnp.logarithm_test_numbers(bits, False)
        t = timer()
        x_c, small_clas, giant_clas, mem_c = shanks_classic(g, h, p, p - 1)
        steps_clas.append(small_clas + giant_clas)
        mem_clas.append(mem_c)
        passed_time = (timer() - t) * 1000
        time_clas.append(passed_time)

        t = timer()
        x_c, oper, mem_c = shanks_interleved(g, h, p, p - 1)
        if x_c != exp:
          print("NI")
        steps_inter.append(oper)
        mem_inter.append(mem_c)
        passed_time = (timer() - t) * 1000
        time_inter.append(passed_time)

        t = timer()
        x_c, oper, mem_c = shanks_two_grumpys_one_baby(g, h, p, p - 1)
        if x_c != exp:
          print("NG")
        steps_grumpy.append(oper)
        mem_grumpy.append(mem_c)
        passed_time = (timer() - t) * 1000
        time_grumpy.append(passed_time)
      seed_steps_class_list.append(avg(steps_clas) / ceil(sqrt(p)))
      seed_steps_inter_list.append(avg(steps_inter) / ceil(sqrt(p)))
      seed_steps_grumpy_list.append(avg(steps_grumpy) / ceil(sqrt(p)))
      seed_times_clasic_list.append(avg(time_clas))
      seed_times_inter_list.append(avg(time_inter))
      seed_times_grumpy_list.append(avg(time_grumpy))
      seed_total_memory_clasic.append(avg(mem_clas) / ceil(sqrt(p)))
      seed_total_memory_inter.append(avg(mem_inter) / ceil(sqrt(p)))
      seed_total_memory_grumpy.append(avg(mem_grumpy) / ceil(sqrt(p)))
    print(f"Avg steps clasic:{avg(seed_steps_class_list):.4f} -- Avg steps inter:{avg(seed_steps_inter_list):.4f} -- Avg steps grumpy:{avg(seed_steps_grumpy_list):.4f}")
    print(f"Avg memory clasic:{avg(seed_total_memory_clasic):.4f} -- Avg memory inter:{avg(seed_total_memory_inter):.4f} -- Avg memory grumpy:{avg(seed_total_memory_grumpy):.4f}")
    print(f"Avg time clasic:{avg(seed_times_clasic_list):.4f} -- Avg time inter:{avg(seed_times_inter_list):.4f} -- Avg time grumpy:{avg(seed_times_grumpy_list):.4f}")
    total_steps_class_list.append(avg(seed_steps_class_list))
    total_steps_inter_list.append(avg(seed_steps_inter_list))
    total_steps_grumpy_list.append(avg(seed_steps_grumpy_list))
    times_clasic_list.append(avg(seed_times_clasic_list))
    times_inter_list.append(avg(seed_times_inter_list))
    times_grumpy_list.append(avg(seed_times_grumpy_list))
    total_memory_clasic.append(avg(seed_total_memory_clasic))
    total_memory_inter.append(avg(seed_total_memory_inter))
    total_memory_grumpy.append(avg(seed_total_memory_grumpy))
    print(f"Finished: {bits}\n")

  print(f"Avg steps clasic:{avg(total_steps_class_list):.4f} -- Avg steps inter:{avg(total_steps_inter_list):.4f} -- Avg steps grumpy:{avg(total_steps_grumpy_list):.4f}")
  print(f"Avg memory clasic:{avg(total_memory_clasic):.4f} -- Avg memory inter:{avg(total_memory_inter):.4f} -- Avg memory grumpy:{avg(total_memory_grumpy):.4f}")
  print(f"Avg time clasic:{avg(times_clasic_list):.4f} -- Avg time inter:{avg(times_inter_list):.4f} -- Avg time grumpy:{avg(times_grumpy_list):.4f}")

def shanks_worst_tests(numOfTests: int, numOfBits: list, randSeed = 0):
  seed(randSeed)
  total_steps_class_list, times_clasic_list, total_memory_clasic = [], [], []
  total_steps_grumpy_list, times_grumpy_list, total_memory_grumpy = [], [], []
  for bits in numOfBits:
    time_clas, time_grumpy = [], []
    steps_clas, steps_grumpy = [], []
    mem_clas, mem_grumpy = [], []
    print("Started: ", bits)
    if bits >= 50:
      numOfTests = 50
    elif bits >= 40:
      numOfTests = 100
    for iter in range(numOfTests):
      p, g, exp, h = gnp.worst_shanks_numbers(bits)
      t = timer()
      x_c, small_clas, giant_clas, mem_c = shanks_classic(g, h, p, p - 1)
      steps_clas.append(small_clas + giant_clas)
      mem_clas.append(mem_c)
      passed_time = (timer() - t) * 1000
      time_clas.append(passed_time)

      t = timer()
      x_c, oper, mem_c = shanks_two_grumpys_one_baby(g, h, p, p - 1)
      if x_c != exp:
        print("NG")
      steps_grumpy.append(oper)
      mem_grumpy.append(mem_c)
      passed_time = (timer() - t) * 1000
      time_grumpy.append(passed_time)
    print(f"Avg steps clasic:{avg(steps_clas) / ceil(sqrt(p)):.4f} -- Avg steps grumpy:{avg(steps_grumpy) / ceil(sqrt(p)):.4f}")
    print(f"Avg memory clasic:{avg(mem_clas) / ceil(sqrt(p)):.4f} -- Avg memory grumpy:{avg(mem_grumpy) / ceil(sqrt(p)):.4f}")
    print(f"Avg time clasic:{avg(time_clas):.4f} -- Avg time grumpy:{avg(time_grumpy):.4f}")
    total_steps_class_list.append(avg(steps_clas) / ceil(sqrt(p)))
    total_steps_grumpy_list.append(avg(steps_grumpy) / ceil(sqrt(p)))
    times_clasic_list.append(avg(time_clas))
    times_grumpy_list.append(avg(time_grumpy))
    total_memory_clasic.append(avg(mem_clas) / ceil(sqrt(p)))
    total_memory_grumpy.append(avg(mem_grumpy) / ceil(sqrt(p)))
    print(f"Finished: {bits}\n")

  print(f"Avg steps clasic:{avg(total_steps_class_list):.4f} -- Avg steps grumpy:{avg(total_steps_grumpy_list):.4f}")
  print(f"Avg memory clasic:{avg(total_memory_clasic):.4f} -- Avg memory grumpy:{avg(total_memory_grumpy):.4f}")
  print(f"Avg time clasic:{avg(times_clasic_list):.4f} -- Avg time grumpy:{avg(times_grumpy_list):.4f}")

def shanks_worst(numOfTests: int, numOfBits: list, randSeed = 0):
  seed(randSeed)
  for bits in numOfBits:
    print("Start: ", bits)
    its = 0
    for i in range(numOfTests):
      # print(f"Iteration:{i}", end=": ")
      p, g, exp, h = gnp.worst_shanks_numbers(bits)
      x, sm, gs, l = shanks_classic(g, h, p, p - 1)
      if gs == sm or gs == sm - 1:
        its += 1
      else:
        print(sm, gs)
    print(f"Fin: {bits}, its:{its}\n")