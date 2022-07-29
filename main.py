import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from fractions import Fraction
from random import randint
import random
import math
import time


count_bern = 0

def unif_sample(Low=0, S=2**10, bits=2**32):
    count = 0
    x = np.random.randint(low = 0, high = bits)
    count += 1
    m = x * S
    l = m % bits
    if l < S:
        t = (bits-S) % S
        while l < t :
            x = np.random.randint(low=0, high=bits)
            count += 1
            m = x * S
            l = m % bits

    return (m / bits), count

def trunc_geom():
    p = .3
    r = np.random.rand(1000)
    N = 5
    gen_trunc = np.floor(np.log(1 - r * (1 - (1 - p) ** N)) / np.log(1 - p))

    return gen_trunc

def sample_uniform(m,rng):
    assert isinstance(m,int)
    assert m>0

    return rng.randrange(m)

def sample_bernoulli(p,rng=random.SystemRandom()):
    assert isinstance(p,Fraction)
    assert 0 <= p <= 1

    # global count_bern
    #
    # count_bern += 1

    m=sample_uniform(p.denominator,rng)
    if m < p.numerator:
        return 1
    else:
        return 0

def sample_bernoulli_exp1(x,rng):
    assert isinstance(x,Fraction)
    assert 0 <= x <= 1

    k=1
    while True:
        if sample_bernoulli(x/k,rng)==1:
            k=k+1
        else:
            break
    return k%2

def sample_bernoulli_exp(x,rng = None):
    assert isinstance(x,Fraction)
    assert x >= 0

    if rng is None:
        rng = random.SystemRandom()

    while x>1:
        if sample_bernoulli_exp1(Fraction(1,1),rng)==1:
            x=x-1
        else:
            return 0
    return sample_bernoulli_exp1(x,rng)

def sample_fancy_trun_geometric(x,rng = None):
    assert isinstance(x,Fraction)
    if x==0: return 0
    assert x>0

    if rng is None:
        rng = random.SystemRandom()

    t=x.denominator
    while True:
        u=sample_uniform(t,rng)
        b=sample_bernoulli_exp(Fraction(u,t),rng)
        if b==1:
            return u


# sample from a geometric(1-exp(-x)) distribution
# assumes x is a rational number >= 0
def sample_geometric_exp_slow(x, rng=random.SystemRandom()):
    assert isinstance(x, Fraction)
    assert x >= 0
    k = 0
    while True:
        if sample_bernoulli_exp(x, rng) == 1:
            k = k + 1
        else:
            return k


def sample_geometric_exp_fast(epsgamma, t, rng = None):
    # assert isinstance(x, Fraction)
    # if x == 0: return 0  # degenerate case
    # assert x > 0

    # t = x.denominator

    if rng is None:
        rng = random.SystemRandom()

    v = sample_geometric_exp_slow(Fraction(1, 1), rng)

    while True:
        u = sample_uniform(t, rng)
        b = sample_bernoulli_exp(Fraction(u, t), rng)
        if b == 1:
            break
    capx = u + t * v # geom(1-exp(-1/t))
    return capx // epsgamma # geom(1-exp(-s/t))


def get_1et(gamma):
    return math.exp(gamma)

def get_1ef(gamma):
    return 1/(1-math.exp(-gamma))

def get_et(gamma):
    #return math.exp(gamma)
    return (np.exp(gamma)/2 + (np.exp(-gamma)/2) - gamma*np.exp(-gamma))*np.exp(gamma)
    #return (np.exp(gamma)/2 - (np.exp(-gamma)/2))*np.exp(gamma)

def get_ef(gamma):
    #return 1/(1-math.exp(-gamma))
    return (np.exp(gamma) / 2 - (np.exp(-gamma) / 2) + (gamma) * np.exp(-gamma) )*np.exp(gamma)
    #return (np.exp(gamma)/2 + (np.exp(-gamma)/2))*np.exp(gamma)

def get_egammadot(x):
    return math.exp(x-math.floor(x))

def get_eut(t):
    return (math.e-1)/(t*(math.e-math.e**(1-(1/t))))


def F_Y_shuffle(arr, n):
    for i in range(n - 1, 0, -1):
        j = randint(0, i + 1)
        arr[i], arr[j] = arr[j], arr[i]
    return arr

def gaptopk_secure(Q,K,eps):

        # gamma* by designer, accuracy is constant

        Qhat = []
        gamma = 1

        begin_time=time.time()

        for q in Q:
            Y0 = sample_geometric_exp_fast(eps*gamma, 2*K)
            Qhat.append(q + gamma*Y0)

        gamma0 = gamma
        t = 1
        tie = True
        M = 10

        while tie:
            gammat = gamma/M
            for i in range(len(Qhat)):
                Yt = sample_geometric_exp_fast(eps*gammat, 2*K)
                Qhat[i] = Qhat[i] + gammat*(Yt % M)

            Qhat_sorted = Qhat.copy()
            j = np.argpartition(Qhat_sorted,-(K+2))[-(K+2):]

            for i in range(K+1):
                if Qhat[j[i]] == Qhat[j[i+1]]:
                    break
                tie = False
            t = t+1

        x = F_Y_shuffle(list(range(K+1)), K+1)

        g=[]
        for i in range(K):
            if x[i] < x[i+1]:
                g.append(math.floor(Qhat[j[i]] - Qhat[j[i+1]] - gammat))
            else:
                g.append(math.floor(Qhat[j[i]] - Qhat[j[i+1]]))

        result = []
        for k in range(K):
            result.append((j[i],g[i]))

        return result

if __name__=='__main__':

    Q = list(np.arange(0.0,1000.0,0.1))
    print(gaptopk_secure(Q, 2, 1))



    ###### Bern and Geom plotting below ######

    # graphs for sample_bernoulli_exp with \gamma in [0,1]#

    # num_sample = list()
    # int_num = list()
    # sample_result = list()
    # sample_result_out = list()
    # for gamma in range(1,1000): # range for gamma  to / by 50000
    #     current_count_total=0
    #     first_true = False
    #     for j in range(1000): # how many trials to average per gamma, modify to change smoothness of graph
    #         bern_result = sample_bernoulli_exp(Fraction(gamma,1000))
    #         if bern_result == 0:
    #             current_count_total += count
    #         if bern_result == 1 and first_true is False:
    #             current_count_total += count
    #             sample_result_out.append(current_count_total)
    #             first_true = True
    #         int_num.append(count)
    #         count = 0
    #     num_sample.append(math.log2(np.mean(int_num)))
    #     int_num.clear()

    # plt.plot([x/50 for x in range(1,50)], num_sample)
    # plt.plot([x / 50 for x in range(1, 50)], [math.log2(math.exp(x / 50)) for x in range(1, 50)])
    # plt.legend(['empirical', 'expectation'])
    # plt.xlabel('Gamma')
    # plt.ylabel('log2 of unif rand samples needed')
    # plt.show()

    # plt.plot([x / 50 for x in range(1, 50)], [math.log2(math.exp(x / 50)) for x in range(1, 50)])
    # plt.xlabel('x')
    # plt.ylabel('y = log2(e^x)')
    # plt.show()

    ################################################################################################

    # graphs for sample_bernoulli_exp first T and F


    # sample_result_out = list()
    # for gamma in range(1,100): # range for gamma  to / by 100
    #     sample_result = list()
    #     for j in range(100): # how many trials to average per gamma, modify to change smoothness of graph
    #         current_count_total = 0
    #         first_true = False
    #         while first_true == False:
    #             bern_result = sample_bernoulli_exp(Fraction(gamma, 100))
    #             if bern_result == 1: # get first true =1 or false =0
    #                 current_count_total += count
    #                 count = 0
    #                 sample_result.append(current_count_total)
    #                 first_true = True
    #
    #     sample_result_out.append(np.mean(sample_result))


    # plt.plot([x / 100 for x in range(1, 100)], sample_result_out)
    # plt.plot([x / 100 for x in range(1, 100)], [ np.exp(x/100) for x in range(1, 100)]) # mine
    # # plt.plot([x / 100 for x in range(1, 100)], [ np.exp(x/100)/2 + np.exp(-x/100)/2 -(x/100)*np.exp(-x/100) for x in range(1, 100)]) # kifer
    # # plt.plot([x / 100 for x in range(1, 100)], [ np.exp(2*x/100)/2 + 1/2 - np.exp(x/100) for x in range(1, 100)]) # kifer*egamma
    # plt.legend(['empirical', 'expectation'])
    # plt.xlabel('Gamma')
    # plt.ylabel('Calls to r $Bern(\gamma/k)$ needed to get first True')
    # plt.show()

    # plt.plot([x / 100 for x in range(20, 100)], sample_result_out)
    # plt.plot([x / 100 for x in range(20, 100)], [1/(1-np.exp(-x/100))  for x in range(20, 100)])
    # #plt.plot([x / 100 for x in range(1, 200)], [0.5*np.exp(2*x / 100) - 0.5 + (x / 100) for x in range(1, 200)]) #kifer
    # plt.legend(['empirical', 'expectation'])
    # plt.xlabel('Gamma')
    # plt.ylabel('Calls to r $Bern(\gamma/k)$ needed to get first False')
    # plt.show()

    # plt.plot([x / 50 for x in range(1, 50)], [(np.exp(x/50)-np.exp(-(x/50)))/2 for x in range(1, 50)])
    # plt.xlabel('Gamma')
    # plt.ylabel('expected samples of getting first True')
    # plt.show()

    # plt.plot([x / 50 for x in range(1, 50)], [(np.exp(x/50)+np.exp(-(x/50)))/2 for x in range(1, 50)])
    # plt.xlabel('Gamma')
    # plt.ylabel('expected samples of getting first False')
    # plt.show()

    ################################################################################################

    # graphs for sample_bernoulli_exp with \gamma in (1,\infty)#

    # num_sample = list()
    # int_num = list()
    # sample_result = list()
    # for i in range(51, 301):  # range for gamma
    #     out = 0
    #     for j in range(2000):  # how many trials to average per gamma, modify to change smoothness of graph
    #         out += sample_bernoulli_exp(Fraction(i, 50))
    #         int_num.append(count)
    #         count = 0
    #     num_sample.append(np.mean(int_num))
    #     sample_result.append(out)
    #     int_num.clear()
    #
    # plt.plot([x / 50 for x in range(51, 301)], num_sample)
    # plt.plot([x / 50 for x in range(51, 301)],
    #          [get_ef(1)*(1-math.exp(-1)) + (get_et(1)+get_egammadot(x))/math.e
    #           for x in [y / 50 for y in range(51, 100)]]+
    #
    #          [get_ef(1)*(1-1/math.e) +
    #           ((1-1/math.e)*math.exp(-(math.floor(x)-1))*(math.exp(math.floor(x))*(get_et(1)+get_ef(1))
    #                                                       -math.e*((math.floor(x)-1)*get_et(1)+get_et(1)+get_ef(1))
    #                                                       +(math.floor(x)-1)*get_et(1)-math.exp(math.floor(x)-1)*get_et(1)+get_ef(1))
    #           / (math.e-1)**2)
    #                     +(math.floor(x)*get_et(1)+get_egammadot(x))/ math.exp(math.floor(x))
    #           for x in [y / 50 for y in range(100, 150)]]+
    #
    #          [get_ef(1) * (1 - 1 / math.e) +
    #           ((1 - 1 / math.e) * math.exp(-(math.floor(x) - 1)) * (math.exp(math.floor(x)) * (get_et(1) + get_ef(1))
    #                                                                 - math.e * ((math.floor(x) - 1) * get_et(
    #                       1) + get_et(1) + get_ef(1))
    #                                                                 + (math.floor(x) - 1) * get_et(1) - math.exp(
    #                       math.floor(x) - 1) * get_et(1) + get_ef(1))
    #            / (math.e - 1) ** 2)
    #           + (math.floor(x) * get_et(1) + get_egammadot(x)) / math.exp(math.floor(x))
    #           for x in [y / 50 for y in range(150, 200)]]+
    #
    #          [get_ef(1) * (1 - 1 / math.e) +
    #           ((1 - 1 / math.e) * math.exp(-(math.floor(x) - 1)) * (math.exp(math.floor(x)) * (get_et(1) + get_ef(1))
    #                                                                 - math.e * ((math.floor(x) - 1) * get_et(
    #                       1) + get_et(1) + get_ef(1))
    #                                                                 + (math.floor(x) - 1) * get_et(1) - math.exp(
    #                       math.floor(x) - 1) * get_et(1) + get_ef(1))
    #            / (math.e - 1) ** 2)
    #           + (math.floor(x) * get_et(1) + get_egammadot(x)) / math.exp(math.floor(x))
    #           for x in [y / 50 for y in range(200, 301)]]
    #          )
    # plt.legend(['empirical', 'expectation'])
    # plt.xlabel('gamma')
    # plt.ylabel('Number of Bern(X) called')
    # plt.show()

    #
    # plt.xlabel('gamma')
    # plt.ylabel('expected amount of sample to return when gamma > 1')
    # plt.show() 

    # plt.plot([x / 500 for x in range(10, 5000)], sample_result)
    # plt.ylabel('sum of bernoulli output')
    # plt.show()

    ################################################################################################

    # graphs for sample_fancy_trun_geometric#

    # num_sample = list()
    # int_num = list()
    # sample_result = list()
    # sample_result_out = list()
    # for i in range(1,10): # range for gamma
    #     out=0
    #     first_false = False
    #     for j in range(100000): # how many trials to average per gamma, modify to change smoothness of graph
    #         bern_result = sample_fancy_trun_geometric(Fraction(i+1,i))
    #         out+=1
    #         if bern_result == 0 and first_false == False:
    #             sample_result.append(out)
    #             first_false = True
    #         int_num.append(count)
    #         count = 0
    #     sample_result_out.append(np.average(sample_result))
    #     num_sample.append(math.log2(np.mean(int_num)))
    #     int_num.clear()
    #
    # plt.plot([x for x in range(1,10)], num_sample)
    # plt.xlabel('t')
    # plt.ylabel('log2 of unif rand samples needed')
    # plt.show()

    # plt.plot([x / 50 for x in range(1, 50)], [math.log2(math.exp(x / 50)) for x in range(1, 50)])
    # plt.xlabel('x')
    # plt.ylabel('y = log2(e^x)')
    # plt.show()

    # plt.plot([x / 50000 for x in range(1000, 50000)], sample_result_out)
    # plt.xlabel('Gamma')
    # plt.ylabel('Sample needed to get first False')
    # plt.show()

    # plt.plot([x / 50 for x in range(1, 50)], [(np.exp(x/50)+np.exp(-x/50))/2 for x in range(1, 50)])
    # plt.xlabel('Gamma')
    # plt.ylabel('expected samples of getting first True')
    # plt.show()

    ################################################################################################

    # trange = 300
    # pergamma = 50000
    #
    # num_sample = list()
    # int_num = list()
    # # sample_result = list()
    # for t in range(1, trange):  # range for gamma
    #     # out = 0
    #     for j in range(pergamma):  # how many trials to average per gamma, modify to change smoothness of graph
    #         sample_geometric_exp_fast(Fraction(1, t))
    #         # out += sample_geometric_exp_fast(Fraction(1, t))
    #         int_num.append(count_bern)
    #         count_bern = 0
    #     num_sample.append(np.mean(int_num))
    #     # sample_result.append(out)
    #     int_num.clear()
    #
    # plt.plot([x for x in range(1, trange)], num_sample)
    # # plt.legend(['empirical', 'expectation'])
    # plt.xlabel('t')
    # plt.ylabel('Number of Bern called')
    # plt.show()

    ################################################################################################

    # running time comparison

    # trange = 100
    # pergamma = 10
    #
    # num_sample = list()
    # int_num = list()
    # # sample_result = list()
    # start_time = time.time()
    # for t in range(1, trange):  # range for gamma
    #     # out = 0
    #     for j in range(pergamma):  # how many trials to average per gamma, modify to change smoothness of graph
    #         sample_geometric_exp_slow(Fraction(1, t))
    #         # out += sample_geometric_exp_fast(Fraction(1, t))
    #         int_num.append(count_bern)
    #         count_bern = 0
    #     num_sample.append(np.mean(int_num))
    #     # sample_result.append(out)
    #     int_num.clear()
    #
    # geom_time = time.time() - start_time
    #
    # start_time = time.time()
    #
    # for t in range(1, trange):  # range for gamma
    #     # out = 0
    #     for j in range(pergamma):  # how many trials to average per gamma, modify to change smoothness of graph
    #         sample_bernoulli(Fraction(1, 1+t))
    #         # out += sample_geometric_exp_fast(Fraction(1, t))
    #         int_num.append(count_bern)
    #         count_bern = 0
    #     num_sample.append(np.mean(int_num))
    #     # sample_result.append(out)
    #     int_num.clear()
    #
    # bern_time = time.time() - start_time
    #
    #
    # print(geom_time, bern_time)



    # plt.plot([x for x in range(1, trange)], num_sample)
    # # plt.legend(['empirical', 'expectation'])
    # plt.xlabel('t')
    # plt.ylabel('Number of Bern called')
    # plt.show()

