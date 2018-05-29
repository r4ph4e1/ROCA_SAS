
# coding: utf8
import pdb
from Crypto.PublicKey import RSA
import time
import multiprocessing
from fractions import gcd as gcd
import math
from multiprocessing import Process, Pool
from sage.all_cmdline import *
#from roca_fingerprint.detect import RocaFingerprinter
import itertools
from functools import reduce
import subprocess
import os
import binascii

# Projekterzeichnis

DEBUG = False

dir = "/home/raphael/Desktop/ROCA_SAS/"
start_time = time.time()
"""
Parameter:
    N = p * q
    M = Produkt der ersten n Primzahlen
    m und t = Optimierungs Parameter fuer Coppersmith
"""

def get_end(c, ord, id):
    start = int(c) / int(2)
    end = (c + ord) / 2
    count = end - start
    cpus = multiprocessing.cpu_count()
    div = floor(count / cpus)
    return int(start + div * (id + 1) - 1)


def get_start(c, ord, id):
    start = int(c) / int(2)
    end = (c + ord) / 2
    count = end - start
    cpus = multiprocessing.cpu_count()
    div = floor(count / cpus)
    return int(start + div * id)


def worker(args):
    #{'cpu': rest, 'n': n, 'M_strich': M_strich,'m': m, 't': t, 'c': c, 'ord_new': ord_new}
    id = args['cpu']
    pub_key = args['n']
    M_strich = args['M_strich']
    t = args['t']
    c = args['c']
    ord = args['ord_new']
    m = args['m']
    n = pub_key.n

    beta = 0.5
    X = 2 * pow(n, beta) / M
    start = get_start(c, ord, id)
    end = get_end(c, ord, id)
    ZmodN = Zmod(n)

    for a_strich in xrange(start, end):
        R.<x> = PolynomialRing(ZmodN)

        invers = inverse_mod(int(M_strich), n)
        pol = x + invers * int(Integer(65537).powermod(a_strich, M_strich))

        roots = coppersmith_howgrave_univariate(pol, n, beta, m, t, X)
        for root in roots:
            p = root*M_strich + int(Integer(65537).powermod(a_strich, M_strich))
            if Mod(n, p) == 0:
                print("--- %s seconds ---" % (time.time() - start_time))
                print("Success p: %d " % p)
                break


    return


def coppersmith_howgrave_univariate(pol, modulus, beta, mm, tt, XX):
    """
    Coppersmith revisited by Howgrave-Graham

    finds a solution if:
    * b|modulus, b >= modulus^beta , 0 < beta <= 1
    * |x| < XX
    """
    #
    # init
    #
    #dd = pol.degree()
    dd = 1
    nn = dd * mm + tt

    #
    # checks
    #
    if not 0 < beta <= 1:
        raise ValueError("beta should belongs in (0, 1]")

    if not pol.is_monic():
        raise ArithmeticError("Polynomial must be monic.")

    #
    # calculate bounds and display them
    #
    """
    * we want to find g(x) such that ||g(xX)|| <= b^m / sqrt(n)
    * we know LLL will give us a short vector v such that:
    ||v|| <= 2^((n - 1)/4) * det(L)^(1/n)
    * we will use that vector as a coefficient vector for our g(x)

    * so we want to satisfy:
    2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)

    so we can obtain ||v|| < N^(beta*m) / sqrt(n) <= b^m / sqrt(n)
    (it's important to use N because we might not know b)
    """
    debug = False
    if debug:
        # t optimized?
        print "\n# Optimized t?\n"
        print "we want X^(n-1) < N^(beta*m) so that each vector is helpful"
        cond1 = RR(XX ^ (nn - 1))
        print "* X^(n-1) = ", cond1
        cond2 = pow(modulus, beta * mm)
        print "* N^(beta*m) = ", cond2
        print "* X^(n-1) < N^(beta*m) \n-> GOOD" if cond1 < cond2 else "* X^(n-1) >= N^(beta*m) \n-> NOT GOOD"

        # bound for X
        print "\n# X bound respected?\n"
        print "we want X <= N^(((2*beta*m)/(n-1)) - ((delta*m*(m+1))/(n*(n-1)))) / 2 = M"
        print "* X =", XX
        cond2 = RR(modulus ^ (((2 * beta * mm) / (nn - 1)) - ((dd * mm * (mm + 1)) / (nn * (nn - 1)))) / 2)
        print "* M =", cond2
        print "* X <= M \n-> GOOD" if XX <= cond2 else "* X > M \n-> NOT GOOD"

        # solution possible?
        print "\n# Solutions possible?\n"
        detL = RR(modulus ^ (dd * mm * (mm + 1) / 2) * XX ^ (nn * (nn - 1) / 2))
        print "we can find a solution if 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n)"
        cond1 = RR(2 ^ ((nn - 1) / 4) * detL ^ (1 / nn))
        print "* 2^((n - 1)/4) * det(L)^(1/n) = ", cond1
        cond2 = RR(modulus ^ (beta * mm) / sqrt(nn))
        print "* N^(beta*m) / sqrt(n) = ", cond2
        print "* 2^((n - 1)/4) * det(L)^(1/n) < N^(beta*m) / sqrt(n) \n-> SOLUTION WILL BE FOUND" if cond1 < cond2 else "* 2^((n - 1)/4) * det(L)^(1/n) >= N^(beta*m) / sqroot(n) \n-> NO SOLUTIONS MIGHT BE FOUND (but we never know)"

        # warning about X
        print "\n# Note that no solutions will be found _for sure_ if you don't respect:\n* |root| < X \n* b >= modulus^beta\n"

    #
    # Coppersmith revisited algo for univariate
    #

    # change ring of pol and x
    polZ = pol.change_ring(ZZ)
    x = polZ.parent().gen()

    # compute polynomials
    gg = []
    for ii in range(mm):
        for jj in range(dd):
            gg.append((x * XX) ** jj * modulus ** (mm - ii) * polZ(x * XX) ** ii)
    for ii in range(tt):
        gg.append((x * XX) ** ii * polZ(x * XX) ** mm)

    # construct lattice B
    BB = Matrix(ZZ, nn)

    for ii in range(nn):
        for jj in range(ii + 1):
            BB[ii, jj] = gg[ii][jj]

    # display basis matrix
    #if debug:
    #    matrix_overview(BB, modulus ^ mm)

    # LLL
    BB = BB.LLL()

    # transform shortest vector in polynomial
    new_pol = 0
    for ii in range(nn):
        new_pol += x ** ii * BB[0, ii] / XX ** ii

    # factor polynomial
    potential_roots = new_pol.roots()
    #print "potential roots:", potential_roots

    # test roots
    roots = []
    for root in potential_roots:
        if root[0].is_integer():
            result = polZ(ZZ(root[0]))
            if gcd(modulus, result) >= modulus ^ beta:
                roots.append(ZZ(root[0]))

    #
    return roots


def lcm(numbers):
    lcm = numbers[0]
    for i in numbers[1:]:
        # print(i)
        lcm = int(lcm * i / gcd(lcm, i))

    return lcm


def ord(i):
    generator = 65537
    ord_pi = []

    for j in range(1, i):
        # print(str(j))
        if generator ** j % i == 1:
            ord_pi.append(j)
            # print("order " + str(i) + " (65537) = " + str(j))
            return j
        else:
            continue


def order(pi):
    ord_pi = []

    # ordPi = ord_pi (65537)
    for i in pi:
        ord_pi.append(ord(i))

    ord_m = lcm(ord_pi)
    if DEBUG:
        print ord_m
    return ord_m


def get_primes(x):
    i = 0
    n = [2]

    while len(n) < x:
        i += 1
        for j in range(2, i):
            if math.fmod(i, j) != 0:
                if j == i - 1:
                    n.append(i)
            else:
                break
    if DEBUG:
        print "First n primes: "
        print(n)
    return n


def calcM(n):
    M = 1
    for i in n:
        M = M * i

    if DEBUG:
        print('-----------M------------------')
        print(M)
    return M


def get_param(key_size):
    if key_size < 510:
        return 0
    elif key_size < 961:
        return {'anz': 39, 'm': 4, 't': 5}
    elif key_size < 992:
        return 0
    elif key_size < 1953:
        return {'anz': 71, 'm': 6, 't': 7}
    elif key_size < 1984:
        return 0
    elif key_size < 3937:
        return {'anz': 126, 'm': 25, 't': 26}
    elif key_size < 3968:
        return 0
    elif key_size < 4097:
        return {'anz': 225, 'm': 7, 't': 8}
    return 0


def prime_factors(n):
    primfac = []
    d = 2
    while d * d <= n:
        while (n % d) == 0:
            primfac.append(d)  # supposing you want multiple factors repeated
            n //= d
        d += 1
    if n > 1:
        primfac.append(n)
    return primfac


# M, Primfaktoren von M und Kandidat für Ordnung
def a2(M, pfo, ord_strich):
    #M_strich = ZZ(M)
    M_strich = M
    print(pfo)
    for p in reversed(pfo):
        # ord_pi teilt nicht ord_strich
        #print("Ord_strich in A2: %d" % ord_strich)
        #print("Ordnung von %d in A2: %d" % (p, ord(p)))
        if ord_strich % ord(p) != 0:

            M_strich /= p
            pfo.remove(p)

    return M_strich, pfo


def choose_divisor(M, Mold, ord, ordold):
    try:
        erg = (math.log(ordold, 2) - math.log(ord, 2)) / (math.log(Mold, 2) - math.log(M, 2))
    except ZeroDivisionError:
        erg = 0

    return erg


def greedy_heuristic(n, M, limes):
    DEBUG = True

    ord_M = order(n)
    pfo = prime_factors(ord_M) # Ordnung der einzelnen Primfaktoren der Ordnung M
    pf_M = n # Primfaktoren von M
    M_old = M
    ord_new = ord_M

    # Fügt die Potenzen der Primfaktoren hinzu zB 2⁴, 3⁴, etc.
    for j in pfo:
        count = 0
        for k in range(0, len(pfo)):
            if pfo[k] == j:
                count += 1
                pfo[k] = pow(j, count)

    runde = 1

    while math.log(M_old, 2) > limes:

        div_dict = {}
        removed = []
        # print("Primfaktor von Kandidat M: " + str(pf_M))
        # print("Primfaktor von Ordnung von M: " + str(pfo))

        # Iteriert durch alle Ordnungen der Primfaktoren von der Ordnung von M
        # Berechnet alle möglichen M_Strichs für die auswahl des besten Kandidaten für M_Strich
        for p in reversed(pfo):
            pf_M_tmp = list(pf_M)
            M_new, pf_M_tmp = a2(M_old, pf_M_tmp, ord_new / p)  # Kandidat für M_strich
            # print("M NEW: " + str(M_new))
            # print(pf_M)
            div = choose_divisor(M_new, M_old, ord_new / p, ord_new)

            div_dict[p] = (div, M_new, pf_M_tmp)
            # print("Div: " + str(div))

            # print(div_dict)

        best_candidate = max(div_dict, key=div_dict.get)
        # print(best_candidate)
        # print(div_dict)

        ord_new /= best_candidate
        M_old = div_dict[best_candidate][1]
        pfo.remove(best_candidate)
        pf_M = div_dict[best_candidate][2]
        if DEBUG:
            print("best candidate:" + str(best_candidate))
            print("M Strich nach Runde %d: %d" % (runde, M_old))
            print("ORD NEW: " + str(ord_new))
            print("PRIME Factors: " + str(pfo))
            print("Bitlength of M_strich %d" % int(M_old).bit_length())
            runde += 1
            print('\n')

    return M_old, ord_new


def get_m(n, limes):
    M = calcM(n)

    greedy_heuristic(n, M, limes)

def parm( n, M_strich, m, t, c, ord_new):
    rest = multiprocessing.cpu_count()

    while rest > 0:
        yield {'cpu': rest, 'n': n, 'M_strich': M_strich,'m': m, 't': t, 'c': c, 'ord_new': ord_new}
        rest -= 1

def fingerprint(M, n):
    try:
        b = Mod(65537, M)
        c = discrete_log(n, b)
        print("The Key is vulnerable to ROCA!")
        return 1
    except ValueError:
        print("The Key is resistant to ROCA!")
        return 0

if __name__ == "__main__":
    with open('tmp.pub', 'r') as f:
        pub_key = RSA.importKey(f.read())

        print "Start Zeit: %f" % start_time
        param = get_param(pub_key.size())

        n = get_primes(param['anz'])
        M = calcM(n)

        # Checks if PubKey is Vulnerable to ROCA
        if fingerprint(M, pub_key.n) == 1:
            limes = math.log(pub_key.n, 2) / 4

            M_strich, ord_new = greedy_heuristic(n, M, limes)
            print("Bitlength of M %d" %int(M).bit_length())
            print("Bitlength of M_strich %d" % int(M_strich).bit_length())

            threads = []
            b = Mod(65537, M_strich)
            c = discrete_log(pub_key.n, b)
        else:
            print("Resistant Key: Terminating execution!")

            #print(pub_key.n)
            #p = Pool()
            #p.map(worker, parm(pub_key, M_strich,  param['m'], param['t'], c, ord_new))

            #worker({'cpu': 0, 'n': pub_key, 'M_strich': M_strich, 'm': param['m'], 't': param['t'], 'c': c, 'ord_new': ord_new})
