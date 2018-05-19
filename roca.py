# coding: utf8
from Crypto.PublicKey import RSA

from fractions import gcd as gcd
import math
import itertools
from functools import reduce
import subprocess
import os
import binascii

# Projekterzeichnis

DEBUG = False

dir = "/home/raphael/Desktop/ROCA_SAS/"

"""
Parameter:
    N = p * q
    M = Produkt der ersten n Primzahlen
    m und t = Optimierungs Parameter fuer Coppersmith
"""


def roca(N, M, m, t):
    # N = p * q
    # M = produkt der ersten n primzahlen
    # m =
    # t =
    # M_strich =

    # TODO Polynom f(x) erstellen,
    # c = diskreter Logarithmus von N zur Basis 65537 mod M
    c_strich = math.log(N, 65537)
    # ord = Ordnung von 65537 zur Basis M
    p = 0

    # Todo: script einbinden
    c = call_fingerprint_tool
    ord = order(n)
    beta = 0.5
    X = 2 * pow(N, beta) / M

    for a in range(c / 2, (c + ord) / 2):
        # TODO k = Coppersmith(f(x), N, beta, m, t, X)
        # TODO p = k * M + (65337^a mod M)

        if N % p == 0:
            return p


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
    generator = 65537
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


# TODO muss noch überarbeitet werden, damit pubkey übergeben werden kann und roca fingerprint ausgeführt werden kann.
def call_fingerprint_tool():
    # path = dir + 'roca_fingerprint/detect.py '
    # subprocess.Popen(['python', path, pubkey], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # subprocess.run('python' + path + 'roca_fingerprint/detect.py ' + pubkey)

    # if not os.path.exists(dir + "tmp.txt"):
    #    open(dir + "tmp.txt", 'w+')

    f = open(dir + "tmp.txt", 'r+')
    c = f.readline()

    print(c)

    return c


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
    M_strich = M

    for p in reversed(pfo):
        # ord_pi teilt nicht ord_strich
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
    ord_M = order(n)
    pfo = prime_factors(ord_M)
    pf_M = n
    M_old = M
    ord_new = ord_M

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

        # TODO Händisch verifizieren ob die Werte nach der ersten Runde passen.

        for p in reversed(pfo):  # 53 ist in der Ordnung dabei und im Paper nicht, warum?
            pf_M_tmp = list(pf_M)
            M_new, pf_M_tmp = a2(M_old, pf_M_tmp, ord_new / p)  # Kandidat für M_strich
            # print("M NEW: " + str(M_new))
            print(pf_M)
            div = choose_divisor(M_new, M_old, ord_new / p, ord_new)

            div_dict[p] = (div, M_new, pf_M_tmp)
            # print("Div: " + str(div))

            # print(div_dict)

        best_candidate = max(div_dict, key=div_dict.get)
        # print(best_candidate)
        #print(div_dict)

        ord_new /= best_candidate
        M_old = div_dict[best_candidate][1]
        pfo.remove(best_candidate)
        pf_M = div_dict[best_candidate][2]
        #print(pf_M)
        print("best candidate:" + str(best_candidate))
        print("M Strich nach Runde %d: %d" % (runde, M_old))
        print("ORD NEW: " + str(ord_new))
        print("PRIME Factors: " + str(pfo))
        runde += 1
        print('\n')



def get_m(n, limes):
    M = calcM(n)

    greedy_heuristic(n, M, limes)


if __name__ == "__main__":
    with open('tmp.pub', 'r') as f:
        pub_key = RSA.importKey(f.read())
        param = get_param(pub_key.size())
        n = get_primes(param['anz'])
        M = calcM(n)

        limes = math.log(pub_key.n, 2) / 4
        greedy_heuristic(n, M, limes)

    # roca(pub_key.n, M, param['m'], param['t'])
