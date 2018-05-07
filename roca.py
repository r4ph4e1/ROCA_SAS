# coding: utf8
from Crypto.PublicKey import RSA

from fractions import gcd as gcd
import math
import itertools
from functools import reduce
import subprocess
import os

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


def a2(M, n, div):
    M_strich = M
    ord_strich = div
    orders = {}

    for p in n:
        orders[p] = ord(p)

        if orders[p] % ord_strich != 0:
            M_strich /= p
    return M_strich

def choose_divisor(M, Mold, ord, ordold):
    return (math.log(ordold, 2) - math.log(ord, 2)) / (math.log(Mold, 2) - math.log(M, 2))

def greedy_heuristic(n, M):
    ord_M = order(n)
    pfo = prime_factors(ord_M)

    M_new = a2(M, n, ord_M/max(pfo))
    print(M_new)
    ord_new = ord_M / max(pfo)

    div = choose_divisor(M_new, M, ord_new, ord_M)
    print(div)




def get_m(n, limes):
    M = calcM(n)

    greedy_heuristic(n, M)



''' for po in reversed(pf):
        M_strich = M
        for p in reversed(n):
            orders[p] = ord(p)
            if M_strich % po != 0:
                M_strich = M_strich / po
                #print(str(M_strich) + " / " + str(p))
                if math.log(M_strich, 2) < limes:
                    ergs.add((M_strich, po))
                    break
    for erg in ergs:
        print(erg)


    #print(M_strich)
    return M_strich
'''


if __name__ == "__main__":
    with open('tmp.pub', 'r') as f:
        pub_key = RSA.importKey(f.read())
        param = get_param(pub_key.size())
        n = get_primes(param['anz'])


        # print pub_key.n
        limes = math.log(pub_key.n, 2) / 4
        M = get_m(n, limes)
        #y = 83*53*41*29*37*23*17*99
        #print(M/y)

    #roca(pub_key.n, M, param['m'], param['t'])
