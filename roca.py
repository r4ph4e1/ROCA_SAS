import math
from functools import reduce
import subprocess
import os


# Projekterzeichnis
dir = "/home/raphael/Desktop/ROCA_SAS/"

"""
Parameter:
    N = p * q
    M = Produkt der ersten n Primzahlen
    m und t = Optimierungs Parameter fuer Coppersmith
"""
def roca(N, M, m, t):
    # TODO Polynom f(x) erstellen,
    # c = diskreter Logarithmus von N zur Basis 65537 mod M
    # ord = Ordnung von 65537 zur Basis M
    p = 0


    c = call_fingerprint_tool
    ord = order(M, n)
    beta = 0.5
    X = 2 * pow(N, beta) / M

    for a in range(c/2, (c+ord)/2):
        #TODO k = Coppersmith(f(x), N, beta, m, t, X)
        #TODO p = k * M + (65337^a mod M)

        if N % p == 0:
            return p

def dlp():
    return 0

def lcm(numbers):
    lcm = numbers[0]
    for i in numbers[1:]:
        #print(i)
        lcm = int(lcm*i/math.gcd(lcm, i))
    print(lcm)

    return lcm

# to see if the correct prime factos are calculated
def get_primefactors(n):
    result = []
    for i in range(2,n):
        while n % i == 0:
            #print i,"|",n
            n = n/i
            result.append(i)

        if n == 1:
            break

    if n > 1:
        result.append(n)

    print(result)

def order(M, pi):
    generator = 65537
    ord_pi = []

    #ordPi = ord_pi (65537)
    for i in pi:
        for j in range(1, i):
            #print(str(j))
            if generator ** j % i == 1:
                ord_pi.append(j)
                print("ord" + str(i) + "(65537) = " + str(j))
                break
            else:
                continue

    ord_m = lcm(ord_pi)
    print(M)
    print(generator)
    print(ord_m)

    return ord_m

# test for deduplicate pi orders
# to match wolframalpha value
def deduplicate(ord_list):

    dup = set(ord_list)
    ord_list = list(dup)
    return ord_list

def checkprime(x):
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

    return n

def calcM(n):
    M = 1
    for i in range(0, len(n)):
        M = M * n[i]

    print('-----------M------------------')
    print(M)
    return M


# TODO muss noch überarbeitet werden, damit pubkey übergeben werden kann und roca fingerprint ausgeführt werden kann.
def call_fingerprint_tool():
    #path = dir + 'roca_fingerprint/detect.py '
    #subprocess.Popen(['python', path, pubkey], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #subprocess.run('python' + path + 'roca_fingerprint/detect.py ' + pubkey)

    #if not os.path.exists(dir + "tmp.txt"):
    #    open(dir + "tmp.txt", 'w+')


    f = open(dir + "tmp.txt", 'r+')
    c = f.readline()

    print(c)

    return c

if __name__ == "__main__":
    #TODO Parameter für roca übergeben
    roca(N, M, m, t)
