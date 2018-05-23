from sage.doctest.util import Timer
ttry = Timer()
t = Timer()

L = 27771430913146044712156219115012732149015337058745243774375474371978395728107173008782747458575903820497344261101333156469136833289328084229401057505005215261077328417649807720533310592783171487952296983742789708502518237023426083874832018749447215424764928016413509553872836856095214672430
L *= 701 # if 701 is included
g = Mod(65537,L)

pmin = 12*2^1020
pmax = 13*2^1020
proof.arithmetic(False) # do not bother proving primality

t.start()
u = lift(g^randrange(L))
while True:
  p = u + randrange(ceil(pmin/L),floor(pmax/L)) * L
  if p.is_prime(): break
print 'time for first prime',t.stop().cputime

t.start()
u = lift(g^randrange(L))
while True:
  q = u + randrange(ceil(pmin/L),floor(pmax/L)) * L
  if q.is_prime(): break
print 'time for second prime',t.stop().cputime



#n = 3625190506732597997596942780169386867073561084337530686734832416152327831391665940242449464840635004205892996719944945548184675893828257084981274302452421
n = p*q
print 'public key',n

smooth,k,H = 41902660800,5,7*2^461
m = 2*k+1

print 'smooth',smooth
print 'multiplicity',k
print 'lattice rank',m
print 'Hbits',H.nbits()
print 'H',H

def lg(x): return log(1.0 * x) / log(2.0)

lggamma = lg(m)/(2.0*k) + lg(2*H)*(m-1.0)/(2.0*k) + lg(n)*((k+1.0)/(2.0*m)-1)
gap = lg(pmin) - lg(n) - lggamma
print 'gap',gap
print 'guaranteed',gap > 0

def smoothorder(l):
  return smooth % Mod(g,l).multiplicative_order() == 0

v = prod(l for l,e in factor(L) if smoothorder(l))
print 'v',v
u = p % v
print 'p residue class',u

ttry.start()

pmin = max(pmin,ceil(n / pmax))
pmax = min(pmax,floor(n / pmin))
pradius = (pmax - pmin) // 2
print 'n-dependent pradius bits',pradius.nbits()

R.<x> = ZZ[]
u += floor((pmin + H * v - u) / v) * v
# now pmin-v < u-H*v <= pmin

t.start()
wu = lift(u/Mod(v,n))
f = wu+H*x
fpowers = [f^0]
for j in range(k): fpowers += [fpowers[j] * f]
print 'fpowers time',t.stop().cputime,'bits',[fi.nbits() for fi in fpowers[k].coefficients(sparse=False)]

X = matrix([n])
for i in range(1,k):
  t.start()
  X = X.augment(vector(ZZ,i))
  X = X.stack(vector(ZZ,fpowers[i].coefficients(sparse=False)))
  X = X.LLL(delta=0.3,early_red=True,use_siegel=True)
  print 'miniLLL time',t.stop().cputime,'bits',[[Mji.nbits() for Mji in Mi] for Mi in X]
  X = n * X

t.start()
for i in range(k,m):
  X = X.augment(vector(ZZ,i))
  X = X.stack(vector(ZZ,[0]*(i - k) + (H^(i-k)*fpowers[k]).coefficients(sparse=False)))

X = X.LLL(delta=0.6,early_red=True,use_siegel=True)
print 'bigLLL time',t.stop().cputime,'bits',[[Mji.nbits() for Mji in Mi] for Mi in X]
numLLL = 1

shift1 = matrix([[binomial(j,i) for i in range(m)] for j in range(m)])
shift2 = shift1*shift1

t.start()
factored = false
while True:
  # search u-Hv,...,u+Hv in steps of v
  # i.e. search u+vs with s being -H,...,H
  
  M0 = X[0]
  
  Q0 = sum(ZZ(z/H^i)*x^i for i,z in enumerate(M0))
  Qroots = Q0.roots(ring=ZZ)
  
  for r,multiplicity in Qroots:
    if u+v*r > 0:
      g = gcd(n,u+v*r)
      if g > 1 and g < n:
        if not factored:
          print '----- successful factorization',[g,n/g]
          factored = True
          # could abort at this point but want to benchmark failure case

  u += 2*H*v
  if u-H*v > pmax: break
  X *= shift2
  X = X.LLL(delta=0.6,early_red=True,use_siegel=True)
  numLLL += 1
print 'scan time',t.stop().cputime,'numLLL',numLLL

timetry = ttry.stop().cputime
print 'avgyearsx4',timetry * smooth / (86400 * 365.25)

