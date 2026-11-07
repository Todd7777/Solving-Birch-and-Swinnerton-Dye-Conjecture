import sys, time
sys.path.insert(0, '/home/user/Solving-Birch-and-Swinnerton-Dye-Conjecture/research/kato-sha')
from kato_lib import *
from sage.all import EllipticCurve, Qp, prime_range
cls = load_classes()
print(len(cls), cls[0], cls[-1]['label'])
ss = []
for c in cls[:60]:
    E = EllipticCurve(c['ainvs'])
    for p in prime_range(3, 30):
        if c['N'] % p and E.ap(p) == 0:
            ss.append((c['label'], p))
print(ss[:30])
