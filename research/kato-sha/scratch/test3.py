import sys, time
sys.path.insert(0, '/home/user/Solving-Birch-and-Swinnerton-Dye-Conjecture/research/kato-sha')
from kato_lib import *
from sage.all import EllipticCurve, Qp, prime_range
for lab, p in [('1483a1',5), ('1613a1',5), ('794a1', 11), ('1531a1', 3)]:
    E = EllipticCurve(lab)
    ms = PrimitiveModSym(E)
    print(lab, 'p', p, 'ap', E.ap(p), 'g', ms.g, ms.sanity_checks())
    t1 = theta_X_coeffs(ms, p, 1)
    T1 = X_to_T(t1, 4)
    print('  theta1 T-coeffs', T1, 'vals', [vp_int(x,p) for x in T1], '=> L^- == -theta1 mod p for j<=p-1')
    t2 = theta_X_coeffs(ms, p, 2)
    Q2, rem2 = divide_by_cyclotomic(t2, p, 1)
    print('  theta2/Phi_p remainder zero?', all(x == 0 for x in rem2))
    G2 = X_to_T(Q2, 4)
    print('  G2 T-coeffs', G2, 'vals', [vp_int(x,p) for x in G2], '=> L^+ == G2 mod p')
    if p <= 11:
        t3 = theta_X_coeffs(ms, p, 3)
        Q3, rem3 = divide_by_cyclotomic(t3, p, 2)
        print('  theta3/Phi_p^2 remainder zero?', all(x == 0 for x in rem3))
        G3 = X_to_T(Q3, 4)
        print('  G3 T-coeffs', G3, 'mod p^2:', [x % p**2 for x in G3], '; -theta1 mod p:', [(-x) % p for x in T1], ' G3 mod p:', [x % p for x in G3])
    # Sage supersingular series for comparison: p * coeff_2(L_alpha) == L^-_2 + alpha L^+_2 (when lower coeffs of L^+- vanish)
    if p <= 7:
        L = E.padic_lseries(p)
        t=time.time()
        s = L.series(5 if p==5 else 4, prec=4)
        print('  sage L_alpha series', s, '%.1fs'%(time.time()-t))
        c2 = s[2]
        K = c2.parent()
        alpha = K.gen()
        assert alpha**2 - E.ap(p)*alpha + p == 0
        pred = (K(-T1[2]) + alpha*K(G2[2]))/K(p) * K(ms.g) / E.real_components()
        print('  sage c2 =', c2, '\n  predicted (L^-_2 + alpha L^+_2)/p * g / c_inf =', pred.add_bigoh(c2.precision_absolute()) , ' diff', (c2 - pred))
