import sys, time
sys.path.insert(0, '/home/user/Solving-Birch-and-Swinnerton-Dye-Conjecture/research/kato-sha')
from kato_lib import *
from sage.all import EllipticCurve, Qp
for lab, ps in [('11a1',[3,7]), ('389a1',[3,5,7,11,13])]:
    E = EllipticCurve(lab)
    ms = PrimitiveModSym(E)
    print(lab, 'g =', ms.g, 'n_manin', ms.n_manin, 'build %.2f content %.2f'%(ms.t_build, ms.t_content))
    print('  sanity', ms.sanity_checks())
    for p in ps:
        ap = int(E.ap(p))
        for n in [2,3]:
            M = n+3
            co = ordinary_riemann(ms, p, ap, n, 4, M)
            vals = [vp_int(c, p) for c in co]
            precs = [ordinary_precision(p, n, j) for j in range(5)]
            L = E.padic_lseries(p)
            s = L.series(n, prec=5)
            # compare: Sage series coefficient * c_inf / g should equal ours
            cinf = E.real_components()
            K = Qp(p, M+5)
            ok = []
            for j in range(min(5, s.prec())):
                sj = s[j]
                ours = K(co[j])
                theirs = K(sj.lift() if hasattr(sj,'lift') else sj) * cinf / K(ms.g)
                pr = min(precs[j], sj.precision_absolute() + 0 if sj.precision_absolute() < 10**6 else precs[j], M)
                diff = (ours - theirs)
                ok.append(diff.valuation() >= min(precs[j], M) - 0 if diff != 0 else True)
            print(f'  p={p} n={n} ap={ap} coeffs mod p^{M}: {co} vals {vals} prec {precs[:5]}  sage: {s}  agree(mod prec): {ok}')
