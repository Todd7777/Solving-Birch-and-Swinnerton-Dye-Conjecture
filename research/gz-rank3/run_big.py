"""Gross-Zagier certificate attempt for the curated rank-5 and rank-7 curves (large conductors).

Modular symbols at these levels are infeasible (PARI already overflows a 2 GB stack near N ~ 1.2e4),
so L(E^D,1) != 0 is established by a rigorous ball enclosure of the twisted L-series.  Everything is
resumable: finished steps and partial sums (exact dyadic balls) live in big_state.json; a_n are cached
as .npy files in a scratch directory (recomputed if missing).

Usage: python3 research/gz-rank3/run_big.py [name ...]     (names: rank5, rank7)
"""
import json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
import numpy as np                                                    # noqa: E402
from bsd import gz_certificate as gzc                                 # noqa: E402
from bsd.pari_session import get_pari                                 # noqa: E402

STATE = os.path.join(HERE, "big_state.json")
OUT = os.path.join(HERE, "big_certificates.jsonl")
CACHE = os.environ.get("GZ_CACHE", "/tmp/claude-0/-home-user/49d47a5d-6794-5fd0-9790-0bd7cf4fd119/scratchpad/gz")
os.makedirs(CACHE, exist_ok=True)
pari = get_pari()
pari.allocatemem(2 * 2**30, silent=True)

CURVES = {
    # smallest known conductor of rank 5 (Buhler-Gross-Zagier / Brumer-McGuinness type example); D = -7: twist has 2-Selmer rank 0
    "rank5": dict(label="19047851a1", ainvs=[0, 0, 1, -79, 342], D=-7, t_lp=-40.0, t_tw=-20.0, third=True),
    # smallest known conductor of rank 7; D = -23 is useless (E^-23 has rank 2, L(E^-23,1) = 0), D = -55: twist rank 0
    "rank7": dict(label="rank7-382623908456", ainvs=[0, 0, 0, -10012, 346900], D=-55, t_lp=-30.0, t_tw=-3.0, third=False),
}


def load_state():
    if os.path.exists(STATE):
        return json.load(open(STATE))
    return {}


def save_state(st):
    tmp = STATE + ".tmp"
    with open(tmp, "w") as f:
        json.dump(st, f, indent=1)
    os.replace(tmp, STATE)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def run(name):
    cfg = CURVES[name]
    st = load_state()
    S = st.setdefault(name, {})
    cd = gzc.curve_data(cfg["ainvs"])
    E = cd.pop("_E")
    N, D = cd["conductor"], cfg["D"]
    assert cd["root_number"] == -1
    assert gzc.is_heegner_discriminant(D, N), D
    S["curve"] = cd
    Et = pari.ellinit(pari.elltwist(E, D))
    S["twist_root_number"] = int(pari.ellrootno(Et))
    if "twist_2descent" not in S:
        r = pari.ellrank(pari.ellminimalmodel(Et))
        S["twist_2descent"] = dict(rank_lower=int(r[0]), rank_upper=int(r[1]),
                                   note="PARI ellrank on E^D (predictor only: rank 0 makes L(E^D,1) != 0 expected)")
    alpha = 2 * math.pi / math.sqrt(N)
    M1 = gzc._choose_terms(alpha, cfg["t_lp"], 1)
    MD = gzc._choose_terms(2 * math.pi / (abs(D) * math.sqrt(N)), cfg["t_tw"], 0)
    M3 = gzc._choose_terms(alpha, -30.0, 3) if cfg["third"] else 0
    Mmax = max(M1, MD, M3)
    S.update(M_lprime=M1, M_twist=MD, M_third=M3)
    save_state(st)
    log(name, "N", N, "D", D, "M1", M1, "MD", MD, "M3", M3)

    npy = os.path.join(CACHE, f"an_{name}_{Mmax}.npy")
    if os.path.exists(npy):
        an = np.load(npy, mmap_mode="r")
        log("loaded a_n cache", npy)
    else:
        t = time.time()
        an = gzc.an_numpy(E, Mmax, "int32", progress=lambda p: log("  sieve at p =", p))
        chk = gzc.an_list(E, 20000)
        assert [int(x) for x in an[:20000]] == chk, "sieve disagrees with ellan"
        np.save(npy, an)
        S["an_seconds"] = round(time.time() - t, 1)
        save_state(st)
        log("a_n computed", S["an_seconds"], "s")

    def ckpt(key):
        def f(n, state):
            st2 = load_state()
            st2[name][key + "_partial"] = dict(n=n, S=state, time=time.time())
            save_state(st2)
            log(f"  {key} checkpoint n={n}")
        return f

    # twisted value L(E^D,1)
    if "twist" not in S:
        part = S.get("twist_partial")
        res = gzc.twist_value_ball(an, N, D, cfg["t_tw"], M=MD, start_n=part["n"] if part else 1,
                                   S0=part["S"] if part else None, checkpoint=ckpt("twist"))
        st = load_state(); S = st[name]
        S["twist"] = dict(gzc._strip(res), state=gzc.arb_to_state(res["value"]))
        S.pop("twist_partial", None)
        save_state(st)
        log("L(E^D,1) =", res["ball"]["ball"], "excludes 0:", not res["ball"]["contains_zero"])
    # L'(E,1)
    if "lprime" not in S:
        part = S.get("lprime_partial")
        res = gzc.lprime_ball(an, N, cfg["t_lp"], M=M1, start_n=part["n"] if part else 1,
                              S0=part["S"] if part else None, checkpoint=ckpt("lprime"))
        st = load_state(); S = st[name]
        S["lprime"] = dict(gzc._strip(res), state=gzc.arb_to_state(res["value"]))
        S.pop("lprime_partial", None)
        save_state(st)
        log("L'(E,1) =", res["ball"]["ball"])
    if cfg["third"] and "third" not in S:
        res = gzc.third_derivative_ball(an, N, -30.0, M=M3)
        st = load_state(); S = st[name]
        S["third"] = dict(gzc._strip(res), state=gzc.arb_to_state(res["value"]))
        save_state(st)
        log("L'''(E,1)/3! =", res["ball"]["ball"])
    # volume, height bound, Manin constant, inequality
    vol = gzc.period_area_ball(cd["ainvs"])
    if "height" not in S:
        S["height"] = gzc.height_lower_bound_sage(cd["ainvs"])
        save_state(st)
    mc = gzc.manin_constant_status(cd)
    with gzc._prec(256):
        lp_v = gzc.arb_from_state(S["lprime"]["state"])
        tw_v = gzc.arb_from_state(S["twist"]["state"])
        gz = gzc.gz_inequality(D, lp_v, tw_v, vol["value"], S["height"]["value"], mc["c_max"])
    nonzero = not tw_v.contains(0)
    cond = []
    if mc["c_max"] is None:
        cond.append("Manin constant c_E <= 10^%.1f (Manin's conjecture: c_E = 1)" % gz["largest_admissible_manin_constant_log10"])
    if not nonzero:
        cond.append("L(E^D,1) != 0 not established")
    ok = gz["inequality_holds"] and nonzero and not cond
    cert = dict(label=cfg["label"], curve=cd, heegner=dict(D=D, class_number=int(pari.qfbclassno(D)), u_K=1,
                twist_root_number=S["twist_root_number"], twist_2descent=S["twist_2descent"]),
                lprime={k: v for k, v in S["lprime"].items() if k != "state"},
                twist={k: v for k, v in S["twist"].items() if k != "state"},
                third_derivative={k: v for k, v in S.get("third", {}).items() if k != "state"} or None,
                volume=gzc._strip(vol), height_lower_bound=S["height"], manin_constant=mc, gross_zagier=gz,
                twist_nonvanishing=dict(ok=nonzero, how="rigorous ball enclosure of L(E^D,1) (modular symbols infeasible)"),
                conclusion=dict(certified=bool(ok), ord_lower_bound=3 if (gz["inequality_holds"] and nonzero) else 1,
                                conditional_on=cond,
                                statement=("L'(E,1) = 0 exactly, hence ord >= 3" if ok else
                                           ("L'(E,1) = 0 and ord >= 3, CONDITIONAL on: " + "; ".join(cond))
                                           if gz["inequality_holds"] and nonzero else "not certified")))
    with open(OUT, "a") as f:
        f.write(json.dumps(cert, default=str) + "\n")
    log(name, cert["conclusion"]["statement"], "margin 10^%.1f" % gz["log10_margin"])


if __name__ == "__main__":
    for nm in (sys.argv[1:] or ["rank5", "rank7"]):
        run(nm)
