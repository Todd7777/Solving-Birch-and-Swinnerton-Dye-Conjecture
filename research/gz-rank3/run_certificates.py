"""Run the Gross-Zagier rank-3 certificate on every rank-3 curve of conductor < 20000 (Cremona).

Resumable: one JSON line per (label, D) in certificates.jsonl; finished items are skipped.
Usage:  python3 research/gz-rank3/run_certificates.py
"""
import json, os, sys, time, traceback
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
from bsd.ecdata import load_range                      # noqa: E402
from bsd.gz_certificate import certify_rank3, heegner_discriminants  # noqa: E402
from bsd.pari_session import get_pari                  # noqa: E402

OUT = os.path.join(HERE, "certificates.jsonl")
# PARI 2.17.2 msinit/msfromell overflow a 2 GB stack at these levels (observed in the first pass)
PARI_MODSYM_OVERFLOW = {"12279a1", "16811a1", "18562c1", "18745a1"}
pari = get_pari()
pari.allocatemem(2 * 2**30, silent=True)


def done_keys():
    keys = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            try:
                r = json.loads(line)
                keys.add((r["label"], r["heegner"]["D"]))
            except Exception:
                pass
    return keys


def rank3_curves():
    cs = [c for c in load_range(1, 19999) if c.rank == 3]
    return cs


def extra_checks(c):
    """mwrank rank proof + Sage CS bound + ellweilcurve Smith invariants (consistency of c_E = 1)."""
    import sage.all as S
    E = S.EllipticCurve(S.QQ, c.ainvs)
    out = {}
    t = time.time()
    out["mwrank_rank_proof"] = int(E.rank(proof=True, only_use_mwrank=True))
    out["mwrank_seconds"] = round(time.time() - t, 2)
    sat = E.saturation([E(P) for P in [(S.QQ(str(x)), S.QQ(str(y))) for x, y in c.gens]])
    out["cremona_gens_saturation_index"] = int(sat[1])
    Ep = pari.ellinit(c.ainvs)
    try:
        vE, vS = pari.ellweilcurve(Ep)
        out["ellweilcurve_smith"] = [[str(a) for a in s] for s in vS]
    except Exception as exc:  # noqa: BLE001
        out["ellweilcurve_smith"] = f"failed: {exc}"
    return out


def main():
    done = done_keys()
    for c in rank3_curves():
        N = c.conductor
        Ds = heegner_discriminants(N, 12)
        gens = [[str(x), str(y)] for x, y in c.gens]
        # primary certificate: smallest odd D with S_D != 0; secondary: next such D
        if (c.label, "primary") in done:
            continue
        prim_done = [k for k in done if k[0] == c.label]
        if len(prim_done) >= 2:
            continue
        print(f"=== {c.label} N={N} Heegner D: {Ds}", flush=True)
        try:
            mode = "eclib" if c.label in PARI_MODSYM_OVERFLOW else True
            cert = certify_rank3(c.ainvs, None, height_bound="sage", gens=gens, label=c.label, modsym=mode)
            ms = cert.get("modular_symbols", {})
            nonzero = [d["D"] for d in ms.get("exact", []) if d["S_D"] != "0"]
            ec = ms.get("eclib_cross_check") or {}
            if not nonzero and "values" in ec:
                nonzero = [int(d) for d, v in ec["values"].items() if v != "0"]
            cert["role"] = "primary"
            cert["extra"] = extra_checks(c)
            cert["cremona"] = dict(rank=c.rank, torsion=c.torsion_order, tamagawa_product=c.tamagawa_product,
                                   omega=c.omega, L_value=c.L_value, regulator=c.regulator, sha_an=c.sha_an)
            cert.pop("references", None)
            if (c.label, cert["heegner"]["D"]) not in done:
                with open(OUT, "a") as f:
                    f.write(json.dumps(cert, default=str) + "\n")
            print(c.label, cert["conclusion"]["statement"][:200], "margin 10^%.1f" % cert["gross_zagier"]["log10_margin"], flush=True)
            if len(nonzero) >= 2 and (c.label, nonzero[1]) not in done:
                cert2 = certify_rank3(c.ainvs, nonzero[1], height_bound="sage", gens=gens, label=c.label,
                                      want_third=False, modsym=mode)
                cert2["role"] = "secondary (independent D)"
                cert2.pop("references", None)
                with open(OUT, "a") as f:
                    f.write(json.dumps(cert2, default=str) + "\n")
                print(c.label, "second D", nonzero[1], cert2["conclusion"]["statement"][:120],
                      "margin 10^%.1f" % cert2["gross_zagier"]["log10_margin"], flush=True)
        except Exception:
            traceback.print_exc()


if __name__ == "__main__":
    main()
