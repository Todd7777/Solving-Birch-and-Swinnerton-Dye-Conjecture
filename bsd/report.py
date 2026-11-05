"""Markdown / JSON reporting."""
from __future__ import annotations

import json
from collections import Counter


def _fmt(x, nd=10):
    if x is None:
        return "—"
    if isinstance(x, float):
        return f"{x:.{nd}g}"
    return str(x)


def table_header() -> str:
    return ("| curve | N | r_an | rigorous ord | 2-descent rank bounds | #tors | ∏c_p | Ω_E | Reg | L^(r)(E,1)/r! | Ш_an | square? | status |\n"
            "|---|---|---|---|---|---|---|---|---|---|---|---|---|")


def table_row(rec: dict) -> str:
    inv, an, des = rec["invariants"], rec["analytic"], rec["descent"]
    rig = f"≥{an['rigorous_order_lower']}" if not an["rigorous_order_exact"] else f"= {an['analytic_rank']}"
    status = "theorem" if rec["rank_status"].startswith("THEOREM") else "open"
    sq = "" if rec["sha_is_square"] is None else ("yes" if rec["sha_is_square"] else "**NO**")
    err = " ⚠" if rec["errors"] else ""
    return (f"| {rec['label']} | {inv['conductor']} | {an['analytic_rank']} | {rig} | [{des['rank_lower']},{des['rank_upper']}] | "
            f"{inv['torsion_order']} | {inv['tamagawa_product']} | {_fmt(inv['omega'])} | {_fmt(des['regulator'])} | "
            f"{_fmt(an['leading_term'])} | {_fmt(rec['sha_an'], 8)} | {sq} | {status}{err} |")


def records_table(records: list[dict]) -> str:
    return "\n".join([table_header()] + [table_row(r) for r in records])


def curve_section(rec: dict) -> str:
    inv, an, des = rec["invariants"], rec["analytic"], rec["descent"]
    lines = [f"### {rec['label']}", ""]
    lines.append(f"Minimal model a-invariants `{inv['ainvs']}`, conductor N = {inv['conductor']}, "
                 f"Δ_min = {inv['discriminant']}, root number {inv['root_number']}, torsion {inv['torsion_structure'] or 'trivial'}"
                 + (f", CM by the order of discriminant {inv['cm_discriminant']}" if inv["cm_discriminant"] else "") + ".")
    if inv["local_data"]:
        ld = ", ".join(f"p={d['p']}: {d['reduction']}, Kodaira {d['kodaira']}, ord_p Δ = {d['ord_disc']}, c_p = {d['tamagawa']}"
                       for d in inv["local_data"])
        lines.append(f"Bad primes — {ld}.")
    lines.append(f"Rational isogeny degrees from E: {inv['isogeny_degrees']} (E[ℓ] reducible for ℓ in {inv['reducible_primes']}).")
    lines.append("")
    lines.append(f"**Analytic side.** ord_{{s=1}} L(E,s) = {an['analytic_rank']} numerically; "
                 f"L^({an['analytic_rank']})(E,1)/{an['analytic_rank']}! = {_fmt(an['leading_term'], 15)}. "
                 f"Exact L(E,1)/Ω⁺ = {an['L1_over_omega_plus_exact']} (modular symbols). {an['certification']}")
    lines.append("")
    lines.append(f"**Algebraic side.** 2-descent: {des['rank_lower']} ≤ rank ≤ {des['rank_upper']}"
                 + (f" (Ш[2]/2Ш[4] has 2-rank {des['sha2_mod_4_rank']})" if des['sha2_mod_4_rank'] not in (None, -1) else "")
                 + f". Saturated generators ({des['gens_source']}, saturated at primes < {des['saturation_bound']}): "
                 f"{des['gens'] or 'none'}; regulator {_fmt(des['regulator'], 15)}.")
    if des["notes"]:
        lines.append("Notes: " + "; ".join(des["notes"]))
    lines.append("")
    lines.append(f"**BSD quotient.** Ω_E = {_fmt(inv['omega'], 15)}, ∏c_p = {inv['tamagawa_product']}, #tors = {inv['torsion_order']}, "
                 f"so L^(r)(E,1)/r! ÷ (Ω·Reg·∏c_p/#tors²) = **{_fmt(rec['sha_an'], 12)}**"
                 + (f" → #Ш_an = {rec['sha_an_int']}" + (" (a perfect square, as Cassels' theorem requires)" if rec['sha_is_square'] else " **(NOT a square!)**") if rec["sha_an_int"] is not None else ""))
    lines.append("")
    lines.append(f"**Status.** {rec['rank_status']} Ш: {rec['sha_status']}.")
    if rec["p_part"]:
        proven = [str(p["p"]) for p in rec["p_part"] if p["status"] == "proven"]
        notc = [str(p["p"]) for p in rec["p_part"] if p["status"] == "not covered"]
        if proven or notc:
            lines.append(f"p-part of the BSD formula proven by the mechanised theorems at p ∈ {{{', '.join(proven)}}}"
                         + (f"; not covered by them at p ∈ {{{', '.join(notc)}}}" if notc else "") + ".")
    if rec.get("padic"):
        fin = [str(d["p"]) for d in rec["padic"] if d.get("kato", {}).get("sha_p_finite")]
        units = [f"p={d['p']}: {d['sha_p']}" for d in rec["padic"] if d.get("sha_p")]
        if fin:
            lines.append(f"**Ш[p^∞] is finite for p ∈ {{{', '.join(fin)}}}** (Kato 2004 + Mazur's control theorem + "
                         f"certified vanishing order of the p-adic L-function; cf. Stein–Wuthrich 2013).")
        if units:
            lines.append("p-adic BSD check (L_p^{(r)}(0)/r! ÷ p-adic regulator, conjecturally #Ш): " + "; ".join(units) + ".")
    if rec.get("sage"):
        sg = rec["sage"]
        if sg.get("saturation"):
            st = sg["saturation"]
            lines.append(f"Sage provable saturation: index of PARI basis = {st.get('index')}, provable height lower bound "
                         f"{st.get('height_lower_bound')}: {st.get('status')}.")
        if sg.get("prove_bsd"):
            pb = sg["prove_bsd"]
            lines.append(f"**Sage `prove_BSD`:** {pb['status']}." + (f" Proof sketch: `{pb['proof_sketch']}`" if pb.get("proof_sketch") else ""))
        if sg.get("error"):
            lines.append(f"Sage stage error: {sg['error']}")
    if rec.get("cremona"):
        c = rec["cremona"]
        lines.append(f"Cross-check with Cremona's tables (independent software): rank {c['rank']}, Reg {c['regulator']}, Ш_an {c['sha_an']}"
                     + (" — **mismatch: " + "; ".join(c["mismatches"]) + "**" if c["mismatches"] else " — agree."))
    if rec["errors"]:
        lines.append("⚠ " + " | ".join(rec["errors"]))
    lines.append("")
    return "\n".join(lines)


def summary(records: list[dict]) -> str:
    n = len(records)
    ranks = Counter(r["analytic"]["analytic_rank"] for r in records)
    shas = Counter(r["sha_an_int"] for r in records if r["sha_an_int"] is not None)
    nonsq = [r["label"] for r in records if r["sha_is_square"] is False]
    errs = [r["label"] for r in records if r["errors"]]
    mism = [r["label"] for r in records if r.get("cremona") and r["cremona"]["mismatches"]]
    lines = [f"- curves processed: {n}",
             "- analytic rank distribution: " + ", ".join(f"rank {k}: {v}" for k, v in sorted(ranks.items())),
             "- Ш_an distribution: " + ", ".join(f"{k}: {v}" for k, v in sorted(shas.items(), key=lambda kv: (kv[0] is None, kv[0]))),
             f"- Ш_an not a perfect square: {len(nonsq)}" + (f" ({', '.join(nonsq[:20])})" if nonsq else ""),
             f"- curves with errors/warnings: {len(errs)}" + (f" ({', '.join(errs[:20])}{'…' if len(errs) > 20 else ''})" if errs else ""),
             f"- mismatches with Cremona's tables: {len(mism)}" + (f" ({', '.join(mism[:20])})" if mism else "")]
    return "\n".join(lines)


def write_jsonl(path: str, records: list[dict]):
    with open(path, "w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")


def read_jsonl(path: str) -> list[dict]:
    out = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out
