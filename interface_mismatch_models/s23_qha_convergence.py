"""STAGE 23: is the refined QHA grid actually converged?

The whole point of s22 is that the s8 grid gave a thermal expansion that depended on the
F(c) fit order by a factor of two.  A denser, centred grid is only an improvement if the
fitted minimum has stopped caring about the fit.  This script is the gate.

Three tests, all of which must pass before s22_qha_grid.json is used for anything:

  1. FIT-ORDER INDEPENDENCE.  c(T) from polynomial fits of order 2, 3, 4 (and a local
     3-point parabola through the sampled minimum, which uses no global fit at all).
     Delta c(200->500 K) must agree across these to a few per cent, not a factor of two.
  2. WINDOW INDEPENDENCE.  Refit using only the points within +-1 %, +-1.5 % and +-2 % of
     the minimum.  A converged fit does not care where the window is cut.
  3. RESIDUALS.  The quadratic residual must be small on this grid (it was 0.40 meV/atom
     on the old one).
  4. DYNAMICAL-STABILITY CONTAMINATION.  Some grid points carry imaginary modes -- this
     is EXPECTED, because the grid deliberately samples compressed and expanded cells that
     are not all stable, and it was never checked before.  Diagnosed: they are entirely in
     the lowest acoustic branch at LARGE |q| (0.43-0.62, zone-boundary, NOT a near-Gamma
     artefact), only in Bi2Te3, and they appear at the compressed edge a = 4.27 and grow
     steeply with expansion.  On the OLD grid fc = 1.05 carried 677 imaginary modes at all
     five a-points -- an unstable, over-expanded cell whose F_vib is meaningless, sitting
     inside the fit.  That is a second, independent reason the old quadratic was biased.
     The requirement is therefore NOT "zero imaginary modes anywhere" (that would reject
     any grid wide enough to bracket the minimum) but:
        (a) points near the minimum must be essentially clean, and
        (b) the fitted Delta c must not move when every contaminated point is DISCARDED.

Prints the old grid alongside for comparison.
"""
import argparse, json, sys, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument('--amin', type=float, default=0.0,
                help='drop in-plane columns with a <= AMIN.  The production setting used by '
                     's24_models_FINAL.py is --amin 4.28, which removes the dynamically '
                     'unstable a = 4.27 column (0 of 17 clean c points).  Default 0 keeps '
                     'every column, which is the harder test and FAILS test 4 for the '
                     'reason documented in s24_FINAL_NOTES.md.')
AMIN = ap.parse_args().amin

def load(fn):
    raw = json.load(open(fn)); G = {}
    for k, v in raw.items():
        cat, a, fc = k.split('|')
        if float(a) <= AMIN: continue
        G[(cat, float(a), float(fc))] = v
    return G, sorted(set(k[1] for k in G)), sorted(set(k[2] for k in G))

def cmin(G, fcs, cat, a, ti, order=2, window=None):
    """minimising c of F(c) at fixed (a, T).  window = fractional half-width about the
    SAMPLED minimum; None uses every point."""
    cs = np.array([G[(cat, a, fc)]['cQL'] for fc in fcs])
    Fs = np.array([G[(cat, a, fc)]['E'] + G[(cat, a, fc)]['F'][ti] for fc in fcs]) / 15.0
    if window is not None:
        c0 = cs[np.argmin(Fs)]
        m = np.abs(cs/c0 - 1.0) <= window
        if m.sum() < order + 1: return np.nan, np.nan
        cs, Fs = cs[m], Fs[m]
    if order == 0:                      # local 3-point parabola, no global fit
        i = int(np.argmin(Fs))
        if i == 0 or i == len(Fs)-1: return float(cs[i]), np.nan
        x, y = cs[i-1:i+2], Fs[i-1:i+2]
        p = np.polyfit(x, y, 2)
        return float(-p[1]/(2*p[0])), float(np.polyval(p, -p[1]/(2*p[0])))
    p = np.polyfit(cs, Fs, order)
    if order == 2:
        c = -p[1]/(2*p[0])
    else:
        r = np.roots(np.polyder(p)); r = r[np.isreal(r)].real
        r = r[(r > cs.min()) & (r < cs.max())]
        c = float(r[np.argmin(np.polyval(p, r))]) if len(r) else float(cs[np.argmin(Fs)])
    return float(c), float(np.polyval(p, c))

def cT(G, agrid, fcs, cat, ti, order=2, window=None):
    """c(T) for one compound, at the a that minimises the MEAN free energy of the stack."""
    Ft, cB = [], []
    for a in agrid:
        fB = cmin(G, fcs, 'Bi', a, ti, order, window)[1]
        fS = cmin(G, fcs, 'Sb', a, ti, order, window)[1]
        Ft.append(0.5*(fB + fS))
        cB.append(cmin(G, fcs, cat, a, ti, order, window)[0])
    p = np.polyfit(agrid, Ft, 2); aT = -p[1]/(2*p[0])
    return float(np.polyval(np.polyfit(agrid, cB, 2), aT)), float(aT)

PASS = True
for fn, lab in (('s8_qha_grid.json', 'OLD grid (s8)'),
                ('s22_qha_grid.json', 'REFINED grid (s22)')):
    try:
        G, agrid, fcs = load(fn)
    except FileNotFoundError:
        print(f"\n=== {lab}: not present, skipped ==="); continue
    cs = [G[('Bi', agrid[len(agrid)//2], fc)]['cQL'] for fc in fcs]
    print(f"\n=== {lab}: a = {agrid}" + (f" (a<={AMIN} dropped)" if AMIN else " (all columns)")
          + f", {len(fcs)} c points, span {100*(max(cs)-min(cs))/np.mean(cs):.2f} % ===")

    print("  TEST 1 -- fit-order independence, Delta c_QL(Bi2Te3, 200->500 K):")
    dcs = {}
    for order, nm in ((0, 'local 3-pt'), (2, 'quadratic'), (3, 'cubic'), (4, 'quartic')):
        try:
            c0 = cT(G, agrid, fcs, 'Bi', 0, order)[0]; c1 = cT(G, agrid, fcs, 'Bi', 3, order)[0]
            dcs[nm] = c1 - c0
            print(f"    {nm:11s}: c {c0:.4f} -> {c1:.4f}   Delta = {c1-c0:+.4f} A")
        except Exception as e:
            print(f"    {nm:11s}: failed ({e})")
    ref = [v for k, v in dcs.items() if k != 'local 3-pt']
    spread = (max(ref) - min(ref))/abs(np.mean(ref)) if ref else 9.9
    ok1 = spread < 0.15
    print(f"    -> orders 2-4 spread = {100*spread:.1f} % of the mean  "
          f"[{'PASS' if ok1 else 'FAIL'}, need < 15 %]")

    print("  TEST 2 -- window independence (cubic), Delta c_QL:")
    ws = {}
    for w in (0.010, 0.015, 0.020, None):
        try:
            c0 = cT(G, agrid, fcs, 'Bi', 0, 3, w)[0]; c1 = cT(G, agrid, fcs, 'Bi', 3, 3, w)[0]
            ws[w] = c1 - c0
            print(f"    +-{'all ' if w is None else f'{100*w:.1f}%'}: Delta = {c1-c0:+.4f} A")
        except Exception as e:
            print(f"    {w}: failed ({e})")
    vv = [v for v in ws.values() if np.isfinite(v)]
    spread2 = (max(vv) - min(vv))/abs(np.mean(vv)) if len(vv) > 1 else 9.9
    ok2 = spread2 < 0.20
    print(f"    -> window spread = {100*spread2:.1f} %  [{'PASS' if ok2 else 'FAIL'}, need < 20 %]")

    print("  TEST 3 -- quadratic fit residual (300 K, Bi2Te3, mid a):")
    a0 = agrid[len(agrid)//2]
    Fs = np.array([G[('Bi', a0, fc)]['E'] + G[('Bi', a0, fc)]['F'][1] for fc in fcs])/15.0
    cc = np.array([G[('Bi', a0, fc)]['cQL'] for fc in fcs])
    r2 = Fs - np.polyval(np.polyfit(cc, Fs, 2), cc)
    ok3 = np.abs(r2).max()*1e3 < 0.10
    print(f"    max|residual| = {np.abs(r2).max()*1e3:.4f} meV/atom  "
          f"[{'PASS' if ok3 else 'FAIL'}, need < 0.10]")

    print("  TEST 4 -- contamination by dynamically unstable grid points:")
    nim = sum(v['nimag'] for v in G.values())
    NMODE = 16*16*6*15
    worst = max(v['nimag'] for v in G.values())
    # how contaminated is the neighbourhood of the minimum?
    near = [v['nimag'] for k, v in G.items()
            if abs(G[k]['cQL']/G[(k[0], k[1], min(fcs, key=lambda f: abs(f-1.005)))]['cQL'] - 1) < 0.01]
    print(f"    total imaginary modes on the grid: {nim}; worst single point "
          f"{worst} of {NMODE} ({100*worst/NMODE:.3f} %)")
    print(f"    worst point within +-1 % of the minimum: {max(near)} "
          f"({100*max(near)/NMODE:.3f} %)")
    ok4a = max(near)/NMODE < 1e-3
    # refit using ONLY points that carry no imaginary mode at all
    def cT_clean(cat, ti, order=3):
        Ft, cB = [], []
        for a in agrid:
            fcl = [fc for fc in fcs if G[('Bi', a, fc)]['nimag'] == 0
                                    and G[('Sb', a, fc)]['nimag'] == 0]
            if len(fcl) < order + 1: return np.nan
            fB = cmin(G, fcl, 'Bi', a, ti, order); fS = cmin(G, fcl, 'Sb', a, ti, order)
            Ft.append(0.5*(fB[1] + fS[1])); cB.append(cmin(G, fcl, cat, a, ti, order)[0])
        p = np.polyfit(agrid, Ft, 2); aT = -p[1]/(2*p[0])
        return float(np.polyval(np.polyfit(agrid, cB, 2), aT))
    try:
        dc_clean = cT_clean('Bi', 3) - cT_clean('Bi', 0)
        dc_all = cT(G, agrid, fcs, 'Bi', 3, 3)[0] - cT(G, agrid, fcs, 'Bi', 0, 3)[0]
        rel = abs(dc_clean - dc_all)/abs(dc_all)
        print(f"    Delta c, all points      : {dc_all:+.4f} A")
        print(f"    Delta c, CLEAN points only: {dc_clean:+.4f} A   "
              f"(shift {100*rel:.1f} %)")
        ok4b = rel < 0.10
    except Exception as e:
        print(f"    clean-only refit failed: {e}"); ok4b = False; rel = 9.9
    ok4 = ok4a and ok4b
    print(f"    -> [{'PASS' if ok4 else 'FAIL'}] near-minimum clean (<0.1 %) and "
          f"clean-only refit within 10 %")
    if 'REFINED' in lab:
        PASS = ok1 and ok2 and ok3 and ok4

print(f"\n{'='*64}\nREFINED GRID {'IS USABLE' if PASS else 'IS NOT USABLE -- do not proceed'}\n{'='*64}")
sys.exit(0 if PASS else 1)
