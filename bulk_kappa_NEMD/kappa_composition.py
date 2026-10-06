#!/usr/bin/env python
"""
Composition-dependent bulk thermal conductivity of (Bi(1-x)Sb(x))2Te3 at 300 K, from the
NEMD runs in this folder.  Walks <comp>/<dir>/L_*/ , reports kappa per run, extrapolates
kappa_inf for each (composition, direction), and writes kappa_composition.csv.

Three kappa are reported per run, differing ONLY in how the heat flux and the geometry are
taken.  They are not alternatives to choose between after the fact -- the first two contain
known errors and are printed so the correction is auditable:

  k_pub  the recipe used for the published tables: the local virial flux jp divided by a
         bin volume computed as L*(1-2*0.02-2*0.18)/8, and the cross-section A_cross as
         printed in nemd_setup.txt.  Both are wrong (see README.md): the bin width is 4/3
         too large, and A_cross carries a 2/sqrt(3) factor for the hexagonal cell.
  k_jp   the same local-flux route with both geometry errors fixed: each bin's own width
         from the GROUP BOUNDARIES table, and the true cross-section recomputed from the
         NPT-equilibrated cell in thermo.out.
  k_th   the thermostat route: the flux is the energy actually added at the source and
         removed at the sink (the last two columns of compute.out, differentiated in
         time), divided by the true cross-section.  No bin volume enters at all, which is
         why this is the route to quote.

Needs only numpy.  compute.out and thermo.out may be plain or gzipped; numpy reads .gz
directly, so nothing has to be unpacked.

    python kappa_composition.py              # table to stdout + kappa_composition.csv
    python kappa_composition.py --no-csv
"""
import argparse, csv, os, re
import numpy as np

eV   = 1.602176634e-19
amu  = 1.66053906660e-27
JPc  = (eV**1.5)*(amu**-0.5)      # GPUMD's jp unit -> SI
NUM  = re.compile(r'-?\d+\.\d+')

COMPS = ['Bi2Te3', 'BiSbTe20', 'BiSbTe40', 'BiSbTe60', 'BiSbTe80', 'Sb2Te3']
XSB   = {'Bi2Te3': 0, 'BiSbTe20': 20, 'BiSbTe40': 40,
         'BiSbTe60': 60, 'BiSbTe80': 80, 'Sb2Te3': 100}
DIRS  = ['X', 'Y', 'Z']
LENS  = ['L_0250A', 'L_0500A', 'L_0750A', 'L_1000A']
FRACW, FRACS = 0.02, 0.18          # wall / source fractions of the published recipe

# kappa_perp(Bi2Te3) and kappa_par(Sb2Te3) as printed in the paper, for comparison only
PAPER = {'Bi2Te3': (1.37, 0.88), 'Sb2Te3': (2.02, 0.69)}

HERE = os.path.dirname(os.path.abspath(__file__))


def find(run, name):
    """path to <name> or <name>.gz inside run, or None."""
    for p in (os.path.join(run, name), os.path.join(run, name + '.gz')):
        if os.path.exists(p):
            return p
    return None


def parse_setup(run):
    """L_transport, the A_cross as printed, and the 8 group boundaries."""
    p = find(run, 'nemd_setup.txt')
    if p is None:
        return None
    L = A = None
    bins = {}
    txt = open(p, errors='replace').read().splitlines()
    for s in txt:
        if L is None and 'L_transport' in s:
            L = float(NUM.search(s).group())
        if A is None and 'A_cross' in s and 'A2' in s:
            A = float(NUM.search(s).group())
    intab = False
    for s in txt:
        if 'GROUP BOUNDARIES' in s:
            intab = True
            continue
        if not intab:
            continue
        if 'L_EFF' in s or s.strip().startswith('==='):
            break
        m = re.match(r'\s*(\d+)\s+(\S+)', s)
        if not m:
            continue
        g, role = int(m.group(1)), m.group(2)
        f = NUM.findall(s)
        if len(f) < 3 or 'wall' in role:
            continue
        if 1 <= g <= 8 and g not in bins:
            bins[g] = (float(f[0]), float(f[1]))
    if L is None or A is None or len(bins) < 8:
        return None
    return dict(L=L, A_used=A, bins=bins)


def true_geom(run, L_setup):
    """(A_true, L_actual) from the NPT-equilibrated cell vectors in thermo.out.

    The transport axis is identified as the cell length closest to L_transport, so this
    does not depend on how the axes happen to be ordered."""
    p = find(run, 'thermo.out')
    if p is None:
        return None, None
    th = np.loadtxt(p)
    if th.ndim == 1 or th.shape[1] < 18:
        return None, None
    t = th[-200:] if len(th) > 200 else th
    a, b, c = t[:, 9:12].mean(0), t[:, 12:15].mean(0), t[:, 15:18].mean(0)
    V = abs(np.dot(a, np.cross(b, c)))
    lens = np.array([np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)])
    i = int(np.argmin(abs(lens - L_setup)))
    return V/lens[i], lens[i]


def analyse(run, direction):
    s = parse_setup(run)
    if s is None:
        return None
    f = find(run, 'compute.out')
    if f is None:
        return None
    dat = np.loadtxt(f)
    if dat.ndim != 2 or len(dat) < 500:
        return None
    n = len(dat)
    st = n//2                                    # second half = steady state
    A_true, L_act = true_geom(run, s['L'])
    if A_true is None:                           # no thermo.out: fall back, flagged below
        A_true, L_act = s['A_used'], s['L']
    sc = L_act/s['L']                            # NPT rescaling of the bin grid
    fit = [2, 3, 4, 5]                           # G3..G6, away from both thermostats

    xc = np.array([(s['bins'][g][0] + s['bins'][g][1])/2*sc*1e-10 for g in range(1, 9)])
    T = dat[st:, 1:9].mean(0)
    grad = np.polyfit(xc[fit], T[fit], 1)[0]

    Lbin_true = np.mean([(s['bins'][g][1] - s['bins'][g][0]) for g in range(2, 8)])*sc*1e-10
    Lbin_pub = s['L']*(1 - 2*FRACW - 2*FRACS)/8*1e-10

    blk = {'X': 0, 'Y': 1, 'Z': 2}[direction]
    jp = dat[st:, 9 + blk*9 + 1: 9 + blk*9 + 9].mean(0)
    ajp = jp[fit].mean()
    Q_jp_true = ajp*JPc/((A_true*1e-20)*Lbin_true)
    Q_jp_pub = ajp*JPc/((s['A_used']*1e-20)*Lbin_pub)

    t = np.arange(n)*1e-12                       # compute every 1 ps
    Ps = np.polyfit(t[st:], dat[st:, -2]*eV, 1)[0]
    Pk = np.polyfit(t[st:], dat[st:, -1]*eV, 1)[0]
    imb = (abs(Pk) - abs(Ps))/((abs(Pk) + abs(Ps))/2)
    Q_th = 0.5*(abs(Ps) + abs(Pk))/(A_true*1e-20)

    return dict(nrows=n, L_nm=L_act/10, A_used=s['A_used'], A_true=A_true,
                grad=grad, imb=imb,
                k_pub=abs(Q_jp_pub/grad), k_jp=abs(Q_jp_true/grad), k_th=abs(Q_th/grad))


def fsc(Ls, ks):
    """kappa_inf = 1/intercept of 1/kappa vs 1/L, with the fit's R^2."""
    Ls, ks = np.asarray(Ls, float), np.asarray(ks, float)
    if len(Ls) < 2:
        return np.nan, np.nan
    c = np.polyfit(1/Ls, 1/ks, 1)
    r2 = 1 - np.sum((1/ks - np.polyval(c, 1/Ls))**2)/np.sum((1/ks - (1/ks).mean())**2)
    return 1/c[1], r2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--no-csv', action='store_true')
    args = ap.parse_args()

    DB, rows = {}, []
    print('='*104)
    print('PER-RUN kappa (W/m/K)   k_pub = published recipe | k_jp = geometry-corrected jp '
          '| k_th = thermostat')
    print('='*104)
    print(f"{'comp':<9}{'d':<2}{'len':<9}{'rows':>6}{'L_nm':>7}{'A_used':>9}{'A_true':>9}"
          f"{'k_pub':>8}{'k_jp':>8}{'k_th':>8}{'imbal':>8}")
    for c in COMPS:
        for d in DIRS:
            for L in LENS:
                run = os.path.join(HERE, c, d, L)
                r = analyse(run, d) if os.path.isdir(run) else None
                if r is None:
                    print(f"{c:<9}{d:<2}{L:<9}{'  -- not present / unreadable':<40}")
                    continue
                DB[(c, d, L)] = r
                rows.append(dict(composition=c, x_Sb_percent=XSB[c], direction=d, length=L,
                                 rows=r['nrows'], L_nm=round(r['L_nm'], 2),
                                 A_used_A2=round(r['A_used'], 1),
                                 A_true_A2=round(r['A_true'], 1),
                                 dTdx_K_per_m=f"{r['grad']:.6e}",
                                 kappa_pub=round(r['k_pub'], 4),
                                 kappa_jp=round(r['k_jp'], 4),
                                 kappa_th=round(r['k_th'], 4),
                                 source_sink_imbalance=round(r['imb'], 4)))
                print(f"{c:<9}{d:<2}{L:<9}{r['nrows']:>6}{r['L_nm']:>7.1f}{r['A_used']:>9.0f}"
                      f"{r['A_true']:>9.0f}{r['k_pub']:>8.3f}{r['k_jp']:>8.3f}"
                      f"{r['k_th']:>8.3f}{r['imb']:>+8.2f}")

    print()
    print('='*104)
    print('FINITE-SIZE EXTRAPOLATION   kappa_inf = 1/intercept of 1/kappa vs 1/L')
    print('='*104)
    print(f"{'comp':<9}{'xSb%':>5}{'dir':>4}{'n':>3}{'kinf_pub':>10}{'kinf_jp':>9}"
          f"{'kinf_th':>9}{'R2(jp)':>8}{'R2(th)':>8}")
    SUM = {}
    for c in COMPS:
        for d in DIRS:
            Ls = [DB[(c, d, L)]['L_nm'] for L in LENS if (c, d, L) in DB]
            if len(Ls) < 2:
                print(f"{c:<9}{XSB[c]:>5}{d:>4}{len(Ls):>3}   (insufficient lengths)")
                continue
            kp = [DB[(c, d, L)]['k_pub'] for L in LENS if (c, d, L) in DB]
            kj = [DB[(c, d, L)]['k_jp'] for L in LENS if (c, d, L) in DB]
            kt = [DB[(c, d, L)]['k_th'] for L in LENS if (c, d, L) in DB]
            ip, _ = fsc(Ls, kp)
            ij, r2j = fsc(Ls, kj)
            it, r2t = fsc(Ls, kt)
            SUM[(c, d)] = (ip, ij, it)
            print(f"{c:<9}{XSB[c]:>5}{d:>4}{len(Ls):>3}{ip:>10.3f}{ij:>9.3f}{it:>9.3f}"
                  f"{r2j:>8.3f}{r2t:>8.3f}")

    print()
    print('='*104)
    print('COMPOSITION DEPENDENCE AT 300 K   in-plane = (X+Y)/2, cross-plane = Z')
    print('='*104)
    print(f"{'comp':<9}{'xSb%':>5}{'par_pub':>9}{'par_th':>8}{'perp_pub':>10}{'perp_th':>9}"
          f"{'anis_pub':>10}{'anis_th':>9}{'paper_par':>11}{'paper_perp':>12}")
    series = []
    for c in COMPS:
        X, Y, Z = SUM.get((c, 'X')), SUM.get((c, 'Y')), SUM.get((c, 'Z'))
        if not (X and Y and Z):
            print(f"{c:<9}{XSB[c]:>5}   incomplete")
            continue
        par_pub, par_th = (X[0] + Y[0])/2, (X[2] + Y[2])/2
        perp_pub, perp_th = Z[0], Z[2]
        ref = PAPER.get(c, (np.nan, np.nan))
        series.append(dict(composition=c, x_Sb_percent=XSB[c],
                           kappa_par_pub=round(par_pub, 3), kappa_par_th=round(par_th, 3),
                           kappa_perp_pub=round(perp_pub, 3), kappa_perp_th=round(perp_th, 3),
                           anisotropy_pub=round(par_pub/perp_pub, 2),
                           anisotropy_th=round(par_th/perp_th, 2),
                           paper_par=ref[0], paper_perp=ref[1]))
        print(f"{c:<9}{XSB[c]:>5}{par_pub:>9.3f}{par_th:>8.3f}{perp_pub:>10.3f}"
              f"{perp_th:>9.3f}{par_pub/perp_pub:>10.2f}{par_th/perp_th:>9.2f}"
              f"{ref[0]:>11.2f}{ref[1]:>12.2f}")
    print("\nQuote the _th columns.  paper_* are the published values, shown for comparison "
          "only;\nthey were computed with the k_pub recipe and both of its geometry errors.")

    if not args.no_csv and rows:
        with open(os.path.join(HERE, 'kappa_composition.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        with open(os.path.join(HERE, 'kappa_vs_composition.csv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(series[0]))
            w.writeheader()
            w.writerows(series)
        print("\nwrote kappa_composition.csv (per run) and kappa_vs_composition.csv (series)")


if __name__ == '__main__':
    main()
