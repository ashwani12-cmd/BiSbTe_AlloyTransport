"""STAGE 19: the full set of analytic interface models on the CORRECTED cell.

s9 gave DMM only.  This adds the two things needed to bracket the NEMD result:

  * **phonon radiation limit** -- every incident phonon transmitted (alpha = 1).
    A hard upper bound on any harmonic model.  Computed from BOTH sides; the smaller
    one is the bound that matters.
  * **AMM** (acoustic mismatch model) -- specular transmission set by the acoustic
    impedances Z = rho*v, with the proper angle integration over the transmission cone
    and Snell refraction, per polarisation, using sound speeds ALONG z (the transport
    direction) taken from the slope of the NEP acoustic branches at Gamma.

Everything is evaluated on the SAME quasi-harmonically expanded, strain-balanced
geometries s9 used, so DMM / AMM / radiation limit are mutually consistent.

Honest about the approximation levels:
  DMM and the radiation limit use the FULL NEP dispersion.
  AMM's transmission is an elastic-continuum (long-wavelength) quantity; applying it
  across the whole spectrum is AMM's standard weakness, not a bug here.  It is applied
  as a flux-weighted average transmission times the radiation limit, and the
  per-polarisation numbers are printed so the spread is visible.
"""
import json, numpy as np
from ase import Atoms
from ase.optimize import BFGS
from calorine.calculators import CPUNEP
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

POT = 'nep.txt'; SITE = np.array([1/3, 2/3]); TS = [200, 300, 400, 500]
H, KB, AMU = 6.62607015e-34, 1.380649e-23, 1.66053906660e-27
NBIN = 500; fgrid = np.linspace(0, 5.2, NBIN); dfS = (fgrid[1] - fgrid[0]) * 1e12
MASS = {'Bi': 208.98040, 'Sb': 121.760, 'Te': 127.60}

raw = json.load(open('s8_qha_grid.json'))
G = {}
for k, v in raw.items():
    cat, a, fc = k.split('|'); G[(cat, float(a), float(fc))] = v
agrid = sorted(set(k[1] for k in G)); fcs = sorted(set(k[2] for k in G))

def Fmin_c(cat, a, ti):
    cs = np.array([G[(cat, a, fc)]['cQL'] for fc in fcs])
    Fs = np.array([G[(cat, a, fc)]['E'] + G[(cat, a, fc)]['F'][ti] for fc in fcs]) / 15.0
    p = np.polyfit(cs, Fs, 2); c = -p[1] / (2 * p[0])
    return float(np.polyval(p, c)), float(c)

def hexcell(a, c): return np.array([[a, 0, 0], [-a/2, a*np.sqrt(3)/2, 0], [0, 0, c]])

def r3m(a, cQL, d1, d2, cat):
    zq = np.array([0, d1, d1+d2, d1+2*d2, 2*d1+2*d2]); c = 3*cQL
    fr, sy = [], []; m = 0
    for n in range(3):
        for k in range(5):
            s = (m*SITE) % 1.0
            fr.append([s[0], s[1], (n*cQL+zq[k])/c]); sy.append('Te' if k % 2 == 0 else cat); m += 1
    return Atoms(symbols=sy, scaled_positions=fr, cell=hexcell(a, c), pbc=True)

def phonons(at):
    ph = Phonopy(PhonopyAtoms(symbols=at.get_chemical_symbols(), cell=at.get_cell().array,
                 scaled_positions=at.get_scaled_positions()),
                 supercell_matrix=np.diag([4, 4, 1]), primitive_matrix='auto')
    ph.generate_displacements(distance=0.01)
    F = []
    for s in ph.supercells_with_displacements:
        b = Atoms(symbols=s.symbols, positions=s.positions, cell=s.cell, pbc=True)
        b.calc = CPUNEP(POT); F.append(b.get_forces())
    ph.forces = np.array(F); ph.produce_force_constants(); ph.symmetrize_force_constants()
    return ph

def flux(ph):
    """one-sided phonon flux spectrum Phi(f) [1/(m^2 s) per Hz], full dispersion"""
    ph.run_mesh([24, 24, 8], with_group_velocities=True, is_mesh_symmetry=False)
    md = ph.get_mesh_dict(); f = md['frequencies']; gv = md['group_velocities']
    w = md['weights'].astype(float)
    V = float(abs(np.linalg.det(ph.primitive.cell))) * 1e-30
    W = np.repeat(w[:, None], f.shape[1], axis=1); good = f > 0.05
    idx = np.clip((f / (fgrid[1] - fgrid[0])).astype(int), 0, NBIN - 1)
    vz = gv[:, :, 2] * 100.0                      # THz*Angstrom -> m/s
    m = good & (vz > 0)
    return (np.bincount(idx[m], weights=vz[m] * W[m], minlength=NBIN) / (V * w.sum() * dfS),
            float(f.max()), int((f < -0.05).sum()), f.size)

def sound_speeds_z(ph, c_hex_ang):
    """v along CARTESIAN z for the 3 acoustic branches, labelled LA/TA by EIGENVECTOR.

    Two traps, both hit on the way to this version:
    1. primitive_matrix='auto' puts phonopy in the 5-atom RHOMBOHEDRAL cell, so [0,0,q]
       is not along z and its period is not the hexagonal c.  q must be built from the
       Cartesian k:  q_j = (a_j . k)/2pi over the PRIMITIVE vectors a_j, and then
       v = domega/dk = 2*pi*df/dk.  Getting this wrong gave 16-26 km/s.
    2. Sorting the three acoustic branches by frequency does NOT identify LA.  In Sb2Te3
       the degenerate TA pair lies ABOVE the LA singleton, so "largest = LA" labels it
       backwards and AMM then pairs LA on one side with TA on the other.  Polarisation is
       taken from the eigenvector: LA if the displacement is mostly along z.
    """
    prim = np.array(ph.primitive.cell)            # Angstrom, rows = lattice vectors
    az = prim[:, 2] * 1e-10                       # z-components, m
    d = c_hex_ang / 3 * 1e-10                     # one QL repeat, m
    ks = np.linspace(0.0, 0.03 * 2 * np.pi / d, 8)[1:]
    qs = [list(az * k / (2 * np.pi)) for k in ks]
    ph.run_qpoints(qs, with_eigenvectors=True)
    qd = ph.get_qpoints_dict()
    fr = qd['frequencies']; ev = qd['eigenvectors']
    nb = fr.shape[1]

    # follow the 3 lowest branches at the SMALLEST k and classify them there
    i0 = 0
    order0 = np.argsort(fr[i0])[:3]
    pol = []
    for b in order0:
        e = ev[i0][:, b].reshape(-1, 3)
        pz = float(np.sum(np.abs(e[:, 2]) ** 2))
        pxy = float(np.sum(np.abs(e[:, :2]) ** 2))
        pol.append('LA' if pz > pxy else 'TA')
    if pol.count('LA') != 1:
        # degenerate mixing at tiny k: fall back to the most-z-polarised of the three
        zf = []
        for b in order0:
            e = ev[i0][:, b].reshape(-1, 3)
            zf.append(np.sum(np.abs(e[:, 2]) ** 2) /
                      max(np.sum(np.abs(e) ** 2), 1e-30))
        pol = ['TA'] * 3; pol[int(np.argmax(zf))] = 'LA'

    ac = np.sort(fr, axis=1)[:, :3]
    v, r2 = [], []
    for b in range(3):
        y = ac[:, b] * 1e12
        p = np.polyfit(ks, y, 1)
        v.append(2 * np.pi * abs(p[0]))
        yy = np.polyval(p, ks)
        r2.append(1 - np.sum((y - yy) ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-30))
    # branch b of the sorted set corresponds to order0[b] at the smallest k
    lab = {}
    ta = [v[i] for i in range(3) if pol[i] == 'TA']
    la = [v[i] for i in range(3) if pol[i] == 'LA'][0]
    lab = dict(vTA1=float(min(ta)), vTA2=float(max(ta)), vLA=float(la),
               r2_min=float(min(r2)), pol_order=pol)
    return lab

def amm_transmission(v1, v2, Z1, Z2, n=4001):
    """angle-averaged AMM energy transmission 1->2 for one polarisation.
    alpha(th1) = 4 Z1 Z2 cos1 cos2 / (Z1 cos1 + Z2 cos2)^2 ; Snell sin2 = (v2/v1) sin1.
    Normalised by int cos sin dth so that alpha==1 everywhere gives 1."""
    thc = np.arcsin(min(1.0, v1 / v2)) if v2 > v1 else np.pi / 2
    th = np.linspace(0, thc, n)
    s2 = np.clip((v2 / v1) * np.sin(th), -1, 1)
    c1, c2 = np.cos(th), np.cos(np.arcsin(s2))
    al = 4 * Z1 * Z2 * c1 * c2 / (Z1 * c1 + Z2 * c2) ** 2
    num = np.trapezoid(al * np.cos(th) * np.sin(th), th)
    den = np.trapezoid(np.cos(th) * np.sin(th), np.linspace(0, np.pi / 2, n))
    return float(num / den), float(np.degrees(thc))

def dT_lambda_over_delta(tau):
    """Modified-Landauer interface temperature drop, as a fraction of the reservoir
    difference Delta = T_e,1 - T_e,2.  Shi et al., as used by Chowdhury et al.
    (ACS Appl. Mater. Interfaces 13, 4636 (2021)) for THIS interface:

        T_lambda,1 = T_e,1 - tau*Delta/2
        T_lambda,2 = T_e,2 + (1-tau)*Delta/2
        => dT_lambda = Delta/2,  independent of tau

    The ORIGINAL Landauer/DMM form divides the same heat flux by the full Delta, so
    G_modified = G_original / 0.5 = 2 x G_original exactly.  The modified form is the
    standard fix for the original's pathology of a finite resistance at an IMAGINARY
    interface inside a bulk.  Written out rather than hard-coding "x2" so that the
    assumption is visible and breaks loudly if the convention is ever revised.
    """
    T1 = -tau / 2.0                       # T_lambda,1 - T_e,1, in units of Delta
    T2 = (1.0 - tau) / 2.0                # T_lambda,2 - T_e,2, in units of Delta
    return 1.0 + T1 - T2                  # (T_lambda,1 - T_lambda,2)/Delta


def dndT(f, T):
    x = np.clip(H * f * 1e12 / (KB * T), 1e-12, 500); e = np.exp(x)
    return (H * f * 1e12 / (KB * T * T)) * e / (e - 1) ** 2

g0 = json.load(open('s2_epitaxial.json'))
d1B, d2B = sorted(g0['Bi2Te3']['spacings'])[:2]
d1S, d2S = sorted(g0['Sb2Te3']['spacings'])[:2]

NEMD = {200: 58.09, 300: 57.75, 400: 69.63, 500: 75.57}      # seed 1, kapitza.py
PAPER = {200: 34.75, 300: 56.94, 400: 60.76, 500: 98.60}

out = []
print(f"{'T':>4} {'G_DMM':>8} {'DMM_mod':>8} {'G_AMM':>8} {'G_rad,Bi':>9} {'G_rad,Sb':>9} "
      f"{'a_DMM':>7} {'a_AMM':>7} {'NEMD':>7} {'paper':>7}")
print("-" * 92)
for ti, T in enumerate(TS):
    Ft = []
    for a in agrid:
        fB, _ = Fmin_c('Bi', a, ti); fS, _ = Fmin_c('Sb', a, ti); Ft.append(0.5 * (fB + fS))
    p = np.polyfit(agrid, Ft, 2); aT = -p[1] / (2 * p[0])
    cB = np.polyval(np.polyfit(agrid, [Fmin_c('Bi', a, ti)[1] for a in agrid], 2), aT)
    cS = np.polyval(np.polyfit(agrid, [Fmin_c('Sb', a, ti)[1] for a in agrid], 2), aT)

    rec = {'T': T, 'a': float(aT), 'cQL_Bi': float(cB), 'cQL_Sb': float(cS)}
    phs, fl, sp, rho = {}, {}, {}, {}
    for tag, cQL, d1, d2, cat in (('Bi', cB, d1B, d2B, 'Bi'), ('Sb', cS, d1S, d2S, 'Sb')):
        at = r3m(aT, cQL, d1, d2, cat); at.calc = CPUNEP(POT)
        BFGS(at, logfile=None).run(fmax=1e-4, steps=300)
        cc = float(np.array(at.get_cell())[2, 2])
        ph = phonons(at); phs[tag] = ph
        fl[tag] = flux(ph)[0]
        sp[tag] = sound_speeds_z(ph, cc)
        V = float(abs(np.linalg.det(np.array(at.get_cell())))) * 1e-30
        rho[tag] = sum(MASS[s] for s in at.get_chemical_symbols()) * AMU / V

    # ---- DMM: detailed-balance transmission from the two flux spectra ----
    P1, P2 = fl['Bi'], fl['Sb']
    den = P1 + P2
    aD = np.where(den > 0, P2 / np.where(den > 0, den, 1), 0.0)
    hw = H * fgrid * 1e12
    wgt = hw * P1 * dndT(fgrid, T)
    G_DMM = float(np.sum(wgt * aD) * dfS) / 1e6
    G_radB = float(np.sum(hw * P1 * dndT(fgrid, T)) * dfS) / 1e6
    G_radS = float(np.sum(hw * P2 * dndT(fgrid, T)) * dfS) / 1e6
    aD_eff = G_DMM / G_radB
    # modified Landauer (the convention Chowdhury et al. use for this interface)
    fr = [dT_lambda_over_delta(t) for t in (0.0, 0.25, 0.5, 0.75, 1.0)]
    assert max(abs(x - 0.5) for x in fr) < 1e-12, f"dT_lambda/Delta is not 0.5: {fr}"
    G_DMM_mod = G_DMM / fr[0]

    # ---- AMM: per-polarisation impedance transmission, angle averaged ----
    amm = {}
    for pol in ('vTA1', 'vTA2', 'vLA'):
        v1, v2 = sp['Bi'][pol], sp['Sb'][pol]
        assert 800 < v1 < 6000 and 800 < v2 < 6000, \
            f"sound speed out of physical range: {pol} {v1:.0f} / {v2:.0f} m/s"
        Z1, Z2 = rho['Bi'] * v1, rho['Sb'] * v2
        g, thc = amm_transmission(v1, v2, Z1, Z2)
        amm[pol] = dict(v_Bi=v1, v_Sb=v2, Z_Bi=Z1, Z_Sb=Z2, gamma=g, theta_c_deg=thc)
    aA = float(np.mean([amm[p]['gamma'] for p in ('vTA1', 'vTA2', 'vLA')]))
    G_AMM = aA * G_radB

    print(f"{T:4d} {G_DMM:8.2f} {G_DMM_mod:8.2f} {G_AMM:8.2f} {G_radB:9.2f} {G_radS:9.2f} "
          f"{aD_eff:7.3f} {aA:7.3f} {NEMD[T]:7.2f} {PAPER[T]:7.2f}")
    rec.update(G_DMM=G_DMM, G_DMM_modified_Landauer=G_DMM_mod, G_AMM=G_AMM, G_rad_Bi=G_radB, G_rad_Sb=G_radS,
               alpha_DMM_eff=aD_eff, alpha_AMM=aA, AMM_per_pol=amm,
               rho_Bi=rho['Bi'], rho_Sb=rho['Sb'],
               sound_Bi=sp['Bi'], sound_Sb=sp['Sb'],
               NEMD_seed1=NEMD[T], paper=PAPER[T])
    out.append(rec)

json.dump(out, open('s19_models.json', 'w'), indent=2)
print("\n--- sound speeds along z (m/s), 200 K ---")
r = out[0]
for pol in ('vTA1', 'vTA2', 'vLA'):
    d = r['AMM_per_pol'][pol]
    print(f"  {pol:5s} Bi2Te3 {d['v_Bi']:7.1f}  Sb2Te3 {d['v_Sb']:7.1f}   "
          f"Z_Bi {d['Z_Bi']/1e6:6.3f}  Z_Sb {d['Z_Sb']/1e6:6.3f} MRayl   "
          f"gamma {d['gamma']:.4f}  theta_c {d['theta_c_deg']:.1f} deg")
print(f"  rho: Bi2Te3 {r['rho_Bi']:.1f}  Sb2Te3 {r['rho_Sb']:.1f} kg/m3")
for tag in ('Bi','Sb'):
    print(f"  {tag}2Te3 branch labels (lowest 3 at smallest k): {r['sound_'+tag]['pol_order']}"
          f"  -> TA {r['sound_'+tag]['vTA1']:.1f}/{r['sound_'+tag]['vTA2']:.1f}, "
          f"LA {r['sound_'+tag]['vLA']:.1f} m/s")
print(f"  acoustic-branch linearity: R2_min Bi {r['sound_Bi']['r2_min']:.6f}  "
      f"Sb {r['sound_Sb']['r2_min']:.6f}  (must be ~1 or the slope is not a sound speed)")
print("\n  cross-check vs elastic constants: v = sqrt(C/rho)")
print(f"    LA_z needs C33 = rho*v^2 = {r['rho_Bi']*r['AMM_per_pol']['vLA']['v_Bi']**2/1e9:.1f} GPa (Bi2Te3)"
      f", {r['rho_Sb']*r['AMM_per_pol']['vLA']['v_Sb']**2/1e9:.1f} GPa (Sb2Te3)")
print(f"    TA_z needs C44 = {r['rho_Bi']*r['AMM_per_pol']['vTA1']['v_Bi']**2/1e9:.1f} GPa (Bi2Te3)"
      f", {r['rho_Sb']*r['AMM_per_pol']['vTA1']['v_Sb']**2/1e9:.1f} GPa (Sb2Te3)")
print("\n  G_DMM        = original Landauer, referenced to the RESERVOIR temperature "
      "difference")
print("  G_DMM_mod    = modified Landauer (Shi et al.), referenced to the local modal "
      "temperatures")
print("                 next to the interface -- exactly 2x, and the convention "
      "Chowdhury et al.")
print("                 (ACS AMI 13, 4636 (2021)) use for this same Bi2Te3-Sb2Te3 "
      "interface, where")
print("                 they report ~75 MW/m2K saturating above 40 K.")
print("\nwrote s19_models.json")
