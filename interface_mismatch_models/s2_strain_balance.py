"""STAGE 2: the strain-balanced common in-plane lattice constant.
For a coherent interface both sides share one 'a'. Scan a; at each a relax c and all
internal coordinates (epitaxial constraint: a,b fixed, c free); the balanced a minimises
the TOTAL energy of the stack. 30 QLs each side, 5 atoms per QL -> equal atom counts,
so the stack energy is the plain average of the two per-atom energies."""
import json, numpy as np
from ase.io import read, write
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter
from calorine.calculators import CPUNEP
POT='nep.txt'

class EpiFilter(FrechetCellFilter):
    """relax c and internals, keep a and b frozen"""
    def __init__(self, at):
        super().__init__(at, mask=[False,False,True,False,False,False])

def relax_at_a(at0, a_target):
    at = at0.copy()
    c = np.array(at.get_cell()); s = at.get_scaled_positions()
    scale = a_target/np.linalg.norm(c[0])
    c[0]*=scale; c[1]*=scale
    at.set_cell(c, scale_atoms=False); at.set_scaled_positions(s)
    at.calc = CPUNEP(POT)
    BFGS(EpiFilter(at), logfile=None).run(fmax=1e-4, steps=400)
    return at

bi0 = read('s1_Bi2Te3_free.xyz'); sb0 = read('s1_Sb2Te3_free.xyz')
agrid = np.arange(4.22, 4.42, 0.01)
EB, ES = [], []
for a in agrid:
    b = relax_at_a(bi0, a); s = relax_at_a(sb0, a)
    EB.append(b.get_potential_energy()/len(b)); ES.append(s.get_potential_energy()/len(s))
EB, ES = np.array(EB), np.array(ES)
Etot = 0.5*(EB+ES)
# parabolic refinement around the minimum
i = Etot.argmin(); p = np.polyfit(agrid[i-2:i+3], Etot[i-2:i+3], 2)
a_bal = -p[1]/(2*p[0])
pB = np.polyfit(agrid, EB, 2); pS = np.polyfit(agrid, ES, 2)
a_Bi_free = -pB[1]/(2*pB[0]); a_Sb_free = -pS[1]/(2*pS[0])
print(f"free minima from the scan:  Bi2Te3 a = {a_Bi_free:.4f}   Sb2Te3 a = {a_Sb_free:.4f}")
print(f"misfit = {100*(a_Bi_free-a_Sb_free)/a_Sb_free:.2f} %")
print(f"\nSTRAIN-BALANCED a = {a_bal:.4f} A")
print(f"   strain on Bi2Te3 = {100*(a_bal-a_Bi_free)/a_Bi_free:+.2f} %")
print(f"   strain on Sb2Te3 = {100*(a_bal-a_Sb_free)/a_Sb_free:+.2f} %")

bi = relax_at_a(bi0, a_bal); sb = relax_at_a(sb0, a_bal)
res={}
for nm, at in (('Bi2Te3', bi), ('Sb2Te3', sb)):
    c = np.array(at.get_cell())
    s3 = np.sort(at.get_scaled_positions()[:,2])*c[2,2]
    d = np.diff(np.concatenate([s3,[s3[0]+c[2,2]]]))
    st = at.get_stress(voigt=True)/0.0062415091
    res[nm]=dict(a=float(a_bal), c=float(c[2,2]), cQL=float(c[2,2]/3),
                 gap=float(d.max()), spacings=sorted(set(np.round(d,4).tolist())),
                 E_per_atom=float(at.get_potential_energy()/len(at)),
                 sigma_zz_GPa=float(st[2]), max_force=float(np.abs(at.get_forces()).max()))
    write(f's2_{nm}_epitaxial.xyz', at)
    r=res[nm]
    print(f"\n{nm} at the balanced a:")
    print(f"   c_QL = {r['cQL']:.4f} A   gap = {r['gap']:.4f} A   spacings {r['spacings']}")
    print(f"   E = {r['E_per_atom']:.6f} eV/atom   sigma_zz = {r['sigma_zz_GPa']:.5f} GPa"
          f"   max|F| = {r['max_force']:.2e} eV/A")
res['a_balanced']=float(a_bal); res['a_Bi_free']=float(a_Bi_free); res['a_Sb_free']=float(a_Sb_free)
json.dump(res, open('s2_epitaxial.json','w'), indent=2)
