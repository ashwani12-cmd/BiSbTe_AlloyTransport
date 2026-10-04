"""STAGE 1: NEP-equilibrium R-3m bulk Bi2Te3 and Sb2Te3, then the epitaxially strained
versions at a common in-plane lattice constant.

R-3m: 15-atom hexagonal cell, 3 QLs, registry advancing +1 per plane across EVERY plane
including the van der Waals gaps (that continuity is what the published cell got wrong)."""
import json, numpy as np
from ase import Atoms
from ase.io import write
from ase.optimize import BFGS
from ase.filters import FrechetCellFilter, UnitCellFilter
from ase.constraints import FixSymmetry
from calorine.calculators import CPUNEP

POT='nep.txt'; SITE=np.array([1/3,2/3])

def hexcell(a,c): return np.array([[a,0,0],[-a/2,a*np.sqrt(3)/2,0],[0,0,c]])

def r3m(a, cQL, d1, d2, cat):
    """15-atom R-3m cell: 3 QLs, registry advances +1 on every plane, gaps included."""
    zq=np.array([0,d1,d1+d2,d1+2*d2,2*d1+2*d2]); c=3*cQL
    frac,syms=[],[]; m=0
    for n in range(3):
        for k in range(5):
            s=(m*SITE)%1.0
            frac.append([s[0],s[1],(n*cQL+zq[k])/c])
            syms.append('Te' if k%2==0 else cat); m+=1
    return Atoms(symbols=syms,scaled_positions=frac,cell=hexcell(a,c),pbc=True)

def measure(at,cat):
    """read back a,cQL,d1,d2 from a relaxed cell"""
    c=np.array(at.get_cell()); a=np.linalg.norm(c[0])
    s3=np.sort(at.get_scaled_positions()[:,2])*c[2,2]
    d=np.diff(np.concatenate([s3,[s3[0]+c[2,2]]]))
    gap=d.max(); cQL=c[2,2]/3
    dd=np.sort(d)[:4]   # the four smallest are the QL-internal spacings (per QL)
    return dict(a=float(a),c=float(c[2,2]),cQL=float(cQL),gap=float(gap),
                spacings=[float(x) for x in np.unique(np.round(d,4))])

out={}
# --- free relaxation: each material's own equilibrium
for cat,a0,cQL0,d10,d20 in [('Bi',4.38,10.16,1.75,2.03),('Sb',4.26,10.15,1.68,1.92)]:
    at=r3m(a0,cQL0,d10,d20,cat); at.calc=CPUNEP(POT)
    opt=BFGS(FrechetCellFilter(at),logfile=None); opt.run(fmax=1e-4,steps=500)
    nm=f"{'Bi2Te3' if cat=='Bi' else 'Sb2Te3'}"
    m=measure(at,cat); m['E_per_atom']=float(at.get_potential_energy()/len(at))
    m['max_force']=float(np.abs(at.get_forces()).max())
    m['stress_GPa']=[float(x) for x in at.get_stress(voigt=True)/0.0062415091]
    out[nm]=m; write(f's1_{nm}_free.xyz',at)
    print(f"{nm} free NEP equilibrium:")
    print(f"   a = {m['a']:.4f} A   c = {m['c']:.4f} A   c_QL = {m['cQL']:.4f}   gap = {m['gap']:.4f}")
    print(f"   E = {m['E_per_atom']:.6f} eV/atom   max|F| = {m['max_force']:.2e}   "
          f"|stress| = {np.abs(m['stress_GPa']).max():.4f} GPa")
    print(f"   distinct plane spacings: {m['spacings']}")
json.dump(out,open('s1_bulk_free.json','w'),indent=2)
