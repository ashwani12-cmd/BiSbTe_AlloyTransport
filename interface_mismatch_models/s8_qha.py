"""STAGE 8: NEP + calorine + phonopy quasi-harmonic results for the CORRECTED cell
at 200/300/400/500 K.

For each material we need F(a, c, T) = E_static(a,c) + F_vib(a,c,T), with internal
coordinates relaxed at every (a,c). The superlattice is coherent, so at each T the common
in-plane a is the one minimising the stack free energy; each side then takes its own c(T).
That gives thermal expansion, and hence a T-dependent geometry to feed the DMM -- the
honest analogue of what the NPT stage of the NEMD does.
"""
import json, numpy as np
from ase import Atoms
from ase.optimize import BFGS
from calorine.calculators import CPUNEP
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
POT='nep.txt'; SITE=np.array([1/3,2/3])
TS=[200,300,400,500]

def hexcell(a,c): return np.array([[a,0,0],[-a/2,a*np.sqrt(3)/2,0],[0,0,c]])
def r3m(a,cQL,d1,d2,cat):
    zq=np.array([0,d1,d1+d2,d1+2*d2,2*d1+2*d2]); c=3*cQL
    fr,sy=[],[]; m=0
    for n in range(3):
        for k in range(5):
            s=(m*SITE)%1.0; fr.append([s[0],s[1],(n*cQL+zq[k])/c])
            sy.append('Te' if k%2==0 else cat); m+=1
    return Atoms(symbols=sy,scaled_positions=fr,cell=hexcell(a,c),pbc=True)

def point(a,cQL,d1,d2,cat):
    """relax internals at fixed (a,c); return E_static, free energies F_vib(T), geometry"""
    at=r3m(a,cQL,d1,d2,cat); at.calc=CPUNEP(POT)
    BFGS(at,logfile=None).run(fmax=1e-4,steps=300)          # atoms only, cell fixed
    E=at.get_potential_energy()
    ph=Phonopy(PhonopyAtoms(symbols=at.get_chemical_symbols(),cell=at.get_cell().array,
               scaled_positions=at.get_scaled_positions()),
               supercell_matrix=np.diag([4,4,1]),primitive_matrix='auto')
    ph.generate_displacements(distance=0.01)
    F=[]
    for s in ph.supercells_with_displacements:
        b=Atoms(symbols=s.symbols,positions=s.positions,cell=s.cell,pbc=True)
        b.calc=CPUNEP(POT); F.append(b.get_forces())
    ph.forces=np.array(F); ph.produce_force_constants(); ph.symmetrize_force_constants()
    ph.run_mesh([16,16,6]); ph.run_thermal_properties(temperatures=TS)
    tp=ph.get_thermal_properties_dict()
    # phonopy free energy is kJ/mol per UNIT CELL -> eV per cell
    Fvib=np.array(tp['free_energy'])*0.0103642695/ (1.0)       # kJ/mol -> eV
    nim=int((ph.mesh.frequencies<-0.05).sum())
    return E, Fvib, len(at), nim, tp

g=json.load(open('s2_epitaxial.json'))
base={'Bi':(g['Bi2Te3']['cQL'],sorted(g['Bi2Te3']['spacings'])[0],sorted(g['Bi2Te3']['spacings'])[1]),
      'Sb':(g['Sb2Te3']['cQL'],sorted(g['Sb2Te3']['spacings'])[0],sorted(g['Sb2Te3']['spacings'])[1])}
agrid=np.array([4.27,4.30,4.3162,4.33,4.36])
store={}
for cat in ('Bi','Sb'):
    cQL0,d1,d2=base[cat]
    for a in agrid:
        for fc in (0.97,0.99,1.00,1.01,1.03,1.05):
            E,Fv,n,nim,_=point(a,cQL0*fc,d1,d2,cat)
            store[(cat,round(float(a),4),round(fc,3))]=dict(E=E,F=Fv.tolist(),n=n,
                                                            cQL=cQL0*fc,nimag=nim)
    print(f"  {cat}: {len(agrid)*6} (a,c) points done")
json.dump({f"{k[0]}|{k[1]}|{k[2]}":v for k,v in store.items()},open('s8_qha_grid.json','w'),indent=1)
print("wrote s8_qha_grid.json")
