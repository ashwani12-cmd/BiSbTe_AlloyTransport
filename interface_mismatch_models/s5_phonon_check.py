"""STAGE 5: dynamical stability of the corrected epitaxial bulks + the DMM conductance
they imply. The published Bi2Te3 side carried 2.10 % imaginary modes; this must not."""
import sys, json, numpy as np
sys.path.insert(0,'../INTERFACE_DMM')
from ase.io import read
from ase import Atoms
from calorine.calculators import CPUNEP
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
POT='nep.txt'; H,KB=6.62607015e-34,1.380649e-23
NBIN=500; fgrid=np.linspace(0,5.2,NBIN); dfS=(fgrid[1]-fgrid[0])*1e12

def flux(at,dim,mesh):
    ph=Phonopy(PhonopyAtoms(symbols=at.get_chemical_symbols(),cell=at.get_cell().array,
               scaled_positions=at.get_scaled_positions()),
               supercell_matrix=np.diag(dim),primitive_matrix='auto')
    ph.generate_displacements(distance=0.01)
    F=[]
    for s in ph.supercells_with_displacements:
        a=Atoms(symbols=s.symbols,positions=s.positions,cell=s.cell,pbc=True)
        a.calc=CPUNEP(POT); F.append(a.get_forces())
    ph.forces=np.array(F); ph.produce_force_constants(); ph.symmetrize_force_constants()
    ph.run_mesh(mesh,with_group_velocities=True,is_mesh_symmetry=False)
    md=ph.get_mesh_dict(); f=md['frequencies']; gv=md['group_velocities']
    w=md['weights'].astype(float); V=float(abs(np.linalg.det(ph.primitive.cell)))*1e-30
    W=np.repeat(w[:,None],f.shape[1],axis=1); good=f>0.05
    idx=np.clip((f/(fgrid[1]-fgrid[0])).astype(int),0,NBIN-1)
    vz=gv[:,:,2]*100.0; m=good&(vz>0)
    Phi=np.bincount(idx[m],weights=vz[m]*W[m],minlength=NBIN)/(V*w.sum()*dfS)
    return Phi,float(f.min()),float(f.max()),int((f<-0.05).sum()),f.size

res={}
for nm,fn in (('Bi2Te3','s2_Bi2Te3_epitaxial.xyz'),('Sb2Te3','s2_Sb2Te3_epitaxial.xyz')):
    at=read(fn)
    Phi,fmin,fmax,nim,tot=flux(at,[4,4,1],[24,24,8])
    res[nm]=dict(Phi=Phi,fmin=fmin,fmax=fmax,nimag=nim,tot=tot)
    print(f"{nm} (corrected, epitaxial): f = {fmin:.3f} .. {fmax:.3f} THz   "
          f"imaginary {nim}/{tot} = {100*nim/tot:.3f} %")

P1,P2=res['Bi2Te3']['Phi'],res['Sb2Te3']['Phi']
den=P1+P2; al=np.where(den>0,P2/np.where(den>0,den,1),0.0)
def dndT(f,T):
    x=np.clip(H*f*1e12/(KB*T),1e-12,500); e=np.exp(x); return (H*f*1e12/(KB*T*T))*e/(e-1)**2
print(f"\n{'T(K)':>5} {'G_DMM corrected':>16} {'G_DMM published':>16} {'paper NEMD':>12}")
PUB={200:24.61,300:13.06,400:13.11,500:13.13}; PAP={200:34.75,300:56.94,400:60.76,500:98.60}
out={}
for T in (200,300,400,500):
    G=float(np.sum((H*fgrid*1e12)*P1*al*dndT(fgrid,T))*dfS)/1e6
    out[T]=G
    print(f"{T:5d} {G:16.2f} {PUB[T]:16.2f} {PAP[T]:12.2f}")
json.dump(dict(G_DMM={str(k):v for k,v in out.items()},
               imag={k:100*v['nimag']/v['tot'] for k,v in res.items()},
               fmax={k:v['fmax'] for k,v in res.items()}),
          open('s5_phonons.json','w'),indent=2)
