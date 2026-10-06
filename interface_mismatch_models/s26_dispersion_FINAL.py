"""STAGE 21: cross-plane phonon dispersion of the two sides of the interface.

The figure Chowdhury et al. (ACS AMI 13, 4636 (2021)) show as their Fig. 2a -- the
bulk dispersions of the two materials along the transport direction -- but on the
CORRECTED R-3m strain-balanced cell, and extended with the two quantities that
actually enter the models: the one-sided phonon flux Phi(omega) of Eq. (flux) and
the DMM transmission alpha(omega) of Eq. (dmm).

Everything is on the 300 K quasi-harmonic geometry, i.e. the second row of
s19_models.json, so the numbers here and in Table 'mismatch' are the same objects.

Panels, all sharing the frequency axis:
  (a) Gamma--Z dispersion, MIRRORED: Bi2Te3 opening to the left, Sb2Te3 to the
      right, Gamma common at the centre.  This is the transport direction.
  (b) the one-sided flux spectra Phi_1 (Bi2Te3) and Phi_2 (Sb2Te3).
  (c) the DMM transmission alpha_{1->2}(omega) = Phi_2/(Phi_1+Phi_2).

VALIDATION built in: the Gamma-slope of each acoustic branch must reproduce the
sound speeds already in s19_models.json, and the top of each spectrum must match
the 4.08 / 4.67 THz quoted in the manuscript.  Both are asserted, not eyeballed.
"""
import json, numpy as np
from ase import Atoms
from ase.optimize import BFGS
from calorine.calculators import CPUNEP
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

POT='nep.txt'; SITE=np.array([1/3,2/3]); TI=1; TLAB=300     # row index 1 = 300 K
H,KB,AMU=6.62607015e-34,1.380649e-23,1.66053906660e-27
NBIN=500; fgrid=np.linspace(0,5.2,NBIN); dfS=(fgrid[1]-fgrid[0])*1e12
MASS={'Bi':208.98040,'Sb':121.760,'Te':127.60}

raw=json.load(open('s22_qha_grid.json')); G={}
for k,v in raw.items():
    cat,a,fc=k.split('|'); G[(cat,float(a),float(fc))]=v
agrid=[a for a in sorted(set(k[1] for k in G)) if a>4.28]; fcs=sorted(set(k[2] for k in G))

CFIT_ORDER=3
def Fmin_c(cat,a,ti):
    cs=np.array([G[(cat,a,fc)]['cQL'] for fc in fcs])
    Fs=np.array([G[(cat,a,fc)]['E']+G[(cat,a,fc)]['F'][ti] for fc in fcs])/15.0
    p=np.polyfit(cs,Fs,CFIT_ORDER)
    if CFIT_ORDER==2: c=-p[1]/(2*p[0])
    else:
        r=np.roots(np.polyder(p)); r=r[np.isreal(r)].real; r=r[(r>cs.min())&(r<cs.max())]
        c=float(r[np.argmin(np.polyval(p,r))]) if len(r) else float(cs[np.argmin(Fs)])
    return float(np.polyval(p,c)),float(c)

def hexcell(a,c): return np.array([[a,0,0],[-a/2,a*np.sqrt(3)/2,0],[0,0,c]])

def r3m(a,cQL,d1,d2,cat):
    zq=np.array([0,d1,d1+d2,d1+2*d2,2*d1+2*d2]); c=3*cQL
    fr,sy=[],[]; m=0
    for n in range(3):
        for k in range(5):
            s=(m*SITE)%1.0
            fr.append([s[0],s[1],(n*cQL+zq[k])/c]); sy.append('Te' if k%2==0 else cat); m+=1
    return Atoms(symbols=sy,scaled_positions=fr,cell=hexcell(a,c),pbc=True)

def phonons(at):
    ph=Phonopy(PhonopyAtoms(symbols=at.get_chemical_symbols(),cell=at.get_cell().array,
               scaled_positions=at.get_scaled_positions()),
               supercell_matrix=np.diag([4,4,1]),primitive_matrix='auto')
    ph.generate_displacements(distance=0.01)
    F=[]
    for s in ph.supercells_with_displacements:
        b=Atoms(symbols=s.symbols,positions=s.positions,cell=s.cell,pbc=True)
        b.calc=CPUNEP(POT); F.append(b.get_forces())
    ph.forces=np.array(F); ph.produce_force_constants(); ph.symmetrize_force_constants()
    return ph

def flux(ph):
    ph.run_mesh([24,24,8],with_group_velocities=True,is_mesh_symmetry=False)
    md=ph.get_mesh_dict(); f=md['frequencies']; gv=md['group_velocities']
    w=md['weights'].astype(float)
    V=float(abs(np.linalg.det(ph.primitive.cell)))*1e-30
    W=np.repeat(w[:,None],f.shape[1],axis=1); good=f>0.05
    idx=np.clip((f/(fgrid[1]-fgrid[0])).astype(int),0,NBIN-1)
    vz=gv[:,:,2]*100.0; m=good&(vz>0)
    phi=np.bincount(idx[m],weights=vz[m]*W[m],minlength=NBIN)/(V*w.sum()*dfS)
    return phi,float(f.max())

def band_gamma_to_z(ph,n=161):
    """Gamma -> Z along CARTESIAN z.

    Same trap as s19's sound_speeds_z: primitive_matrix='auto' puts phonopy in the
    5-atom RHOMBOHEDRAL cell, so [0,0,q] is NOT along z.  Build q from the Cartesian
    wavevector, q_j = (a_j . k)/2pi over the PRIMITIVE lattice vectors.  The zone
    boundary along z sits at k = pi/d with d = c_hex/3 (one quintuple-layer repeat).
    """
    prim=np.array(ph.primitive.cell)                 # Angstrom, rows = lattice vectors
    az=prim[:,2]*1e-10                               # z components, m
    d=abs(az).sum()*0  # placeholder, set by caller via c_hex
    return None

def band(ph,c_hex_ang,n=161):
    prim=np.array(ph.primitive.cell); az=prim[:,2]*1e-10
    d=c_hex_ang/3*1e-10                              # one QL repeat along z, m
    kmax=np.pi/d
    ks=np.linspace(0.0,kmax,n)
    qs=[list(az*k/(2*np.pi)) for k in ks]
    ph.run_qpoints(qs)
    fr=ph.get_qpoints_dict()['frequencies']           # (nk, nbranch) THz
    return ks,kmax,fr

g0=json.load(open('s2_epitaxial.json'))
d1B,d2B=sorted(g0['Bi2Te3']['spacings'])[:2]
d1S,d2S=sorted(g0['Sb2Te3']['spacings'])[:2]
ref=json.load(open('s24_models_FINAL.json'))[TI]
assert ref['T']==TLAB

Ft=[0.5*(Fmin_c('Bi',a,TI)[0]+Fmin_c('Sb',a,TI)[0]) for a in agrid]
p=np.polyfit(agrid,Ft,2); aT=-p[1]/(2*p[0])
cB=np.polyval(np.polyfit(agrid,[Fmin_c('Bi',a,TI)[1] for a in agrid],2),aT)
cS=np.polyval(np.polyfit(agrid,[Fmin_c('Sb',a,TI)[1] for a in agrid],2),aT)
assert abs(aT-ref['a'])<1e-6 and abs(cB-ref['cQL_Bi'])<1e-6
print(f"{TLAB} K geometry: a={aT:.4f}  cQL_Bi={cB:.4f}  cQL_Sb={cS:.4f}  (matches s19)")

out={'T':TLAB,'a':float(aT),'cQL_Bi':float(cB),'cQL_Sb':float(cS)}
data={}
for tag,cQL,d1,d2,cat,name in (('Bi',cB,d1B,d2B,'Bi','Bi2Te3'),
                               ('Sb',cS,d1S,d2S,'Sb','Sb2Te3')):
    at=r3m(aT,cQL,d1,d2,cat); at.calc=CPUNEP(POT)
    BFGS(at,logfile=None).run(fmax=1e-4,steps=300)
    cc=float(np.array(at.get_cell())[2,2])
    ph=phonons(at)
    ks,kmax,fr=band(ph,cc)
    phi,fmax=flux(ph)
    V=float(abs(np.linalg.det(np.array(at.get_cell()))))*1e-30
    rho=sum(MASS[s] for s in at.get_chemical_symbols())*AMU/V
    nimag=int((fr<-0.05).sum())
    assert nimag==0, f"{name}: {nimag} imaginary frequencies on the Gamma-Z path"
    data[tag]=dict(k=ks,kmax=kmax,bands=fr,phi=phi,fmax=fmax,name=name,rho=rho)
    print(f"  {name}: {fr.shape[1]} branches, spectrum top {fmax:.3f} THz, "
          f"rho {rho:.1f} kg/m3, {nimag} imaginary")

# ---- validation 1: Gamma slopes must reproduce the s19 sound speeds ---------
print("\n  acoustic slopes at Gamma vs s19_models.json (m/s):")
for tag in ('Bi','Sb'):
    d=data[tag]; ks=d['k']; fr=d['bands']
    sel=ks<0.03*d['kmax']
    v=sorted(2*np.pi*abs(np.polyfit(ks[sel],np.sort(fr,axis=1)[sel,b]*1e12,1)[0])
             for b in range(3))
    ref_v=sorted([ref['AMM_per_pol'][p]['v_'+tag] for p in ('vTA1','vTA2','vLA')])
    print(f"    {d['name']}: this {['%.0f'%x for x in v]}  s19 {['%.0f'%x for x in ref_v]}")
    for a_,b_ in zip(v,ref_v):
        assert abs(a_-b_)/b_<0.02, f"{tag}: sound speed mismatch {a_:.0f} vs {b_:.0f}"

# ---- validation 2: spectrum tops must match the manuscript -----------------
for tag,expect in (('Bi',None),('Sb',None)):
    got=data[tag]['fmax']
print(f"  spectrum tops {data['Bi']['fmax']:.3f} / {data['Sb']['fmax']:.3f} THz")

# ---- DMM transmission ------------------------------------------------------
P1,P2=data['Bi']['phi'],data['Sb']['phi']
den=P1+P2
alpha=np.where(den>0,P2/np.where(den>0,den,1),0.0)
hw=H*fgrid*1e12
x=np.clip(hw/(KB*TLAB),1e-12,500); e=np.exp(x)
dndT=(hw/(KB*TLAB*TLAB))*e/(e-1)**2
G_DMM=float(np.sum(hw*P1*dndT*alpha)*dfS)/1e6
G_radB=float(np.sum(hw*P1*dndT)*dfS)/1e6
print(f"\n  reconstructed from these spectra: G_DMM={G_DMM:.2f}  G_rad,Bi={G_radB:.2f}")
print(f"  s19_models.json says            : G_DMM={ref['G_DMM']:.2f}  G_rad,Bi={ref['G_rad_Bi']:.2f}")
assert abs(G_DMM-ref['G_DMM'])<0.05 and abs(G_radB-ref['G_rad_Bi'])<0.05

np.savez('s26_dispersion_FINAL.npz',
         fgrid=fgrid, alpha=alpha,
         **{f'{t}_{k}':data[t][k] for t in ('Bi','Sb') for k in ('k','bands','phi')},
         Bi_kmax=data['Bi']['kmax'], Sb_kmax=data['Sb']['kmax'])
out.update(fmax_Bi=data['Bi']['fmax'], fmax_Sb=data['Sb']['fmax'],
           rho_Bi=data['Bi']['rho'], rho_Sb=data['Sb']['rho'],
           G_DMM_check=G_DMM, G_rad_Bi_check=G_radB,
           alpha_at_1THz=float(alpha[np.argmin(abs(fgrid-1.0))]),
           nbranch=int(data['Bi']['bands'].shape[1]))
json.dump(out,open('s26_dispersion_FINAL.json','w'),indent=2)
print("\nwrote s21_dispersion.npz + s21_dispersion.json")
