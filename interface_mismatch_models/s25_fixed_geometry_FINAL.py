"""CONTROL for s19: are the harmonic models falling with T because of the QHA
expansion, or would they fall anyway?

The literature (Swartz-Pohl; Duda 2011; Reddy 2005; Li 2022) is unanimous that a
harmonic AMM/DMM at FIXED geometry rises as T^3 and then SATURATES -- it never
decreases.  Our Table 'mismatch' decreases by 6-7% from 200 to 500 K, and the text
attributes that to quasi-harmonic softening.  That is a causal claim, so test it:
freeze the geometry at the 200 K QHA solution, let ONLY dn/dT vary, and re-integrate.

Expected if the text is right: flat-to-rising.  If it still falls, the stated
mechanism is wrong.
"""
import json, numpy as np, sys
sys.path.insert(0, '.')
from ase import Atoms
from ase.optimize import BFGS
from calorine.calculators import CPUNEP
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

POT='nep.txt'; SITE=np.array([1/3,2/3]); TS=[200,300,400,500]
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
    return np.bincount(idx[m],weights=vz[m]*W[m],minlength=NBIN)/(V*w.sum()*dfS)

def sound_speeds_z(ph,c_hex_ang):
    prim=np.array(ph.primitive.cell); az=prim[:,2]*1e-10
    d=c_hex_ang/3*1e-10
    ks=np.linspace(0.0,0.03*2*np.pi/d,8)[1:]
    ph.run_qpoints([list(az*k/(2*np.pi)) for k in ks],with_eigenvectors=True)
    qd=ph.get_qpoints_dict(); fr=qd['frequencies']; ev=qd['eigenvectors']
    order0=np.argsort(fr[0])[:3]; pol=[]
    for b in order0:
        e=ev[0][:,b].reshape(-1,3)
        pol.append('LA' if np.sum(np.abs(e[:,2])**2)>np.sum(np.abs(e[:,:2])**2) else 'TA')
    if pol.count('LA')!=1:
        zf=[np.sum(np.abs(ev[0][:,b].reshape(-1,3)[:,2])**2)/
            max(np.sum(np.abs(ev[0][:,b].reshape(-1,3))**2),1e-30) for b in order0]
        pol=['TA']*3; pol[int(np.argmax(zf))]='LA'
    ac=np.sort(fr,axis=1)[:,:3]; v=[]
    for b in range(3):
        v.append(2*np.pi*abs(np.polyfit(ks,ac[:,b]*1e12,1)[0]))
    ta=[v[i] for i in range(3) if pol[i]=='TA']; la=[v[i] for i in range(3) if pol[i]=='LA'][0]
    return dict(vTA1=float(min(ta)),vTA2=float(max(ta)),vLA=float(la))

def amm_transmission(v1,v2,Z1,Z2,n=4001):
    thc=np.arcsin(min(1.0,v1/v2)) if v2>v1 else np.pi/2
    th=np.linspace(0,thc,n)
    s2=np.clip((v2/v1)*np.sin(th),-1,1)
    c1,c2=np.cos(th),np.cos(np.arcsin(s2))
    al=4*Z1*Z2*c1*c2/(Z1*c1+Z2*c2)**2
    return float(np.trapezoid(al*np.cos(th)*np.sin(th),th)/
                 np.trapezoid(np.cos(th)*np.sin(th),np.linspace(0,np.pi/2,n)))

def dndT(f,T):
    x=np.clip(H*f*1e12/(KB*T),1e-12,500); e=np.exp(x)
    return (H*f*1e12/(KB*T*T))*e/(e-1)**2

g0=json.load(open('s2_epitaxial.json'))
d1B,d2B=sorted(g0['Bi2Te3']['spacings'])[:2]
d1S,d2S=sorted(g0['Sb2Te3']['spacings'])[:2]

# ---- geometry FROZEN at the ti=0 (200 K) QHA solution ----
ti=0
Ft=[0.5*(Fmin_c('Bi',a,ti)[0]+Fmin_c('Sb',a,ti)[0]) for a in agrid]
p=np.polyfit(agrid,Ft,2); aT=-p[1]/(2*p[0])
cB=np.polyval(np.polyfit(agrid,[Fmin_c('Bi',a,ti)[1] for a in agrid],2),aT)
cS=np.polyval(np.polyfit(agrid,[Fmin_c('Sb',a,ti)[1] for a in agrid],2),aT)
print(f"frozen geometry (200 K QHA): a={aT:.4f}  cQL_Bi={cB:.4f}  cQL_Sb={cS:.4f}")

fl,sp,rho={},{},{}
for tag,cQL,d1,d2,cat in (('Bi',cB,d1B,d2B,'Bi'),('Sb',cS,d1S,d2S,'Sb')):
    at=r3m(aT,cQL,d1,d2,cat); at.calc=CPUNEP(POT)
    BFGS(at,logfile=None).run(fmax=1e-4,steps=300)
    cc=float(np.array(at.get_cell())[2,2])
    ph=phonons(at); fl[tag]=flux(ph); sp[tag]=sound_speeds_z(ph,cc)
    V=float(abs(np.linalg.det(np.array(at.get_cell()))))*1e-30
    rho[tag]=sum(MASS[s] for s in at.get_chemical_symbols())*AMU/V

P1,P2=fl['Bi'],fl['Sb']; den=P1+P2
aD=np.where(den>0,P2/np.where(den>0,den,1),0.0)
hw=H*fgrid*1e12
amm={}
for pol in ('vTA1','vTA2','vLA'):
    v1,v2=sp['Bi'][pol],sp['Sb'][pol]
    amm[pol]=amm_transmission(v1,v2,rho['Bi']*v1,rho['Sb']*v2)
aA=float(np.mean(list(amm.values())))
print(f"alpha_AMM (geometry-independent here) = {aA:.4f}   per-pol "
      + ", ".join(f"{k} {v:.4f}" for k,v in amm.items()))

out=[]
print(f"\n{'T':>4} {'G_DMM':>8} {'DMM_mod':>8} {'G_AMM':>8} {'G_rad,Bi':>9} {'G_rad,Sb':>9} {'a_DMM':>7}")
print("-"*62)
for T in TS:
    w=hw*P1*dndT(fgrid,T)
    G_DMM=float(np.sum(w*aD)*dfS)/1e6
    G_radB=float(np.sum(hw*P1*dndT(fgrid,T))*dfS)/1e6
    G_radS=float(np.sum(hw*P2*dndT(fgrid,T))*dfS)/1e6
    G_AMM=aA*G_radB
    print(f"{T:4d} {G_DMM:8.2f} {2*G_DMM:8.2f} {G_AMM:8.2f} {G_radB:9.2f} {G_radS:9.2f} {G_DMM/G_radB:7.3f}")
    out.append(dict(T=T,G_DMM=G_DMM,G_DMM_modified_Landauer=2*G_DMM,G_AMM=G_AMM,
                    G_rad_Bi=G_radB,G_rad_Sb=G_radS,alpha_DMM_eff=G_DMM/G_radB,alpha_AMM=aA))
json.dump(dict(frozen_at_T=200,a=float(aT),cQL_Bi=float(cB),cQL_Sb=float(cS),
               sound=sp,rho=rho,rows=out),
          open('s25_fixed_geometry_FINAL.json','w'),indent=2)

q=json.load(open('s24_models_FINAL.json'))
print("\n--- QHA-expanded (s19) vs frozen geometry (this) ---")
print(f"{'T':>4} {'DMM_qha':>8} {'DMM_fix':>8} | {'AMM_qha':>8} {'AMM_fix':>8} | {'rad_qha':>8} {'rad_fix':>8}")
for a_,b_ in zip(q,out):
    print(f"{a_['T']:4d} {a_['G_DMM']:8.2f} {b_['G_DMM']:8.2f} | "
          f"{a_['G_AMM']:8.2f} {b_['G_AMM']:8.2f} | {a_['G_rad_Bi']:8.2f} {b_['G_rad_Bi']:8.2f}")
def pc(v): return 100*(v[-1]/v[0]-1)
print(f"\n200->500 K change:")
print(f"  DMM  QHA {pc([r['G_DMM'] for r in q]):+6.2f}%   frozen {pc([r['G_DMM'] for r in out]):+6.2f}%")
print(f"  AMM  QHA {pc([r['G_AMM'] for r in q]):+6.2f}%   frozen {pc([r['G_AMM'] for r in out]):+6.2f}%")
print(f"  rad  QHA {pc([r['G_rad_Bi'] for r in q]):+6.2f}%   frozen {pc([r['G_rad_Bi'] for r in out]):+6.2f}%")
print("\nwrote s19c_fixed_geometry.json")
