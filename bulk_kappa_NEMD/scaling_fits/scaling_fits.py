"""
Finite-size scaling 1/kappa vs 1/L for all 18 (composition, direction) series at 300 K.
Reports kappa_inf, its standard error, R^2, and the effective mean free path from the
slope, for each of the three flux/geometry routes.

    1/kappa(L) = 1/kappa_inf + Lambda/(kappa_inf * L)

so a straight line in (1/L, 1/kappa): the intercept is 1/kappa_inf and slope/intercept is
Lambda.  Ordinary least squares on the n = 4 lengths (2 dof); SE(kappa_inf) is the
intercept standard error carried through 1/b by the delta method.

Reads the runs in ../<composition>/<direction>/L_*/ through ../kappa_composition.py, so
the intercepts here and the table there are the same numbers by construction.  Writes
scaling_fits.csv, kappa_per_run.csv and the 6x3 panel figure Fig_scaling_fits.{png,pdf}.
The per-series figures Fig_scaling_fit_<comp>_<dir>.{png,pdf} shipped alongside are the
same 18 fits drawn one at a time; see README.md.
"""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.dirname(OUT))
from kappa_composition import COMPS,DIRS,LENS,XSB,analyse

LABEL={'Bi2Te3':r'Bi$_2$Te$_3$','BiSbTe20':r'20% Sb','BiSbTe40':r'40% Sb',
       'BiSbTe60':r'60% Sb','BiSbTe80':r'80% Sb','Sb2Te3':r'Sb$_2$Te$_3$'}
DIRNAME={'X':r'$x$ (in-plane)','Y':r'$y$ (in-plane)','Z':r'$z$ (cross-plane)'}

def fit_with_se(L_nm,k):
    """Least squares 1/k = a*(1/L) + b.  Returns dict with kappa_inf, its SE, R^2, MFP."""
    x=1.0/np.asarray(L_nm,float); y=1.0/np.asarray(k,float)
    n=len(x)
    A=np.vstack([x,np.ones(n)]).T
    coef,_,_,_=np.linalg.lstsq(A,y,rcond=None)
    a,b=coef
    resid=y-A@coef
    dof=n-2
    s2=float(resid@resid)/dof                       # residual variance
    cov=s2*np.linalg.inv(A.T@A)
    se_a,se_b=np.sqrt(np.diag(cov))
    ss_tot=float(((y-y.mean())**2).sum())
    r2=1.0-float(resid@resid)/ss_tot if ss_tot>0 else np.nan
    kinf=1.0/b
    se_kinf=se_b/b**2                               # delta method on 1/b
    mfp=a/b                                         # kappa_inf * slope, in nm
    se_mfp=abs(mfp)*np.sqrt((se_a/a)**2+(se_b/b)**2) if a!=0 else np.nan
    return dict(a=a,b=b,se_a=se_a,se_b=se_b,kinf=kinf,se_kinf=se_kinf,
                r2=r2,mfp=mfp,se_mfp=se_mfp,dof=dof,x=x,y=y)

# ---- gather the runs from the merged tree -------------------------------
DB={}
for c in COMPS:
    for d in DIRS:
        for L in LENS:
            rel=os.path.join(c,d,L)
            r=analyse(os.path.join(os.path.dirname(OUT),rel),d)
            if r is not None:
                DB[(c,d,L)]=dict(r,path=rel,ratio_jp_th=r['k_jp']/r['k_th'])

# ---- fit ---------------------------------------------------------------
rows=[]
FITS={}
for c in COMPS:
    for d in DIRS:
        Ls=[];kt=[];kj=[];kp=[]
        for L in LENS:
            r=DB.get((c,d,L))
            if r is None: continue
            Ls.append(r['L_nm']); kt.append(r['k_th']); kj.append(r['k_jp']); kp.append(r['k_pub'])
        ft=fit_with_se(Ls,kt); fj=fit_with_se(Ls,kj); fp=fit_with_se(Ls,kp)
        FITS[(c,d)]=(ft,Ls,kt)
        rows.append(dict(composition=c,xSb_pct=XSB[c],direction=d,n=len(Ls),
            L_min_nm=round(min(Ls),2),L_max_nm=round(max(Ls),2),
            kappa_inf_thermostat=round(ft['kinf'],4), se_kappa_inf=round(ft['se_kinf'],4),
            rel_se_pct=round(100*ft['se_kinf']/ft['kinf'],2), R2_thermostat=round(ft['r2'],4),
            MFP_nm=round(ft['mfp'],2), se_MFP_nm=round(ft['se_mfp'],2),
            kappa_inf_jp=round(fj['kinf'],4), R2_jp=round(fj['r2'],4),
            kappa_inf_published_recipe=round(fp['kinf'],4), R2_published=round(fp['r2'],4)))

with open(os.path.join(OUT,'scaling_fits.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

# ---- per-run CSV regenerated from the merged tree ----------------------
with open(os.path.join(OUT,'kappa_per_run.csv'),'w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['composition','xSb_pct','direction','length_tag','rows','L_nm',
                'A_used_A2','A_true_A2','kappa_published_recipe','kappa_jp_corrected',
                'kappa_thermostat','ratio_jp_over_thermostat','src_snk_imbalance','path'])
    for c in COMPS:
        for d in DIRS:
            for L in LENS:
                r=DB.get((c,d,L))
                if r is None: continue
                w.writerow([c,XSB[c],d,L,r['nrows'],round(r['L_nm'],3),
                            round(r['A_used'],2),round(r['A_true'],2),round(r['k_pub'],4),
                            round(r['k_jp'],4),round(r['k_th'],4),round(r['ratio_jp_th'],4),
                            round(r['imb'],4),r['path']])

# ---- the 18-panel figure ----------------------------------------------
plt.rcParams.update({'font.size':8,'axes.linewidth':0.8,'mathtext.default':'regular'})
fig,axes=plt.subplots(6,3,figsize=(8.0,12.6),sharex=False)
# Block-averaged total uncertainty (BLOCK_AVERAGING, 2026-09-26). The regression SE
# below is the fit's own statement about its intercept; with four lengths it carries
# 2 dof and for the alloy cross-plane series it is 3-7x SMALLER than the measured
# run-to-run noise. Both are drawn: the thin capped bar is the regression SE, the thick
# pale bar behind it is SE_tot = max(SE_stat, SE_regr), which is what the tables quote.
import csv as _csv
_SE_TOT={}
_bp=os.path.join(OUT,'kappa_inf_uncertainty.csv')
if os.path.exists(_bp):
    for r in _csv.DictReader(open(_bp)):
        _SE_TOT[(r['comp'],r['dir'])]=float(r['se_total'])

for i,c in enumerate(COMPS):
    for j,d in enumerate(DIRS):
        ax=axes[i,j]
        ft,Ls,kt=FITS[(c,d)]
        x,y=ft['x'],ft['y']
        xx=np.linspace(0,x.max()*1.08,100)
        yy=ft['a']*xx+ft['b']
        weak = ft['r2']<0.90
        col='#c0392b' if weak else '#1f4e79'
        ax.plot(xx,yy,'-',lw=1.1,color=col,zorder=1)
        ax.plot(x,y,'o',ms=4.5,mfc='white',mec=col,mew=1.2,zorder=3)
        # SE_tot lives in kappa space; convert to the 1/kappa axis: d(1/k) = se/k^2
        se_b_tot=_SE_TOT.get((c,d),ft['se_kinf'])/ft['kinf']**2
        ax.errorbar(0,ft['b'],yerr=se_b_tot,fmt='none',ecolor=col,alpha=0.30,
                    elinewidth=4.0,capsize=0,zorder=2)
        ax.errorbar(0,ft['b'],yerr=ft['se_b'],fmt='s',ms=4.5,color=col,
                    capsize=2.5,lw=1.0,zorder=4)
        ax.set_xlim(left=0)
        wide=max(ft['se_b'],se_b_tot)
        lo=min(ft['b']-2*wide, y.min()); hi=max(y.max(), (ft['a']*x.max()*1.08+ft['b']))
        pad=0.12*(hi-lo)
        ax.set_ylim(lo-pad, hi+pad)
        txt=(r'$\kappa_\infty$ = %.3f $\pm$ %.3f'%(ft['kinf'],
             _SE_TOT.get((c,d),ft['se_kinf']))+'\n'
             r'$R^2$ = %.3f'%ft['r2']+'\n'
             r'$\Lambda$ = %.0f nm'%ft['mfp'])
        ax.text(0.96,0.06,txt,transform=ax.transAxes,ha='right',va='bottom',
                fontsize=6.8,color=col,
                bbox=dict(fc='white',ec=col,lw=0.5,alpha=0.85,pad=2))
        if weak:
            ax.text(0.04,0.94,'weak fit',transform=ax.transAxes,ha='left',va='top',
                    fontsize=6.5,color='#c0392b',style='italic')
        if i==0: ax.set_title(DIRNAME[d],fontsize=9)
        if j==0: ax.set_ylabel(LABEL[c]+'\n'+r'1/$\kappa$  (mK/W)',fontsize=8)
        if i==5: ax.set_xlabel(r'1/$L$  (nm$^{-1}$)',fontsize=8)
        ax.tick_params(labelsize=7)
fig.suptitle('Finite-size scaling of the NEMD thermal conductivity at 300 K\n'
             r'$1/\kappa(L)=1/\kappa_\infty+\Lambda/(\kappa_\infty L)$'
             '   (thermostat flux, true hexagonal cross-section)',fontsize=9.5,y=0.997)
fig.tight_layout(rect=[0,0,1,0.975])
fig.savefig(os.path.join(OUT,'Fig_scaling_fits.png'),dpi=300)
fig.savefig(os.path.join(OUT,'Fig_scaling_fits.pdf'))
print('wrote scaling_fits.csv, kappa_per_run.csv, Fig_scaling_fits.{png,pdf}')

# ---- console table -----------------------------------------------------
print()
print(f"{'comp':<9}{'dir':<4}{'kappa_inf':>10}{'SE':>8}{'rel%':>7}{'R2':>8}{'MFP_nm':>9}")
for r in rows:
    flag=' <-- weak' if r['R2_thermostat']<0.90 else ''
    print(f"{r['composition']:<9}{r['direction']:<4}{r['kappa_inf_thermostat']:>10.3f}"
          f"{r['se_kappa_inf']:>8.3f}{r['rel_se_pct']:>7.1f}{r['R2_thermostat']:>8.3f}"
          f"{r['MFP_nm']:>9.1f}{flag}")
