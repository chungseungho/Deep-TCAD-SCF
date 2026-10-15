"""Only plots real saved simulations; never manufacture an IV curve."""
from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_csv(path):
    with open(path,newline='') as f: return list(csv.DictReader(f))


def make_modes_figure(cfg,bases,path):
    b=bases[0]
    count=min(4,cfg.modes)
    dy=cfg.width_nm/(cfg.ny+1); dz=cfg.thickness_nm/(cfg.nz+1)
    fig,axs=plt.subplots(1,count,figsize=(3.6*count,3.7),layout='constrained')
    axs=np.ravel(axs)
    for m,ax in enumerate(axs):
        img=(np.abs(b.U[:,m])**2/(dy*dz)).reshape(cfg.ny,cfg.nz).T
        im=ax.imshow(img,extent=(-cfg.width_nm/2,cfg.width_nm/2,-cfg.thickness_nm/2,cfg.thickness_nm/2),
                     origin='lower',aspect='equal',interpolation='bilinear',cmap='viridis')
        ax.set_title(f'{b.valley.name}, mode {m+1}\nE={b.levels[m]:.4f} eV')
        ax.set_xlabel('y (nm)')
        if m==0: ax.set_ylabel('z (nm)')
        fig.colorbar(im,ax=ax,label='|chi|^2 (nm^-2)',shrink=.7)
    fig.savefig(path,dpi=175); plt.close(fig)


def make_idvgs(out):
    out=Path(out); f=out/'idvgs.csv'
    if not f.exists(): return
    rows=read_csv(f)
    if not rows:return
    v=np.array([float(r['vg_v']) for r in rows]); I=np.array([float(r['id_a']) for r in rows])
    ok=np.array([r['status']=='CONVERGED' for r in rows])
    summary={}
    j=out/'iv_metrics.json'
    if j.exists(): summary=json.loads(j.read_text())
    fig,axs=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    axs[0].plot(v,np.abs(I),'o-',ms=3)
    axs[0].set_ylabel('|Id| (A)');axs[0].set_title('Id-Vgs (linear)')
    axs[1].semilogy(v,np.maximum(np.abs(I),1e-30),'o-',ms=3)
    axs[1].set_ylabel('|Id| (A), log');axs[1].set_title('Id-Vgs (log)')
    if len(v)>1:
        axs[2].plot(v,np.gradient(I,v),'o-',ms=3)
    axs[2].set_ylabel('gm (A/V)');axs[2].set_title('Transconductance')
    for ax in axs:
        if len(ok) and not np.all(ok): ax.scatter(v[~ok],np.maximum(np.abs(I[~ok]),1e-30) if ax is axs[1] else np.zeros(np.count_nonzero(~ok)),marker='x',color='red',label='NOT CONVERGED')
        ax.set_xlabel('Vgs (V)');ax.grid(alpha=.25)
        ax.axvspan(2,6,color='orange',alpha=.045)
    vt=summary.get('vth',{}).get('vth_v')
    if vt is not None:
        axs[0].axvline(vt,ls='--',color='black',label=f'Vth={vt:.3f} V')
        axs[1].axvline(vt,ls='--',color='black')
        axs[0].legend()
    fig.suptitle('Mono-nanosheet GAA — fully SCF mode-space NEGF (model validity requires checks)')
    fig.savefig(out/'figures'/'02_idvgs.png',dpi=175);plt.close(fig)


def make_idvds(out):
    out=Path(out); f=out/'idvds.csv'
    if not f.exists():return
    rows=read_csv(f)
    if not rows:return
    fig,axs=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    for vg in sorted(set(float(r['vg_v']) for r in rows)):
        sub=[r for r in rows if abs(float(r['vg_v'])-vg)<1e-7]
        v=np.array([float(r['vd_v']) for r in sub]);I=np.abs([float(r['id_a']) for r in sub])
        axs[0].plot(v,I,'o-',ms=3,label=f'Vgs={vg:g} V')
        axs[1].semilogy(v,np.maximum(I,1e-30),'o-',ms=3,label=f'Vgs={vg:g} V')
    for ax in axs:
        ax.set_xlabel('Vds (V)');ax.set_ylabel('|Id| (A)');ax.grid(alpha=.25)
        ax.legend(fontsize=8)
        ax.axvspan(2,6,color='orange',alpha=.045)
    axs[0].set_title('Id-Vds (linear)');axs[1].set_title('Id-Vds (log)')
    fig.savefig(out/'figures'/'03_idvds.png',dpi=175);plt.close(fig)


def make_device_maps(out,state):
    out=Path(out)
    p=np.load(state)
    phi=p['phi_v'];n=p['electron_density_cm3']
    x=p['x_nm'];y=p['y_nm'];z=p['z_nm']
    fig,ax=plt.subplots(2,2,figsize=(11,6.5),layout='constrained')
    for (a,data,label) in [(ax[0,0],phi[:,:,len(z)//2],'Potential (V), z=0'),
                            (ax[0,1],phi[:,len(y)//2,:],'Potential (V), y=0'),
                            (ax[1,0],np.log10(np.maximum(n[:,:,len(z)//2],1)),'log10 n (cm^-3), z=0'),
                            (ax[1,1],np.log10(np.maximum(n[:,len(y)//2,:],1)),'log10 n (cm^-3), y=0')]:
        yy=y if 'z=0' in label else z
        im=a.pcolormesh(x,yy,data.T,shading='auto',cmap='viridis')
        a.set(xlabel='Transport x (nm)',ylabel='Transverse (nm)',title=label)
        fig.colorbar(im,ax=a)
    fig.savefig(out/'figures'/'04_device_maps.png',dpi=175);plt.close(fig)


def make_scf_history(out):
    out=Path(out);f=out/'scf_history.csv'
    if not f.exists():return
    rows=read_csv(f)
    if not rows:return
    cases={}
    for r in rows: cases.setdefault(r['key'],[]).append(r)
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    # Only last few biases, otherwise plot unreadable.
    for key,rows in list(cases.items())[-5:]:
        ix=[int(r['iteration']) for r in rows]
        axs[0].semilogy(ix,[float(r['delta_phi_v'])+1e-15 for r in rows],'-o',ms=3,label=key)
        axs[1].semilogy(ix,[float(r['charge_residual'])+1e-15 for r in rows],'-o',ms=3,label=key)
    for a,label in zip(axs,['Max Δφ (V)','NEGF charge residual']):
        a.set(xlabel='SCF iteration',ylabel=label);a.grid(alpha=.25);a.legend(fontsize=7)
    fig.savefig(out/'figures'/'05_scf_convergence.png',dpi=175);plt.close(fig)
