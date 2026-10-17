"""Plots are generated only from recorded calculations; failures stay visible."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def state_plot(path,out):
    s=np.load(path);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    x=s['x_nm'];xyz=s['xyz_nm'];n=s['electron_density_cm3'];phi=s['phi_v']
    means=lambda values:np.array([values[xyz[:,0]==xx].mean() for xx in x])
    fig,ax=plt.subplots(2,1,figsize=(8,6),sharex=True,layout='constrained')
    ax[0].plot(x,-means(phi));ax[0].set_ylabel('Mean Ec (eV)')
    ax[1].semilogy(x,np.maximum(means(n),1));ax[1].set_ylabel('Mean electrons (cm$^{-3}$)');ax[1].set_xlabel('Longitudinal coordinate (nm)')
    fig.suptitle(f'{str(s["backend"])} | Vg={float(s["vg_v"]):g} V, Vd={float(s["vd_v"]):g} V | provisional')
    fig.savefig(out/'profiles.png',dpi=160);plt.close(fig)
    fig=plt.figure(figsize=(8,5),layout='constrained');ax=fig.add_subplot(projection='3d')
    select=np.arange(0,len(xyz),max(1,len(xyz)//10000));v=xyz[select]
    dots=ax.scatter(*v.T,c=np.log10(np.maximum(n[select],1)),s=2,cmap='viridis')
    ax.set(xlabel='Longitudinal (nm)',ylabel='Transverse 1 (nm)',zlabel='Transverse 2 (nm)')
    fig.colorbar(dots,ax=ax,label='log10 n (cm$^{-3}$)',shrink=.6);fig.savefig(out/'density3d.png',dpi=160);plt.close(fig)


def convergence_plot(history,out):
    rows=list(csv.DictReader(Path(history).open()));fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    for key in dict.fromkeys(r['key'] for r in rows):
        part=[r for r in rows if r['key']==key]
        ax.semilogy([int(r['iteration']) for r in part],[float(r['phi_residual_v']) for r in part],label=key)
    ax.set(xlabel='SCF iteration',ylabel='Unmixed potential residual (V)');ax.legend(fontsize=7)
    fig.savefig(Path(out)/'scf.png',dpi=160);plt.close(fig)


def compare_plot(file,out):
    s=np.load(file);xyz=s['xyz_nm'];xs=np.unique(xyz[:,0]);fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    for name,label in [('n_ms_cm3','Closed MS density'),('n_negf_cm3','Open NEGF at same field')]:
        ax.semilogy(xs,[max(s[name][xyz[:,0]==x].mean(),1) for x in xs],label=label)
    ax.set(xlabel='Longitudinal coordinate (nm)',ylabel='Mean n (cm$^{-3}$)',title='Frozen-field density comparison; provisional')
    ax.legend();fig.savefig(Path(out)/'ms_negf_density.png',dpi=160);plt.close(fig)
