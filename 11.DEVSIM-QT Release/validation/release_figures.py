import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]

def run():
    out=ROOT/'results/assessment';out.mkdir(exist_ok=True)
    r=np.load(ROOT/'results/lead_reference/data.npz');phi=r['phi'].reshape(20,15);density=r['density'].reshape(20,15)
    fig,axes=plt.subplots(1,2,figsize=(9,4))
    for ax,arr,label in zip(axes,[phi.T,density.T/1e19],['Lead potential (V)','Lead density / 1e19 cm^-3']):
        im=ax.imshow(arr,origin='lower',extent=[-10,10,-7.5,7.5],aspect='equal');ax.set(xlabel='y (nm)',ylabel='z (nm)',title=label);fig.colorbar(im,ax=ax)
    fig.suptitle('Transverse SP electrode: 20 x 15 nm, mean potential = 0');fig.tight_layout();fig.savefig(out/'lead_reference.png',dpi=160);plt.close(fig)
    f=ROOT/'results/raised_refined_eq/states/vg0.400000_vd0.000000.npz'
    if f.exists():
        s=np.load(f);mask=s['silicon_mask'];y=s['y_nm'];j=np.argmin(abs(y));phi=np.full(mask.shape,np.nan);n=phi.copy();phi[mask]=s['phi_v'];n[mask]=np.log10(np.maximum(s['electron_density_cm3'],1.))
        fig,axes=plt.subplots(2,1,figsize=(9,5),sharex=True)
        for ax,arr,title in zip(axes,[phi[:,j,:].T,n[:,j,:].T],['Electrostatic potential (V)','log10 electron density (cm^-3)']):
            im=ax.imshow(arr,origin='lower',extent=[0,54,-7.5,7.5],aspect='auto');ax.set(ylabel='z (nm)',title=title);fig.colorbar(im,ax=ax)
        axes[-1].set_xlabel('x (nm)');fig.suptitle('Raised S/D, VG=0.4 V, VD=0: numerical SCF, accuracy not certified');fig.tight_layout();fig.savefig(out/'raised_equilibrium_fields.png',dpi=160);plt.close(fig)
    rows=json.loads((ROOT/'results/frozen_modes/summary.json').read_text())['rows'][1:]
    fig,ax=plt.subplots(figsize=(7,4))
    labels=[f"{r['requested_modes']}/{r['requested_sd_modes']}" for r in rows]
    for region in ('channel','junction_extension','reservoir'):
        ax.plot(labels,[100*r[region+'_density_relative_l2_change'] for r in rows],'o-',label=region.replace('_',' '))
    ax.axhline(1,color='black',linestyle='--',label='1% target');ax.set(yscale='log',xlabel='Larger requested channel/contact basis',ylabel='Density change from previous basis (%)',title='Fixed potential and fixed full-reference lead μ');ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'mode_projection_convergence.png',dpi=160);plt.close(fig)
if __name__=='__main__':run()
