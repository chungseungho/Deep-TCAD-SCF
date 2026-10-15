"""Convergence comparisons on real SCF potential snapshots only."""
import json
from dataclasses import replace
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from quantum import make_bases,solve_quantum,get_onsite


def verify_snapshot(cfg,axes,phi,vd,out):
    """Transport discretization checks at a *fixed* real SCF potential.

    Changes in self-consistent charge are NOT included; these plots validate
    single-shot transport consistency and are not whole-SCF convergence proofs.
    """
    out.mkdir(parents=True,exist_ok=True)
    energy_steps=[.01,.005,.0025]
    b=make_bases(cfg)
    energy_curr=[]
    for de in energy_steps:
        I,*_=solve_quantum(cfg,phi,b,vd,de,density=False)
        energy_curr.append(abs(I))
    erel=abs(energy_curr[-2]-energy_curr[-1])/max(energy_curr[-1],1e-20)
    fig,ax=plt.subplots(figsize=(6.5,4))
    ax.plot(energy_steps,energy_curr,'o-',lw=2)
    ax.set(xlabel='Energy integration step (eV)',ylabel='|Id| (A)',title='Energy-grid consistency on converged potential')
    ax.grid(alpha=.25);fig.tight_layout();fig.savefig(out/'06_energy_convergence.png',dpi=175);plt.close(fig)

    mode_list=sorted(set([max(2,cfg.modes-2),cfg.modes,cfg.modes+2]))
    mode_curr=[]
    for nm in mode_list:
        c=replace(cfg,modes=nm)
        cur,*_=solve_quantum(c,phi,make_bases(c),vd,.005,density=False)
        mode_curr.append(abs(cur))
    mrel=abs(mode_curr[-2]-mode_curr[-1])/max(mode_curr[-1],1e-20)
    fig,ax=plt.subplots(figsize=(6.5,4))
    ax.plot(mode_list,mode_curr,'o-',lw=2)
    ax.set(xlabel='Retained transverse modes per valley',ylabel='|Id| (A)',
           title='Mode truncation convergence (fixed potential)')
    ax.grid(alpha=.25);fig.tight_layout();fig.savefig(out/'07_mode_convergence.png',dpi=175);plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,4.4))
    for bs in b:
        H,es,ed=get_onsite(phi,bs,cfg.ec0_ev)
        # Plot true local eigenenergies in chosen fixed-mode subspace.
        levels=np.linalg.eigvalsh(H).real-2*bs.hopping_ev
        x=axes[0]
        for k in range(min(3,cfg.modes)):
            ax.plot(x,levels[:,k],lw=1.3,label=f'{bs.valley.name} #{k+1}')
    for xp in (cfg.ls_nm,cfg.ls_nm+cfg.lg_nm):ax.axvline(xp,color='gray',ls=':',lw=1)
    ax.set(xlabel='x (nm)',ylabel='Local subband E (eV)',
           title='Biased transverse subbands (actual SCF potential)')
    ax.legend(ncol=3,fontsize=7);ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'08_subbands_under_bias.png',dpi=175);plt.close(fig)
    status={'potential_source':'ACTUAL DEVSIM-NEGF SCF snapshot',
        'energy_steps_eV':energy_steps,'energy_currents_A':energy_curr,
        'energy_relative_change_0p005_to_0p0025':erel,
        'energy_check_lt_3_percent':bool(erel<.03),
        'modes':mode_list,'mode_currents_A':mode_curr,
        'mode_relative_change_last_two':mrel,
        'mode_check_lt_5_percent':bool(mrel<.05),
        'important':'Fixed potential only; repeat full-SCF with mode/grid changes for quantitative verification.'}
    (out.parent/'quality_checks.json').write_text(json.dumps(status,indent=2))
    return status
