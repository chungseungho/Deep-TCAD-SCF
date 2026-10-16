#!/usr/bin/env python3
"""Homogeneous chain finite-bath -> ideal semi-infinite Landauer limit.

This is a quantum/contact benchmark, not a manufactured transistor I-V curve.
"""
from pathlib import Path
import argparse,csv,json
import numpy as np
from scipy.integrate import quad
from transport.finite_spectral import solve_finite_spectral
from multispace.occupation import fermi
from config import Q,HBAR_EV_S


def run(directory):
    out=Path(directory);out.mkdir(parents=True,exist_ok=True)
    t=.3;ml=.3;mr=.2;temp=300.;spin=2
    # For a perfect chain between matching semi-infinite leads, T=1 in [0,4t].
    ideal=spin*Q/HBAR_EV_S/(2*np.pi)*quad(lambda E:fermi(E,ml,temp)-fermi(E,mr,temp),0,4*t,epsabs=1e-12)[0]
    rows=[]
    for length in [8,24,64]:
        N=2*length+9;H=np.diag(np.full(N,2*t))+np.diag(np.full(N-1,-t),1)+np.diag(np.full(N-1,-t),-1)
        regions=np.r_[np.zeros(length,int),np.ones(9,int),np.full(length,2)]
        for gamma in [.005,.01,.03,.06]:
            C,d,_=solve_finite_spectral(H,regions,gamma,ml,mr,temp)
            currents=-spin*2*Q/HBAR_EV_S*np.imag(-t*np.diag(C,-1))
            current=float(np.mean(currents[length-1:length+9]));error=abs(current/ideal-1)
            rows.append(dict(lead_sites=length,gamma_ev=gamma,current_a=current,ideal_landauer_a=ideal,
                             relative_error=error,continuity_relative=float(np.ptp(currents[length-1:length+9])/max(abs(current),1e-30))))
    with (out/'finite_vs_semiinfinite.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4))
    for length in [8,24,64]:
        r=[v for v in rows if v['lead_sites']==length]
        ax.semilogx([v['gamma_ev'] for v in r],[v['current_a']/ideal for v in r],'o-',label=f'{length} sites/electrode')
    ax.axhline(1,c='black',ls='--',label='ideal semi-infinite Landauer')
    ax.set(xlabel='External bath width (eV)',ylabel='I / ideal current',title='Homogeneous chain contact benchmark (not GAA)')
    ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout();fig.savefig(out/'finite_vs_semiinfinite.png',dpi=160);plt.close(fig)
    summary={'ideal_current_a':ideal,'best_relative_error':min(v['relative_error'] for v in rows),'results':rows}
    (out/'quantum_benchmark.json').write_text(json.dumps(summary,indent=2));return summary


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',default='output/quantum_benchmark');a=p.parse_args()
    print(json.dumps(run(a.out),indent=2))
