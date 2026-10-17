"""Deliverable scientific plots from saved data and explicit failure statuses."""
import sys,json,csv
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analysis.plots import state_plot,convergence_plot,compare_plot


def run(root):
    root=Path(root);out=root/'figures';out.mkdir(exist_ok=True)
    plt.rcParams.update({'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
    state_plot(root/'gaa_ms_fixed_mu_equilibrium/states/vg0.800000_vd0.000000.npz',out/'ms_equilibrium')
    # Explicit margins keep 3D labels inside the exported page.
    s=np.load(root/'gaa_ms_fixed_mu_equilibrium/states/vg0.800000_vd0.000000.npz')
    xyz=s['xyz_nm'];density=s['electron_density_cm3'];select=np.arange(0,len(xyz),max(1,len(xyz)//10000))
    fig=plt.figure(figsize=(8,5));ax=fig.add_axes([.02,.12,.78,.80],projection='3d')
    points=ax.scatter(*xyz[select].T,c=np.log10(np.maximum(density[select],1)),s=2,cmap='viridis')
    ax.set(xlabel='x (nm)',ylabel='y (nm)',zlabel='z (nm)',title='GAA equilibrium MS density; canonical coordinates')
    ax.tick_params(labelsize=8);bar=fig.add_axes([.84,.23,.025,.5]);fig.colorbar(points,cax=bar,label='log10 n (cm$^{-3}$)')
    fig.savefig(out/'ms_equilibrium/density3d.png',dpi=160);plt.close(fig)
    convergence_plot(root/'gaa_ms_fixed_mu_equilibrium/scf_history.csv',out)
    compare_plot(root/'matched_mu_frozen_comparison/comparison.npz',out)
    rows=list(csv.DictReader((root/'stage_a_transfer/idvgs.csv').open()))
    fig,ax=plt.subplots(figsize=(8,2.5),layout='constrained')
    ax.scatter([float(r['vg_v']) for r in rows],np.zeros(len(rows)),marker='x',color='#b73338',s=45)
    ax.set(xlabel='Gate voltage (V)',yticks=[0],yticklabels=['Origin unresolved'],ylim=(-.5,.5),
        title='MS transfer attempts: 0/25 current-bearing SCF solutions (Vd=0.1 V)')
    fig.savefig(out/'transfer_failures.png',dpi=160);plt.close(fig)
    rows=list(csv.DictReader((root/'stage_a_output/idvds.csv').open()))
    fig,ax=plt.subplots(figsize=(9,3),layout='constrained')
    for good,color,label in [(True,'#27876a','Converged equilibrium'),(False,'#b73338','Origin unresolved')]:
        subset=[r for r in rows if (r['numerical_status']=='SCF_CONVERGED')==good]
        ax.scatter([float(r['vd_v']) for r in subset],[float(r['vg_v']) for r in subset],c=color,s=28,label=label)
    ax.set(xlabel='Drain voltage (V)',ylabel='Gate voltage (V)',yticks=[.4,.8,1.2],ylim=(.25,1.4),
        title='MS output attempts: equilibrium solved; all 72 biased points rejected')
    ax.legend(loc='upper center',ncols=2,fontsize=8);fig.savefig(out/'output_status.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    for name,field in [('stress_gate','vg_v'),('stress_drain','vd_v')]:
        rr=list(csv.DictReader((root/name/'stress.csv').open()));values=[]
        for r in rr:
            if r.get('maximum_failed_iterate_field_v_nm'):values.append((float(r[field]),float(r['maximum_failed_iterate_field_v_nm'])))
        if values:
            a=np.array(values);ax.plot(a[:,0],a[:,1],'o-',label=name.replace('_',' '))
    ax.axhline(.5,color='#b73338',ls='--',lw=1,label='Diagnostic threshold')
    ax.set(xlabel='Swept voltage (V)',ylabel='Maximum available iterate field (V/nm)',
        title='Stress diagnostics; rejected points are not SCF device predictions')
    ax.legend();fig.savefig(out/'stress_fields.png',dpi=160);plt.close(fig)
    study=json.loads((root/'frozen_subspace/summary.json').read_text())
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    rows=study['cases'];names=[r['name'].replace('_','\n') for r in rows];x=np.arange(len(rows))
    ax[0].bar(x,[r['electrons'] for r in rows],color='#477a9c');ax[0].set(xticks=x,xticklabels=names,ylabel='Electrons at fixed mu',title='Frozen GAA subspaces; no plateau claim')
    ax[1].bar(x,[r['electrons_above_mu_plus_8kbt'] for r in rows],color='#aa7a3c');ax[1].set(xticks=x,xticklabels=names,ylabel='Electrons above mu+8 kBT',title='Added-mode charge is not just this tail')
    for a in ax:a.tick_params(axis='x',labelsize=8)
    fig.savefig(out/'frozen_mode_population.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
    for p in sorted((root/'analytical').glob('barrier*.npz')):
        s=np.load(p);ax.plot(s['energy_ev'],s['transmission_rgf'],label=p.stem)
        select=np.arange(0,len(s['energy_ev']),10);ax.scatter(s['energy_ev'][select],s['transmission_dense'][select],s=10,facecolors='none',edgecolors='#222')
    ax.set(xlabel='Energy (eV)',ylabel='Transmission per mode',title='Analytical chain/barrier: RGF and independent dense inversion')
    ax.legend();fig.savefig(out/'analytical_transmission.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    data=np.load(root/'matched_mu_frozen_comparison/comparison.npz')
    ax[0].plot(data['energies_ev']-float(data['mu_l_ev']),data['weighted_transmission'],color='#477a9c')
    ax[0].set(title='GAA MS equilibrium H, semi-infinite leads',xlabel='E - muL (eV)',ylabel='Spin/valley summed T')
    data=np.load(root/'finite_current_check/transmission.npz')
    state=np.load(root/'gaa_spectral_regression/states/vg0.800000_vd0.100000.npz')
    ax[1].plot(data['energy_ev']-float(state['mu_l_ev']),data['spin_valley_summed_transmission'],color='#aa7a3c')
    ax[1].axvspan(-.1,0,color='#999',alpha=.15)
    ax[1].set(title='GAA finite-bath H, Vd=0.1 V',xlabel='E - muL (eV)',ylabel='Spin/valley summed T')
    fig.suptitle('Actual computed spectra; fields/contact models differ; no biased MS validation',fontsize=11)
    fig.savefig(out/'gaa_transmission.png',dpi=160);plt.close(fig)
    return out


if __name__=='__main__':print(run(sys.argv[1] if len(sys.argv)>1 else 'results'))
