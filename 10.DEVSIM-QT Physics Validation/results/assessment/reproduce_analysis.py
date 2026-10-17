"""Reproduce plots/CSV from saved measurements only; no missing points filled."""
import sys,json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
OUT=Path(__file__).parent

def load(name):
    paths=sorted((ROOT/'results'/name/'states').glob('*.npz'))
    if not paths:return None
    return np.load(paths[0])

def compare(a,b):
    return {'density_relative_l2':float(np.linalg.norm(a['electron_density_cm3']-b['electron_density_cm3'])/np.linalg.norm(b['electron_density_cm3'])),
            'potential_max_abs_v':float(max(abs(a['phi_v']-b['phi_v'])))}

def main():
    from validation.step10_report import main as report
    report();summary={}
    modepath=ROOT/'results/frozen_modes/summary.json'
    if modepath.exists():
        data=json.loads(modepath.read_text());rows=[]
        for r in data['cases']:
            rows.append({k:r[k] for k in ('channel_modes','requested_sd_modes','mu_ev','electrons','density_relative_l2_to_largest','density_relative_l2_to_previous')})
        with (OUT/'frozen_mode_convergence.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
        fig,ax=plt.subplots(1,2,figsize=(9,4))
        ax[0].plot([r['requested_sd_modes'] for r in rows],[r['electrons'] for r in rows],'o-')
        ax[1].plot([r['requested_sd_modes'] for r in rows[1:]],[100*r['density_relative_l2_to_previous'] for r in rows[1:]],'o-')
        ax[0].set(xlabel='Requested S/D modes',ylabel='Electrons at the saved potential')
        ax[1].set(xlabel='Requested S/D modes',ylabel='Density change to previous (%)');ax[1].axhline(1,color='r',ls='--',label='1% criterion');ax[1].legend()
        fig.suptitle('Frozen-potential mode test; not a new self-consistent solution');fig.tight_layout();fig.savefig(OUT/'frozen_mode_convergence.png',dpi=180);plt.close(fig)
        summary['frozen_mode_plateau']=data['mode_plateau_observed_at_frozen_potential']
    base=load('final_raised_equilibrium');fine=load('final_raised_fine_poisson');restart=load('final_raised_restart');long=load('final_raised_long_contact')
    if base is not None and fine is not None:summary['poisson_mesh_refinement']=compare(base,fine)
    if base is not None and restart is not None:summary['restart_recalculation']=compare(base,restart)
    if base is not None and long is not None:
        from scipy.spatial import cKDTree
        ca=json.loads(str(base['config_json']));cb=json.loads(str(long['config_json']))
        shift=cb['ls']-ca['ls'];xyz=base['xyz_nm'].copy();xyz[:,0]+=shift
        distance,index=cKDTree(long['xyz_nm']).query(xyz)
        assert max(distance)<1e-9
        summary['contact_length_12_to_24_nm']={'density_relative_l2_in_matched_central_domain':float(np.linalg.norm(base['electron_density_cm3']-long['electron_density_cm3'][index])/np.linalg.norm(long['electron_density_cm3'][index])),
            'potential_max_abs_v_in_matched_central_domain':float(max(abs(base['phi_v']-long['phi_v'][index])))}
    points=[]
    for folder in sorted((ROOT/'results').glob('final_uniform*')):
        for p in (folder/'states').glob('*.json'):
            d=json.loads(p.read_text())
            if d['numerical_status']=='SCF_CONVERGED':points.append({'vg_v':d['vg_v'],'vd_v':d['vd_v'],'id_a':d['id_a'],'SCF_CONVERGED':True,'QUANTITATIVE_VALIDATED':False,'run':folder.name})
    if len(points)>1:
        points.sort(key=lambda r:r['vd_v'])
        with (OUT/'idvds_uniform_vg0p4_diagnostic.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,list(points[0]));w.writeheader();w.writerows(points)
        fig,ax=plt.subplots(figsize=(6,4));ax.plot([r['vd_v'] for r in points],[r['id_a']*1e6 for r in points],'o-')
        ax.set(xlabel='Drain bias (V)',ylabel='Drain current (microampere)',title='Uniform 3D GAA, Vg=0.4 V: actual SCF samples\nNumerical/physical convergence not certified')
        fig.tight_layout();fig.savefig(OUT/'idvds_uniform_vg0p4_diagnostic.png',dpi=180);plt.close(fig)
    (OUT/'numerical_comparisons.json').write_text(json.dumps(summary,indent=2)+'\n')
    (OUT/'requested_sweep_status.json').write_text(json.dumps({'original_raised_idvgs_0_to_1p2':'NOT RUN: biased contact/bound-state and mode-convergence gates unresolved','original_raised_idvds_0_to_1p2':'NOT RUN: same unresolved gates','uniform_points':'Only completed solver states are exported; no missing point filled','training_ready':False},indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
