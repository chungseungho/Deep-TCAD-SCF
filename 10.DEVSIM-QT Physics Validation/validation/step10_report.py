"""Export real computed states and conservative validation flags."""
import json,csv
from pathlib import Path
import numpy as np


def main():
    root=Path(__file__).resolve().parents[1];out=root/'results/assessment';out.mkdir(parents=True,exist_ok=True)
    rows=[]
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(8,7),sharex=True)
    for folder in sorted((root/'results').glob('final_*')):
        for path in sorted((folder/'states').glob('*.json')):
            data=json.loads(path.read_text());state=np.load(path.with_suffix('.npz'))
            c=json.loads(str(state['config_json']));name=folder.name
            row={'run':name,'geometry':c.get('geometry_profile','raised'),'vg_v':data['vg_v'],'vd_v':data['vd_v'],
                 'id_a':data['id_a'],'phi_residual_v':data['phi_residual_v'],
                 'density_relative_residual':data['density_relative_residual'],
                 'gauss_relative_error':data['gauss_relative_error'],'current_conservation_relative':data['current_conservation_relative'],
                 'SCF_CONVERGED':data['numerical_status']=='SCF_CONVERGED',
                 'TRANSPORT_CONSERVED':data['current_conservation_abs_a']<max(1e-15,abs(data['id_a'])*1e-6),
                 'NUMERICAL_CONVERGENCE_VERIFIED':False,'PHYSICAL_MODEL_VALID':False,
                 'QUANTITATIVE_VALIDATED':False,'TRAINING_READY':False}
            rows.append(row)
            qx=state['xyz_nm'][:,0];x=np.unique(qx)
            phi=np.array([state['phi_v'][qx==v].mean() for v in x])
            density=np.array([state['electron_density_cm3'][qx==v].mean() for v in x])
            axes[0].plot(x,phi,label=name);axes[1].semilogy(x,np.maximum(density,1.),label=name)
    if rows:
        with (out/'computed_bias_points.csv').open('w',newline='') as file:
            writer=csv.DictWriter(file,list(rows[0]));writer.writeheader();writer.writerows(rows)
        axes[0].set(ylabel='Cross-section mean potential (V)')
        axes[1].set(xlabel='Longitudinal position (nm)',ylabel='Mean density (cm$^{-3}$)')
        axes[0].legend(fontsize=6);fig.tight_layout();fig.savefig(out/'actual_scf_profiles.png',dpi=180)
    plt.close(fig)
    benchmark=root/'results/independent_final/summary.json'
    passed=json.loads(benchmark.read_text())['pass'] if benchmark.exists() else False
    attempts=[]
    for path in sorted((root/'results').glob('final_*/run/execution.json')):
        record=json.loads(path.read_text());attempts.append({'run':path.parents[1].name,**record})
    status={'STEP':10,'INDEPENDENT_QUANTUM_BENCHMARKS_PASS':passed,
        'DEVICE_SIMULATION_VERIFIED':False,'QUANTITATIVE_VALIDATED':False,'TRAINING_READY':False,
        'numerically_converged_real_scf_states':sum(r['SCF_CONVERGED'] for r in rows),
        'reason':'Device mode/contact/grid/energy convergence has not jointly been established; unresolved isolated-state preparation and contact-cut validity remain explicit limitations.',
        'states':rows,'actual_run_records':attempts,
        'units':{'coordinates':'nm','density':'cm^-3','potential':'V','energy':'eV','current':'A'},
        'iv_scope':'Only actual SCF output points exported. No missing bias point is interpolated or fabricated. Requested transfer/output sweeps require validated biased SCF.'}
    (root/'STATUS.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps({'actual_states':len(rows),'SCF_converged':status['numerically_converged_real_scf_states'],'DEVICE_SIMULATION_VERIFIED':False}),flush=True)

if __name__=='__main__':main()
