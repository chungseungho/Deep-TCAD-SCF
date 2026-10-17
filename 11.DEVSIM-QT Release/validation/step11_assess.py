"""Generate measured tables; no synthetic currents or interpolated successes."""
import json,csv,sys,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from config import Config
from validation.release_cli import compare

def run():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scipy.interpolate import RegularGridInterpolator
    out=ROOT/'results/assessment';out.mkdir(exist_ok=True);states=[]
    for f in sorted((ROOT/'results').glob('*/states/*.npz')):
        if f.parent.parent.name in ('raised_eq',):continue # superseded pre-memory-fix run
        with np.load(f,allow_pickle=False) as s:
            c=Config(**json.loads(str(s['config_json'])));states.append({'file':str(f.relative_to(ROOT)),'geometry':c.geometry_profile,'modes':c.modes,'sd_modes':c.sd_modes,'vg_v':float(s['vg_v']),'vd_v':float(s['vd_v']),'id_a':float(s['id_a']),'SCF_CONVERGED':bool(s['converged']),'runtime_s':float(s['runtime_s']),'teacher_ready':False})
    if states:
        with (out/'computed_bias_points.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(states[0]));w.writeheader();w.writerows(states)
    for geometry in ('raised','uniform'):
        # Keep all real measurements in CSV, draw one series per configuration.
        rows=[s for s in states if s['geometry']==geometry and s['SCF_CONVERGED']]
        with (out/f'{geometry}_actual_iv.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(states[0]) if states else ['vg_v','vd_v','id_a']);w.writeheader();w.writerows(rows)
        fig,ax=plt.subplots(figsize=(6,4))
        for modes,sd in sorted({(r['modes'],r['sd_modes']) for r in rows}):
            points=sorted({(r['vd_v'],r['id_a']) for r in rows if (r['modes'],r['sd_modes'])==(modes,sd) and r['vg_v']==.4})
            if points:ax.plot([p[0] for p in points],[p[1]*1e9 for p in points],'o',label=f'{modes}/{sd} requested modes')
        ax.set(xlabel='Drain bias (V)',ylabel='Drain current (nA)',title=f'{geometry} GAA: actual SCF points, not certified I–V');ax.grid(alpha=.25)
        if rows:ax.legend()
        fig.tight_layout();fig.savefig(out/f'{geometry}_actual_iv.png',dpi=150);plt.close(fig)
    frozen=ROOT/'results/frozen_modes/summary.json'
    if frozen.exists():
        rows=json.loads(frozen.read_text())['rows'];flat=[]
        for r in rows:
            flat.append({k:v for k,v in r.items() if not isinstance(v,(list,dict))})
        keys=list(dict.fromkeys(k for row in flat for k in row))
        with (out/'mode_projection_convergence.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(flat)
    a=ROOT/'results/raised_eq_final/states/vg0.400000_vd0.000000.npz';b=ROOT/'results/raised_refined_eq/states/vg0.400000_vd0.000000.npz'
    if a.exists() and b.exists():(out/'fully_scf_mode_comparison.json').write_text(json.dumps(compare([a,b]),indent=2))
    b=ROOT/'results/raised_resume/states/vg0.400000_vd0.000000.npz'
    if a.exists() and b.exists():
        report=compare([a,b]);row=report['comparisons'][0];report['scope']='equilibrium checkpoint reproducibility, not transport-mode convergence'
        report['status']='PASS' if row['barrier_phi_max_difference_v']<1e-6 and max(row[k] for k in ('channel_density_l2','junction_extension_density_l2','reservoir_density_l2'))<1e-4 else 'FAIL'
        (out/'restart_comparison.json').write_text(json.dumps(report,indent=2))
    matching=[]
    for row in states:
        if not row['SCF_CONVERGED']:continue
        f=ROOT/row['file'];s=np.load(f);lead=f.parent.parent/'lead_reference.npz'
        if not lead.exists():continue
        r=np.load(lead);c=Config(**json.loads(str(s['config_json'])));xyz=s['xyz_nm'];x=s['x_nm'];rphi=r['phi_v'].reshape(len(r['y_nm']),len(r['z_nm']))
        fphi=RegularGridInterpolator((r['y_nm'],r['z_nm']),rphi)
        vals={}
        for name,xx,bias in [('source',x[0],0.),('drain',x[-1],float(s['vd_v']))]:
            mask=np.isclose(xyz[:,0],xx);actual=s['phi_v'][mask];target=fphi(xyz[mask,1:])+bias+c.potential_reference_v
            vals[name+'_endpoint_phi_error_max_v']=float(max(abs(actual-target)))
        matching.append({'state':row['file'],**vals})
    (out/'device_lead_matching.json').write_text(json.dumps({'scope':'cell-averaged first/last quantum slice vs translational SP reference; not constant-flatness criterion','measurements':matching},indent=2))
    have={(round(s['vg_v'],6),round(s['vd_v'],6)) for s in states if s['geometry']=='raised' and s['SCF_CONVERGED']}
    req=[{'sweep':'IdVgs','vg_v':round(i*.05,6),'vd_v':vd} for vd in (.05,.1) for i in range(25)]+[{'sweep':'IdVds','vg_v':vg,'vd_v':round(i*.05,6)} for vg in (.4,.8,1.2) for i in range(25)]
    for row in req:
        row['status']='NUMERICALLY_COMPUTED' if (row['vg_v'],row['vd_v']) in have else 'NOT_RUN';row['teacher_ready']=False
    (out/'requested_sweep_status.json').write_text(json.dumps(req,indent=2))
    # Separate source changes in validation/report scripts from production.
    inventory=[]
    excluded={'tests','validation','results','__pycache__'}
    for p in (ROOT/'results').glob('*/source_manifest.json'):
        old=json.loads(p.read_text());changed=[];missing=[]
        for rel,h in old.items():
            if any(x in excluded for x in Path(rel).parts) or not rel.endswith('.py'):continue
            current=ROOT/rel
            if not current.exists():missing.append(rel)
            elif hashlib.sha256(current.read_bytes()).hexdigest()!=h:changed.append(rel)
        inventory.append({'run':p.parent.name,'production_files_changed_since_run_start':changed,'production_files_missing':missing,'snapshot':str(p.parent.relative_to(ROOT)/'source_snapshot.zip')})
    (out/'production_provenance.json').write_text(json.dumps(inventory,indent=2))
    return states
if __name__=='__main__':print(json.dumps(run(),indent=2))
