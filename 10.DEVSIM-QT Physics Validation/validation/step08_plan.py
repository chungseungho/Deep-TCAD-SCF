"""Generate a conditional HPC plan. Does not submit or execute any calculations."""
import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
from config import Config


def run(out='hpc_plan'):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);rows=[]
    base=Config(backend='negf',newton_start=1,max_iterations=80)
    # These are exploration settings, not a certified production configuration.
    for m,sd in [(36,72),(48,96),(60,120)]:
        for L in (12.,24.,36.):
            c=replace(base,modes=m,sd_modes=sd,ls=L,ld=L,reservoir_left_nm=L/2,reservoir_right_nm=L/2)
            name=f'negf_m{m}_sd{sd}_L{L:g}';path=out/(name+'.json');path.write_text(json.dumps(c.to_dict(),indent=2))
            rows.append({'name':name,'config':str(path),'phase':'frozen/equilibrium first','executed':False})
    for label,settings in [('dx025',{'dx':.25}),('dx0125',{'dx':.125}),('yz05',{'dy':.5,'dz':.5}),
        ('yz025',{'dy':.25,'dz':.25}),('poisson_mid',{'mesh_h':2.5,'channel_h':1.2,'junction_h':.9}),
        ('poisson_fine',{'mesh_h':1.8,'channel_h':.8,'junction_h':.6}),
        ('scf_tight',{'phi_tol':2e-5,'charge_tol':4e-4,'current_tol':4e-4}),
        ('energy_tight',{'energy_step':.005,'energy_order':12,'energy_rtol':4e-6,'energy_atol':2e-8,'energy_max_nodes':480000})]:
        c=replace(base,modes=60,sd_modes=120,**settings).validate();name='conditional_'+label;path=out/(name+'.json')
        path.write_text(json.dumps(c.to_dict(),indent=2));rows.append({'name':name,'config':str(path),'phase':'ONLY after mode/contact gate','executed':False})
    (out/'plan.json').write_text(json.dumps({'executed':False,'jobs_submitted':0,'cases':rows,
        'gate':'3 successive refinements stable in regional density, charge, barrier and independently integrated nonzero current; valid occupations required',
        'configuration_certified':False},indent=2))
    print('Generated',len(rows),'unexecuted configurations')

if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else 'hpc_plan')
