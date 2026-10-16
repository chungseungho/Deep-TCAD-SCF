"""Plan and execute fresh GAA points for joint mesh/mode/contact studies.

Failure is evidence, not a missing value to interpolate. No plateau is inferred
from two levels or from zero equilibrium current. Default plan is not executed.
"""
import sys,json,subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
from config import Config


def cases():
    base=Config()
    axes={'dx':[.5,.25,.125],'modes':[(2,4),(4,8),(6,12)],
          'transverse':[1.,.5,.25],'reservoir':[4.,6.,8.],
          'contact_length':[12.,24.,36.],
          'poisson_mesh':[(3.5,1.8,1.3),(2.5,1.2,.9),(1.8,.8,.6)]}
    output=[]
    for axis,levels in axes.items():
        for i,value in enumerate(levels):
            settings={}
            if axis=='dx':settings={'dx':value}
            if axis=='modes':settings={'modes':value[0],'sd_modes':value[1]}
            if axis=='transverse':settings={'dy':value,'dz':value}
            if axis=='reservoir':settings={'reservoir_left_nm':value,'reservoir_right_nm':value}
            if axis=='contact_length':settings={'ls':value,'ld':value,'reservoir_left_nm':value/2,'reservoir_right_nm':value/2}
            if axis=='poisson_mesh':settings=dict(zip(['mesh_h','channel_h','junction_h'],value))
            output.append({'name':f'{axis}_{i}','config':replace(base,**settings).to_dict(),'vg':.8,'vd':.1})
    for basis in ('fixed','local'):
        for junction in ('plain','nested'):
            output.append({'name':f'basis_{basis}_{junction}','config':replace(base,basis=basis,junction_basis=junction).to_dict(),'vg':.8,'vd':.1})
    return output


if __name__=='__main__':
    out=Path(sys.argv[1]) if len(sys.argv)>1 else Path('results/numerical_plan');out.mkdir(parents=True,exist_ok=True)
    execute='--execute' in sys.argv
    for case in cases():
        f=out/(case['name']+'.json');f.write_text(json.dumps(case['config'],indent=2))
        if execute:subprocess.run([sys.executable,'run.py','point','--config',str(f),'--out',str(out/case['name']),
            '--vg',str(case['vg']),'--vd',str(case['vd'])],check=False)
    (out/'execution_status.json').write_text(json.dumps({'executed':execute,'cases':len(cases()),
        'plateau_certified':False,'acceptance':'at least three successive levels; density, barrier and nonzero current; all solves valid'},indent=2))
