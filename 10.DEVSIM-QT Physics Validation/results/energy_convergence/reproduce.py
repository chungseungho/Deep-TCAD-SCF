"""Frozen final uniform biased SCF field: independent energy-settings sweep."""
from pathlib import Path
import sys,json
import numpy as np
from dataclasses import replace
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from config import Config
from modespace.hamiltonian import Grid,ModeBuilder
from modespace.junction_basis import JunctionCompatibleBuilder
from device.definition import Frame
from transport.registry import create_transport_backend

def main():
 state=np.load(ROOT/'results/final_uniform_vd001/states/vg0.400000_vd0.010000.npz');base=Config(**json.loads(str(state['config_json'])))
 g=Grid.build(base);rows=[];densities=[]
 for step,order,rtol,atol,poles in ((.01,8,2e-5,1e-7,48),(.005,12,5e-6,2.5e-8,96),(.0025,16,1e-6,1e-8,192)):
  c=replace(base,energy_step=step,energy_order=order,energy_rtol=rtol,energy_atol=atol,fermi_poles=poles,energy_max_nodes=240000)
  cls=JunctionCompatibleBuilder if c.junction_basis=='nested' else ModeBuilder
  builders=[cls(c,g,v) for v in Frame(c.transport_axis).valleys()]
  b=create_transport_backend('negf',c);q=b.solve(state['phi_v'],.01,{'builders':builders,'grid':g,'fixed_mu':float(state['mu_l_ev']),'bias':.01,'detail':False})
  densities.append(q['n']);row={'step_ev':step,'order':order,'rtol':rtol,'atol':atol,'fermi_poles':poles,'id_a':q['id_a'],'gf_nodes':[d['real_energy_nodes'] for d in q['valleys']]};rows.append(row);print(row,flush=True)
 for i,r in enumerate(rows):
  r['density_relative_l2_to_finest']=float(np.linalg.norm(densities[i]-densities[-1])/np.linalg.norm(densities[-1]));r['current_relative_to_finest']=float(abs(r['id_a']/rows[-1]['id_a']-1))
 summary={'scope':'frozen final uniform biased field; not a re-solved SCF or mode/grid/contact certificate','cases':rows,'threshold':.01,'pass':all(r['density_relative_l2_to_finest']<.01 and r['current_relative_to_finest']<.01 for r in rows)}
 (Path(__file__).parent/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
