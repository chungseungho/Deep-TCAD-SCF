"""Separate lead DOS/neutral mu, fixed-phi mode projection and grid experiments."""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from dataclasses import replace
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from config import Config
from modespace.hamiltonian import Grid
from device.definition import DeviceDefinition
from transport.contact_model import solve_reference
from multispace.occupation import fermi

def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);c=Config(backend='negf',contact_model='transverse_sp');rows=[];grids=[]
    for spacing in (1.,.5):
        cc=replace(c,dy=spacing,dz=spacing);g=Grid.build(cc);v=DeviceDefinition().frame.valleys();r=solve_reference(cc,g,v)
        grids.append({'dy_nm':spacing,'dz_nm':spacing,**r.diagnostics});np.savez_compressed(out/f'lead_grid_{spacing:g}.npz',y_nm=r.y,z_nm=r.z,phi_v=r.phi,density_cm3=r.density,mu_ev=r.mu)
        if spacing==1.:
            x,w=leggauss(512);k=(x+1)*np.pi/2;w=w/2;target=c.donor_sd*c.sd_width*c.sd_thickness*c.dx*1e-21
            for count in (4,8,16,24,32,48,64,100,150,300):
                def pop(mu):return [2*vv.multiplicity*(fermi(r.spectra[vv.name][0][:count,None]+2*r.spectra[vv.name][2]*(1-np.cos(k)),mu,c.temperature)@w) for vv in v]
                mu=brentq(lambda m:sum(np.sum(p) for p in pop(m))-target,-2,2)
                n=sum((abs(r.spectra[vv.name][1][:,:count])**2@o for vv,o in zip(v,pop(r.mu))),np.zeros(len(r.phi)))/(c.dx*c.dy*c.dz)*1e21
                rows.append({'modes':count,'fixed_full_mu_ev':r.mu,'neutral_truncated_mu_ev':mu,'neutral_mu_error_ev':abs(mu-r.mu),'fixed_mu_density_relative_l2':float(np.linalg.norm(n-r.density)/np.linalg.norm(r.density)), 'fixed_mu_missing_charge':float(1-n.mean()/c.donor_sd)})
    report={'lead_basis_at_fixed_SP_profile':rows,'grid_fully_selfconsistent_lead':grids,
        'grid_mu_change_ev':abs(grids[1]['zero_reference_mu_ev']-grids[0]['zero_reference_mu_ev']),
        'scope':'2D lead reference; not fully SCF device refinement'}
    (out/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':run(sys.argv[1])
