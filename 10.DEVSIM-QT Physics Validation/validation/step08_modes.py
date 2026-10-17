"""Frozen actual GAA M1-M5+, fixed-mu and fixed-N; reuse one eigensolve per basis."""
import sys,time,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from dataclasses import replace
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
from config import Config,VALLEYS,KB_EV
from modespace.hamiltonian import Grid
from modespace.junction_basis import JunctionCompatibleBuilder
from multispace.occupation import fermi,classify,occupations
from multispace.constraints import solve_number
from multispace.tracking import match_states
from device.regions import mode_regions
from transport.ms_postprocess import transmission,postprocess
from transport.negf_reference import solve_negf,blocks_density
from solver import write_json,atomic_npz
from provenance import source_digest


def density_eigen(s,u,f):
    out=[]
    for i,mask in enumerate(s.grid.mask):
        a=s.bases[i][mask]@u[s.sl(i)]
        out.append((abs(a)**2@f)*4/s.grid.cell_volume*1e21)
    return np.concatenate(out)


def run(source,out,levels=(2,4,6,8,10)):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);saved=np.load(source)
    base=Config(**json.loads(str(saved['config_json'])));grid=Grid.build(base);phi=saved['phi_v'];mu=float(saved['mu_l_ev'])
    energies=np.linspace(mu-.15,mu+.4,221);target=float(base.doping(grid.xyz).sum()*grid.cell_volume*1e-21)
    x=grid.xyz[:,0];regions={'channel':(x>=base.gate_start)&(x<=base.gate_end),
        'left':x<base.reservoir_left_nm,'right':x>base.length-base.reservoir_right_nm,
        'extension':((x>=base.reservoir_left_nm)&(x<base.gate_start))|((x>base.gate_end)&(x<=base.length-base.reservoir_right_nm))}
    rows=[];previous={}
    for m in levels:
        tic=time.time();c=replace(base,modes=m,sd_modes=2*m);systems=[];states=[]
        for v in VALLEYS:
            s=JunctionCompatibleBuilder(c,grid,v).build(phi);e,u=eigh(s.H,check_finite=False)
            systems.append(s);states.append({'e':e,'u':u,'origin':np.zeros(len(e),int),'uncertain_origin':np.zeros(len(e),bool)})
        muN=brentq(lambda mm:sum(4*fermi(s['e'],mm,c.temperature).sum() for s in states)-target,-3,3)
        T=sum((4*transmission(s,c,energies,0) for s in systems),np.zeros(len(energies)))
        # Linear response conductance uses equilibrium T; no fictitious biased MS density.
        conductance=float(np.trapezoid(T*fermi(energies,mu,300)*(1-fermi(energies,mu,300))/(KB_EV*300),energies)*3.8740458649318244e-5)
        for prescription,mm in [('fixed_mu',mu),('fixed_N',muN)]:
            n=sum((density_eigen(s,t['u'],fermi(t['e'],mm,c.temperature)) for s,t in zip(systems,states)),np.zeros(len(grid.xyz)))
            open_n=np.zeros(len(n))
            for s in systems:
                blocks,_,_,_=solve_negf(s,c,mm,mm,0)
                open_n+=blocks_density(s,blocks,4)
            row={'modes':m,'sd_modes':2*m,'prescription':prescription,'mu_ev':mm,
                 'electrons':float(n.sum()*grid.cell_volume*1e-21),'open_electrons':float(open_n.sum()*grid.cell_volume*1e-21),
                 'closed_open_relative_l2':float(np.linalg.norm(n-open_n)/np.linalg.norm(open_n)),
                 'dimensions':[s.dimension for s in systems],
                 'channel_subband_minima_ev':[[float(z) for z in np.min(np.stack(s.subbands[len(grid.x)//2:len(grid.x)//2+1]),axis=0)[:5]] for s in systems],
                 'thermal_occupied_eigenstate_counts':[int(sum(fermi(t['e'],mm,c.temperature)>1e-7)) for t in states],
                 'linear_response_conductance_fixed_contact_mu_s':conductance,
                 'scf':False,'regional':{},'runtime_s':time.time()-tic}
            for name,mask in regions.items():
                row['regional'][name]={'electrons':float(n[mask].sum()*grid.cell_volume*1e-21),
                    'open_electrons':float(open_n[mask].sum()*grid.cell_volume*1e-21),
                    'closed_open_relative_l2':float(np.linalg.norm(n[mask]-open_n[mask])/max(np.linalg.norm(open_n[mask]),1e-300))}
            if prescription in previous:
                pn,pT=previous[prescription]
                row['density_change_from_previous']=float(np.linalg.norm(n-pn)/np.linalg.norm(n))
                row['channel_charge_change_from_previous']=float(abs(n[regions['channel']].sum()-pn[regions['channel']].sum())/n[regions['channel']].sum())
                row['transmission_relative_l2_change']=float(np.linalg.norm(T-pT)/np.linalg.norm(T))
            previous[prescription]=(n,T)
            row['basis_orthogonality_max']=max(float(np.max(abs(s.bases[i][mask].T@s.bases[i][mask]-np.eye(s.bases[i].shape[-1])))) for s in systems for i,mask in enumerate(grid.mask))
            rows.append(row);atomic_npz(out/f'M{m}_{prescription}.npz',xyz_nm=grid.xyz,density_cm3=n,open_density_cm3=open_n,
                energy_ev=energies,transmission=T,mu_ev=mm,phi_v=phi,config_json=json.dumps(c.to_dict()),source_sha256=source_digest())
        write_json(out/'summary.json',{'source_sha256':source_digest(),'source_state':str(source),'cases':rows,'complete':m==levels[-1],
            'plateau_certified':False,'scope':'frozen actual GAA field; T at equilibrium, linear-response integration diagnostic only'})
        print('GAA modes',m,'seconds',time.time()-tic,flush=True)

if __name__=='__main__':run(sys.argv[1],sys.argv[2],tuple(map(int,sys.argv[3].split(','))) if len(sys.argv)>3 else (2,4,6,8,10))
