"""Three separate quantum density paths and physical Poisson fixed points."""
from pathlib import Path
import csv,json,time,os
import numpy as np
from config import KB_EV,Q
from device.regions import diagnose
from coupling.mixing import Anderson
from provenance import source_digest,fingerprint


def atomic_npz(path,**values):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.stem+'.tmp.npz')
    np.savez_compressed(tmp,**values);os.replace(tmp,path)


def write_json(path,value):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=True,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item()))


def append_csv(path,row):
    path=Path(path);exists=path.exists()
    with path.open('a',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(row))
        if not exists:w.writeheader()
        w.writerow(row)


def newton_physical_map(c,p,mapping,backend,builders,grid,fixed_mu,vg,vd,out,key,phi_start,start,tic):
    """Admissible-domain Newton on the original physical Poisson map."""
    from coupling.domain_newton import solve,DomainSearchFailure
    from transport.bound_states import BoundOccupationError
    from multispace.occupation import AmbiguousOccupation
    cache={};last=None;outer=0
    def event(row):
        with (Path(out)/(key+'_newton_trials.jsonl')).open('a') as f:
            f.write(json.dumps(row,default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item())+'\n')
    def evaluate(x):
        if 'x' in cache and np.array_equal(x,cache['x']):return cache['value']
        q=electronic(backend,builders,grid,mapping.potential(x),fixed_mu,vd)
        density=mapping.charge(q['n']);p.set_charge(density);new=p.solve(vg,vd)
        value=(new-x,new,density,q);cache.update(x=x.copy(),value=value);return value
    def callback(x,residual):
        nonlocal outer,last
        outer+=1;res,new,n,q=evaluate(x);phi=mapping.potential(new)
        atomic_npz(Path(out)/'trials'/f'{key}_accepted{outer:03d}.npz',
            trial_node_phi_v=x,trial_density_cm3=q['n'],poisson_image_phi_v=new,
            residual_v=res,numerically_converged=False,vg_v=vg,vd_v=vd,
            config_json=json.dumps(c.to_dict()),source_sha256=source_digest())
        # A Poisson image is another trial. Invalid images are not final states.
        try:qnew=electronic(backend,builders,grid,phi,fixed_mu,vd)
        except (BoundOccupationError,AmbiguousOccupation) as e:
            event({'kind':'physical_image','iteration':outer,'valid':False,'reason':str(e),
                   'bound_audit':getattr(e,'bound_diagnostics',None)})
            return
        raw=mapping.charge(qnew['n'])
        nr=float(np.sqrt(np.dot(p.volume,(raw-n)**2)/max(np.dot(p.volume,raw*raw),1e-30)))
        p.set_charge(raw);pr=float(max(abs(p.solve(vg,vd)-new)))
        metrics=backend.convergence_metrics(q,qnew)
        stable=metrics['classification_stable']; ir=metrics['current_relative_change']
        row={'key':key,'iteration':start+outer,'vg_v':vg,'vd_v':vd,'phi_residual_v':pr,
            'density_relative_residual':nr,'mu_l_ev':qnew['mu'],'mu_change_ev':abs(qnew['mu']-q['mu']),
            'occupation_max_change':metrics['occupation_max_change'],
            'classification_stable':stable,'electrons':float(qnew['n'].sum()*grid.cell_volume*1e-21),
            'id_a':qnew['id_a'],'elapsed_s':time.perf_counter()-tic,'current_relative_change':ir}
        append_csv(Path(out)/'scf_history.csv',row)
        accepted=backend.converged(start+outer,metrics,pr,nr,row['mu_change_ev'])
        last=(n,new,phi,qnew,accepted,start+outer)
        atomic_npz(Path(out)/'checkpoints'/f'{key}.npz',n=n,phi_node_v=new,vg_v=vg,vd_v=vd,
                   iteration=start+outer,source_sha256=source_digest(),config_json=json.dumps(c.to_dict()))
        print(f'{key} Domain-Newton {outer:02d}: phiR={pr:.3e} nR={nr:.3e} currentR={ir:.3e}',flush=True)
    preconditioner=None
    if c.preconditioner=='screened':
        from coupling.mixing import screened_inverse
        density=evaluate(phi_start)[2];chi=np.zeros(len(p.allxyz))
        chi[p.indices['silicon']]=density/(KB_EV*c.temperature)
        preconditioner=screened_inverse(p.K,Q*chi*p.si_volume,p.indices['silicon'],p.boundary)
    try:
        root,info=solve(lambda x:evaluate(x)[0],phi_start,
            invalid=(BoundOccupationError,AmbiguousOccupation),preconditioner=preconditioner,
            maxiter=c.max_iterations-start,tol=c.phi_tol*.1,
            trust_radius=c.newton_trust_radius_v,max_backtracks=c.newton_max_backtracks,
            callback=callback,log=event)
        write_json(Path(out)/(key+'_newton_status.json'),{k:v for k,v in info.items() if k!='events'})
    except DomainSearchFailure as exc:
        event({'kind':'termination','status':'DOMAIN_SEARCH_STALLED','not_a_nonexistence_proof':True})
        if last is None:raise
    return last


def electronic(backend,builders,grid,phi,fixed_mu,vd,detail=False):
    """Execute one transport step through the replaceable backend contract."""
    context={'config':backend.config,'builders':builders,'grid':grid,
             'fixed_mu':fixed_mu,'detail':detail,'bias':vd}
    result=backend.solve(phi,vd,context)
    backend.validate(result)
    return result


def scf(c,p,mapping,builders,grid,fixed_mu,vg,vd,out,key,initial=None,backend=None):
    if backend is None:
        raise ValueError('SCF requires a backend instance from the transport registry')
    out=Path(out);tic=time.perf_counter();n=p.donors.copy() if initial is None else initial.copy()
    mixer=Anderson(p.volume,c.mix,c.history);prev_mu=None;prev_track=None;prev_sig=None;prev_q=None;prev_current=None
    accepted=False;rows=[];phi_initial=None
    backend.prepare({'grid':grid,'mapping':mapping})
    for iteration in range(1,c.max_iterations+1):
        from coupling.scf_controller import evaluate_transport_step
        step=evaluate_transport_step(p,mapping,backend,builders,grid,n,fixed_mu,vg,vd)
        phi_node=step.poisson_potential;phi=step.quantum_potential
        if phi_initial is None:phi_initial=phi.copy()
        q=step.transport_result
        raw=step.mapped_density;res=raw-n
        nr=float(np.sqrt(np.dot(p.volume,res*res)/max(np.dot(p.volume,raw*raw),1e-30)))
        new_phi=step.poisson_image;pr=step.potential_residual
        mr=float(abs(q['mu']-prev_mu)) if prev_mu is not None else float('inf')
        metrics=backend.convergence_metrics(prev_q,q)
        tr=metrics['occupation_max_change']
        classification_stable=metrics['classification_stable']
        ir=metrics['current_relative_change']
        row={'key':key,'iteration':iteration,'vg_v':vg,'vd_v':vd,'phi_residual_v':pr,
             'density_relative_residual':nr,'mu_l_ev':q['mu'],'mu_change_ev':mr,
             'occupation_max_change':tr,'classification_stable':classification_stable,
             'electrons':float(q['n'].sum()*grid.cell_volume*1e-21),'id_a':q['id_a'],
             'elapsed_s':time.perf_counter()-tic,'current_relative_change':ir}
        append_csv(out/'scf_history.csv',row);rows.append(row)
        print(f'{key} {iteration:03d}: phiR={pr:.3e} nR={nr:.3e} mu={q["mu"]:.6f} stable={classification_stable}',flush=True)
        accepted=backend.converged(iteration,metrics,pr,nr,mr)
        atomic_npz(out/'checkpoints'/f'{key}.npz',n=n,phi_node_v=phi_node,vg_v=vg,vd_v=vd,
                   iteration=iteration,source_sha256=source_digest(),config_json=json.dumps(c.to_dict()))
        if accepted:break
        if iteration>=c.newton_start and pr>10*c.phi_tol and iteration<c.max_iterations-4:
            print(f'{key}: switching to Newton-Krylov on original physical map',flush=True)
            refined=newton_physical_map(c,p,mapping,backend,builders,grid,fixed_mu,vg,vd,out,key,phi_node,iteration,tic)
            if refined is not None:n,phi_node,phi,q,accepted,iteration=refined
            break
        target=raw
        if c.preconditioner=='screened':
            chi=raw/(KB_EV*c.temperature);screened=p.solve(vg,vd,screening=(phi_node,chi))
            target=np.maximum(n+res+chi*(screened-phi_node),0.)
        nnext,method,gain=mixer.update(n/c.donor_sd,target/c.donor_sd)
        # Do not replace n with an unused mixed proposal after the last solve.
        if iteration<c.max_iterations:n=nnext*c.donor_sd
        prev_mu=q['mu'];prev_q=q
    # Exactly reevaluate the saved H and raw residual before deciding acceptance.
    final=electronic(backend,builders,grid,phi,fixed_mu,vd,True)
    raw=mapping.charge(final['n']);nr=float(np.sqrt(np.dot(p.volume,(raw-n)**2)/max(np.dot(p.volume,raw*raw),1e-30)))
    p.set_charge(raw);pr=float(max(abs(p.solve(vg,vd)-phi_node)))
    accepted=bool(accepted and nr<c.charge_tol and pr<c.phi_tol)
    p.set_charge(n);restored=p.solve(vg,vd)
    if max(abs(restored-phi_node))>1e-10:raise RuntimeError('Saved Poisson field did not restore')
    gauss=p.gauss_balance();accepted=accepted and gauss['gauss_relative_error']<1e-7
    region=diagnose(grid,final['n'],phi,c)
    from analysis.validity import element_fields,evaluate
    validity=evaluate(c,final['systems'],final['mu'],vd,element_fields(p.mesh_path,p))
    final['config']=c
    backend_final=backend.finalize(final,vd,grid,accepted)
    current=final['id_a'];bonds=final['bonds']
    defect=backend_final['current_conservation_abs_a']
    current_ok=backend_final['current_ok']
    result={'key':key,'backend':backend.name,'vg_v':vg,'vd_v':vd,'id_a':current,
        'numerical_status':'SCF_CONVERGED' if accepted else 'SCF_NOT_CONVERGED',
        'transport_status':backend_final['transport_status'],
        'physical_status':validity['physical_status'],'iterations':iteration,
        'phi_residual_v':pr,'density_relative_residual':nr,'mu_l_ev':final['mu'],'mu_r_ev':final['mu']-vd,
        'current_conservation_abs_a':float(defect),'current_conservation_relative':float(defect/max(abs(current),1e-15)),
        'gauss_relative_error':gauss['gauss_relative_error'],
        'reservoir_screening_pass':region['reservoir_screening_diagnostic_pass'],
        'quantitative_validation_pass':False,**backend.result_metadata(),
        'source_sha256':source_digest(),'configuration_fingerprint':fingerprint(c),'runtime_s':time.perf_counter()-tic}
    data={'xyz_nm':grid.xyz,'x_nm':grid.x,'y_nm':grid.y,'z_nm':grid.z,'silicon_mask':grid.mask,
        'phi_v':phi,'electron_density_cm3':final['n'],'poisson_density_cm3':n,
        'node_xyz_nm':p.xyz['silicon'],'node_phi_v':phi_node,'donors_cm3':p.donors,
        'node_volume_nm3':p.volume,'cell_volume_nm3':grid.cell_volume,
        'doping_quantum_cm3':c.doping(grid.xyz),'bond_current_a':bonds,
        'converged':accepted,'config_json':json.dumps(c.to_dict(),sort_keys=True),**result}
    from device.definition import Frame
    data.update(world_xyz_nm=Frame(c.transport_axis).to_world(grid.xyz),
                world_node_xyz_nm=Frame(c.transport_axis).to_world(p.xyz['silicon']),
                world_transport_axis=c.transport_axis,potential_reference_v=c.potential_reference_v)
    data.update(backend.export_fields(final))
    file=out/'states'/f'{key}.npz';atomic_npz(file,**data)
    write_json(file.with_suffix('.json'),{**result,'regions':region,'constraint':final['constraint'],
        'valleys':final['valleys'],'validity':validity,'gauss':gauss,
        'postprocess':backend_final['postprocess'],'postprocess_error':backend_final['postprocess_error'],
        **backend.result_metadata()})
    if accepted:
        from export_fields import export_state
        export_state(file,out/'fields'/f'{key}_nm.vtu')
        p.export(out/'fields'/f'{key}_electrostatics_cm.vtk')
    return result,n,final
