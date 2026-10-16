"""Three separate quantum density paths and physical Poisson fixed points."""
from pathlib import Path
import csv,json,time,os
import numpy as np
from config import KB_EV,Q
from device.regions import mode_regions,diagnose
from modespace.hamiltonian import density_from_matrix
from multispace.occupation import solve_states
from multispace.constraints import solve_number
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


def newton_physical_map(c,p,mapping,builders,grid,fixed_mu,vg,vd,out,key,phi_start,start,tic):
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
        q=electronic(c,builders,grid,mapping.potential(x),fixed_mu,vd)
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
        try:qnew=electronic(c,builders,grid,phi,fixed_mu,vd)
        except (BoundOccupationError,AmbiguousOccupation) as e:
            event({'kind':'physical_image','iteration':outer,'valid':False,'reason':str(e),
                   'bound_audit':getattr(e,'bound_diagnostics',None)})
            return
        raw=mapping.charge(qnew['n'])
        nr=float(np.sqrt(np.dot(p.volume,(raw-n)**2)/max(np.dot(p.volume,raw*raw),1e-30)))
        p.set_charge(raw);pr=float(max(abs(p.solve(vg,vd)-new)))
        stable,track=state_stability(q,qnew,c)
        ir=abs(qnew['id_a']-q['id_a'])/max(abs(qnew['id_a']),1e-15) if c.backend!='ms' else 0.
        row={'key':key,'iteration':start+outer,'vg_v':vg,'vd_v':vd,'phi_residual_v':pr,
            'density_relative_residual':nr,'mu_l_ev':qnew['mu'],'mu_change_ev':abs(qnew['mu']-q['mu']),
            'occupation_max_change':float(max(abs(qnew['tracking']-q['tracking']))),
            'classification_stable':stable,'electrons':float(qnew['n'].sum()*grid.cell_volume*1e-21),
            'id_a':qnew['id_a'],'elapsed_s':time.perf_counter()-tic,'current_relative_change':ir}
        append_csv(Path(out)/'scf_history.csv',row)
        accepted=bool(pr<c.phi_tol and nr<c.charge_tol and stable and row['mu_change_ev']<c.phi_tol and ir<c.current_tol)
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


def state_stability(old,new,c):
    if c.backend!='ms':return True,[]
    from multispace.tracking import match_states
    from scipy.linalg import block_diag
    stable=True;reports=[]
    for a,b,sa,sb in zip(old['states'],new['states'],old['systems'],new['systems']):
        overlaps=[ua.reshape(-1,ua.shape[-1]).T@ub.reshape(-1,ub.shape[-1]) for ua,ub in zip(sa.bases,sb.bases)]
        match=match_states(a,b,c.ms_degenerate_ev,basis_overlap=block_diag(*overlaps))
        i,j=match['old_indices'],match['new_indices'];active=a['active'][i]|b['active'][j]
        stable=stable and bool(np.all(a['origin'][i][active]==b['origin'][j][active]))
        reports.append(match)
    return stable,reports


def electronic(c,builders,grid,phi,fixed_mu,vd,detail=False):
    systems=[b.build(phi) for b in builders];deg=[2*s.valley.multiplicity for s in systems]
    density=np.zeros(len(grid.xyz));left=density.copy();right=density.copy();unresolved=density.copy()
    states=[];diagnostics=[];matrices=[];bonds=np.zeros(len(grid.x)-1);tracking=[]
    mu=float(fixed_mu);constraint={}
    if c.backend=='ms':
        for s in systems:
            labels=mode_regions(s,c);alternatives=[]
            step=c.ms_boundary_shift_slices*c.dx
            if step:
                for sign in (-1,1):alternatives.append(mode_regions(s,c,sign*step))
            states.append(solve_states(s.H,labels,c,alternatives))
        target=float(c.doping(grid.xyz).sum()*grid.cell_volume*1e-21)
        mu,states,constraint=solve_number(states,deg,target,vd,c,fixed_mu)
        for s,g,state in zip(systems,deg,states):
            u,f=state['u'],state['occupation'];rho=(u*f)@u.conj().T
            density+=density_from_matrix(s,rho,g)
            for origin,out in ((0,left),(2,right),(-1,unresolved)):
                use=state['origin']==origin
                part=(u[:,use]*f[use])@u[:,use].conj().T
                out+=density_from_matrix(s,part,g)
            ortho=float(np.max(abs(u.conj().T@u-np.eye(len(u)))))
            from transport.reservoir import bond_currents
            closed=bond_currents(s,rho,g)
            diagnostics.append({'valley':s.valley.name,'dimension':s.dimension,
                'orthogonality_error':ortho,'projector_sum_error':float(max(abs(state['weights'].sum(axis=1)-1))),
                'closed_stationary_max_current_a':float(max(abs(closed))),
                'pauli_min':float(min(f)),'pauli_max':float(max(f)),
                'origin_counts':{str(o):int(sum(state['origin']==o)) for o in (-1,0,2)},
                'active_states':int(sum(state['active'])),'degenerate_clusters':state['clusters'],
                'gf_evaluations_inside_scf':0})
            tracking.extend(f)
            if detail:matrices.append(rho)
    elif c.backend=='negf':
        from transport.negf import solve_negf,blocks_density
        for s,g in zip(systems,deg):
            blocks,flux,diag,components=solve_negf(s,c,mu,mu-vd,vd,detail)
            density+=blocks_density(s,blocks,g);bonds+=g*flux
            tracking.extend(np.concatenate([np.linalg.eigvalsh(x) for x in blocks]))
            diag['valley']=s.valley.name;diagnostics.append(diag)
            if detail:matrices.append(blocks)
    elif c.backend=='spectral':
        from transport.finite_spectral import solve_finite_spectral
        from transport.reservoir import bond_currents
        for s,g in zip(systems,deg):
            rho,diag,parts=solve_finite_spectral(s.H,s.regions,c.gamma,mu,mu-vd,c.temperature,
                                               split=detail,diagnostic_detail=detail)
            tracking.extend(np.linalg.eigvalsh(rho) if detail else diag.pop('_correlation_spectrum'))
            density+=density_from_matrix(s,rho,g);bonds+=bond_currents(s,rho,g)
            if detail:
                left+=density_from_matrix(s,parts[0],g);right+=density_from_matrix(s,parts[1],g)
                matrices.append(rho)
            diag['valley']=s.valley.name;diagnostics.append(diag)
    else:raise ValueError('Unsupported DEVSIM-QT backend')
    current=float('nan') if c.backend=='ms' else float(bonds[(grid.x[:-1]>=c.ls)&(grid.x[1:]<=c.drain_start)].mean())
    return {'n':density,'systems':systems,'states':states,'matrices':matrices,'mu':mu,
        'tracking':np.asarray(tracking),'bonds':bonds,'id_a':current,'constraint':constraint,
        'density_l':left,'density_r':right,'density_unresolved':unresolved,'valleys':diagnostics}


def scf(c,p,mapping,builders,grid,fixed_mu,vg,vd,out,key,initial=None):
    out=Path(out);tic=time.perf_counter();n=p.donors.copy() if initial is None else initial.copy()
    mixer=Anderson(p.volume,c.mix,c.history);prev_mu=None;prev_track=None;prev_sig=None;prev_q=None;prev_current=None
    accepted=False;rows=[];phi_initial=None
    for iteration in range(1,c.max_iterations+1):
        p.set_charge(n);phi_node=p.solve(vg,vd);phi=mapping.potential(phi_node)
        if phi_initial is None:phi_initial=phi.copy()
        q=electronic(c,builders,grid,phi,fixed_mu,vd)
        raw=mapping.charge(q['n']);res=raw-n
        nr=float(np.sqrt(np.dot(p.volume,res*res)/max(np.dot(p.volume,raw*raw),1e-30)))
        p.set_charge(raw);new_phi=p.solve(vg,vd);pr=float(max(abs(new_phi-phi_node)))
        mr=float(abs(q['mu']-prev_mu)) if prev_mu is not None else float('inf')
        tr=float(max(abs(q['tracking']-prev_track))) if prev_track is not None else float('inf')
        signature=tuple(tuple(s['origin'][s['active']]) for s in q['states'])
        classification_stable=prev_q is not None and state_stability(prev_q,q,c)[0]
        ir=abs(q['id_a']-prev_current)/max(abs(q['id_a']),1e-15) if prev_current is not None and c.backend!='ms' else (0. if c.backend=='ms' else float('inf'))
        if c.backend!='ms':classification_stable=True
        row={'key':key,'iteration':iteration,'vg_v':vg,'vd_v':vd,'phi_residual_v':pr,
             'density_relative_residual':nr,'mu_l_ev':q['mu'],'mu_change_ev':mr,
             'occupation_max_change':tr,'classification_stable':classification_stable,
             'electrons':float(q['n'].sum()*grid.cell_volume*1e-21),'id_a':q['id_a'],
             'elapsed_s':time.perf_counter()-tic,'current_relative_change':ir}
        append_csv(out/'scf_history.csv',row);rows.append(row)
        print(f'{key} {iteration:03d}: phiR={pr:.3e} nR={nr:.3e} mu={q["mu"]:.6f} stable={classification_stable}',flush=True)
        accepted=bool(iteration>=3 and pr<c.phi_tol and nr<c.charge_tol and
                      mr<c.phi_tol and classification_stable and ir<c.current_tol)
        atomic_npz(out/'checkpoints'/f'{key}.npz',n=n,phi_node_v=phi_node,vg_v=vg,vd_v=vd,
                   iteration=iteration,source_sha256=source_digest(),config_json=json.dumps(c.to_dict()))
        if accepted:break
        if iteration>=c.newton_start and pr>10*c.phi_tol and iteration<c.max_iterations-4:
            print(f'{key}: switching to Newton-Krylov on original physical map',flush=True)
            refined=newton_physical_map(c,p,mapping,builders,grid,fixed_mu,vg,vd,out,key,phi_node,iteration,tic)
            if refined is not None:n,phi_node,phi,q,accepted,iteration=refined
            break
        target=raw
        if c.preconditioner=='screened':
            chi=raw/(KB_EV*c.temperature);screened=p.solve(vg,vd,screening=(phi_node,chi))
            target=np.maximum(n+res+chi*(screened-phi_node),0.)
        nnext,method,gain=mixer.update(n/c.donor_sd,target/c.donor_sd)
        # Do not replace n with an unused mixed proposal after the last solve.
        if iteration<c.max_iterations:n=nnext*c.donor_sd
        prev_mu=q['mu'];prev_track=q['tracking'];prev_sig=signature;prev_q=q;prev_current=q['id_a']
    # Exactly reevaluate the saved H and raw residual before deciding acceptance.
    final=electronic(c,builders,grid,phi,fixed_mu,vd,True)
    raw=mapping.charge(final['n']);nr=float(np.sqrt(np.dot(p.volume,(raw-n)**2)/max(np.dot(p.volume,raw*raw),1e-30)))
    p.set_charge(raw);pr=float(max(abs(p.solve(vg,vd)-phi_node)))
    accepted=bool(accepted and nr<c.charge_tol and pr<c.phi_tol)
    p.set_charge(n);restored=p.solve(vg,vd)
    if max(abs(restored-phi_node))>1e-10:raise RuntimeError('Saved Poisson field did not restore')
    gauss=p.gauss_balance();accepted=accepted and gauss['gauss_relative_error']<1e-7
    region=diagnose(grid,final['n'],phi,c)
    from analysis.validity import element_fields,evaluate
    validity=evaluate(c,final['systems'],final['mu'],vd,element_fields(p.mesh_path,p))
    post=[];post_error=None
    if accepted and c.backend=='ms':
        from transport.ms_postprocess import postprocess
        try:
            post=[postprocess(s,c,final['mu'],vd) for s in final['systems']]
            final['id_a']=sum(x['current_a'] for x in post)
            final['bonds']=sum((x['bonds_a'] for x in post),np.zeros(len(grid.x)-1))
        except Exception as exc:post_error=repr(exc)
    current=final['id_a'];bonds=final['bonds']
    central=(grid.x[:-1]>=c.ls)&(grid.x[1:]<=c.drain_start)
    if c.backend=='ms':defect=max(abs(bonds-current)) if post else float('nan')
    else:defect=max(abs(bonds[central]-current))
    current_ok=bool(np.isfinite(current) and defect<max(1e-15,abs(current)*1e-6))
    result={'key':key,'backend':c.backend,'vg_v':vg,'vd_v':vd,'id_a':current,
        'numerical_status':'SCF_CONVERGED' if accepted else 'SCF_NOT_CONVERGED',
        'transport_status':'POST_SCF_FAILED' if post_error else 'CALCULATED' if accepted and current_ok else 'UNVALIDATED',
        'physical_status':validity['physical_status'],'iterations':iteration,
        'phi_residual_v':pr,'density_relative_residual':nr,'mu_l_ev':final['mu'],'mu_r_ev':final['mu']-vd,
        'current_conservation_abs_a':float(defect),'current_conservation_relative':float(defect/max(abs(current),1e-15)),
        'gauss_relative_error':gauss['gauss_relative_error'],
        'reservoir_screening_pass':region['reservoir_screening_diagnostic_pass'],
        'quantitative_validation_pass':False,'gf_inside_scf':c.backend=='negf',
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
    if c.backend in ('ms','spectral'):
        data.update(density_from_l_cm3=final['density_l'],density_from_r_cm3=final['density_r'],
                    density_unresolved_origin_cm3=final['density_unresolved'])
    for i,s in enumerate(final['systems']):
        data[f'valley{i}_subbands_ev']=np.stack([a[:min(c.modes,len(a))] for a in s.subbands])
        if c.backend=='ms':
            for name in ('e','weights','origin','largest_region','occupation','ambiguous','reach','contrast','boundary_unstable','cluster_unstable'):
                data[f'valley{i}_{name}']=final['states'][i][name]
    file=out/'states'/f'{key}.npz';atomic_npz(file,**data)
    write_json(file.with_suffix('.json'),{**result,'regions':region,'constraint':final['constraint'],
        'valleys':final['valleys'],'validity':validity,'gauss':gauss,
        'postprocess':post,'postprocess_error':post_error,
        'density_current_distinction':'closed stationary MS density; independent open scattering current after SCF'})
    if accepted:
        from export_fields import export_state
        export_state(file,out/'fields'/f'{key}_nm.vtu')
        p.export(out/'fields'/f'{key}_electrostatics_cm.vtk')
    return result,n,final
