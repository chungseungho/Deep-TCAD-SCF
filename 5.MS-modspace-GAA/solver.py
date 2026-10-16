"""GF-free electronic density, explicit reservoir battery, and genuine SCF residuals."""
from pathlib import Path
import csv,json,time,os
import numpy as np
from config import KB_EV,KIN,VALLEYS
from modespace.hamiltonian import density_from_matrix
from multispace.occupation import partition_eigenstates,closed_ms_density
from transport.reservoir import reservoirs,solve_steady,bond_currents,contact_currents,reservoir_spacing
from transport.finite_spectral import solve_finite_spectral
from coupling.mixing import Anderson
from scipy.optimize import newton_krylov,NoConvergence


def atomic_npz(path,**values):
    path=Path(path);tmp=path.with_name(path.stem+'.tmp.npz')
    try:
        np.savez_compressed(tmp,**values);os.replace(tmp,path)
    finally: tmp.unlink(missing_ok=True)


def electronic(c,builders,grid,phi,mu,vd,detail=False):
    n=np.zeros(len(grid.xyz));bonds=np.zeros(len(grid.x)-1);contacts=np.zeros(2);parts=[]
    diagnostics=[];states=[];matrix=[];occupations=[];eigenvalues=[]
    for builder in builders:
        sys=builder.build(phi);g=2*builder.v.multiplicity
        baths=reservoirs(sys.H,sys.regions,c.gamma,mu,mu-vd,c.temperature)
        baths.regions=sys.regions
        if c.backend in ('spectral','ms'):
            driven,diag,components=solve_finite_spectral(sys.H,sys.regions,c.gamma,mu,mu-vd,c.temperature,split=detail)
        else: driven,diag,components=solve_steady(sys.H,baths,split=detail)
        if c.backend=='ms':
            rho,state=closed_ms_density(sys.H,sys.regions,mu,mu-vd,c.temperature,c.ms_tie_tol)
            # Transport postprocessing opens the converged H; record mismatch.
            diag['ms_driven_density_matrix_relative_difference']=float(np.linalg.norm(driven-rho)/max(np.linalg.norm(rho),1e-30))
        else:
            rho=driven
            state=partition_eigenstates(sys.H,sys.regions,mu,mu-vd,c.temperature,c.ms_tie_tol)
        state['driven_eigen_occupation']=np.diag(state['u'].T@driven@state['u']).real
        occupations.extend(state['occupation'] if c.backend=='ms' else state['driven_eigen_occupation'])
        eigenvalues.extend(state['e'])
        n+=density_from_matrix(sys,rho,g)
        cur=bond_currents(sys,driven,g)
        if c.backend in ('spectral','ms'):
            left=np.where(sys.grid.x<c.ls)[0][-1]
            right=np.where(sys.grid.x>c.drain_start)[0][0]-1
            con=[cur[left],-cur[right]]
        else: con=contact_currents(sys,driven,baths,g)
        bonds+=cur;contacts+=con
        dL=reservoir_spacing(baths.levels_l,mu,c.temperature)
        dR=reservoir_spacing(baths.levels_r,mu-vd,c.temperature)
        diag.update(valley=builder.v.name,dimension=sys.dimension,spacing_l_ev=dL,spacing_r_ev=dR,
                    gamma_ev=c.gamma,gamma_over_kbt=c.gamma/(KB_EV*c.temperature),
                    spacing_over_gamma=max(dL,dR)/c.gamma,
                    continuum_hierarchy_pass=bool(max(dL,dR)<c.gamma/3 and c.gamma<KB_EV*c.temperature/3))
        diagnostics.append(diag)
        if detail:
            states.append((sys,state));matrix.append((sys,driven,baths))
            parts.append((density_from_matrix(sys,components[0],g),density_from_matrix(sys,components[1],g)))
    # All bonds bounding the undriven C contain the same particle flux.
    central=(grid.x[:-1]>=c.ls-c.dx/2)&(grid.x[1:]<=c.drain_start+c.dx/2)
    id_a=float(np.mean(bonds[central]))
    conservation=max(abs(contacts.sum()),max(abs(bonds[central]-id_a)))
    diag={'valleys':diagnostics,'id_a':id_a,'source_a':float(contacts[0]),'drain_a':float(contacts[1]),
          'current_conservation_abs_a':float(conservation),
          'current_conservation_relative':float(conservation/max(abs(id_a),1e-15)),
          'continuum_hierarchy_pass':all(d['continuum_hierarchy_pass'] for d in diagnostics)}
    return {'n':n,'bonds':bonds,'diag':diag,'states':states,'matrices':matrix,'parts':parts,
            'occupation_tracking':np.array(occupations),'eigenvalue_tracking':np.array(eigenvalues)}


def append_csv(path,row):
    path=Path(path);exists=path.exists()
    with path.open('a',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(row));
        if not exists: w.writeheader()
        w.writerow(row)


def newton_refine(c,p,mapping,builders,grid,mu,vg,vd,out,key,last,tic):
    """Matrix-free Newton-Krylov on the ORIGINAL Poisson/quantum potential map.

    No fitted capacitance, broadened tolerance or extra physical density is
    introduced. Jacobian-vector products are finite differences of real
    quantum evaluations. Armijo line search controls nonlinear steps.
    """
    start=last[3]['iteration'];previous_current=last[3]['id_a'];counter=0;cache={}
    previous_occup=last[2]['occupation_tracking'];previous_energy=last[2]['eigenvalue_tracking']
    def evaluate(x):
        if 'x' in cache and np.array_equal(cache['x'],x): return cache['value']
        q=electronic(c,builders,grid,mapping.potential(x),mu,vd)
        n=mapping.charge(q['n']);p.set_charge(n);new=p.solve(vg,vd)
        value=(new-x,new,n,q)
        cache.update(x=x.copy(),value=value);return value
    def fun(x): return evaluate(x)[0]
    def callback(x,f):
        nonlocal counter,previous_current,previous_occup,previous_energy,last
        counter+=1;res,new,n,q=evaluate(x)
        qnew=electronic(c,builders,grid,mapping.potential(new),mu,vd)
        nn=mapping.charge(qnew['n']);nr=np.sqrt(np.dot(p.volume,(nn-n)**2)/max(np.dot(p.volume,nn*nn),1e-30))
        cur=qnew['diag']['id_a'];ir=abs(cur-q['diag']['id_a'])/max(abs(cur),1e-15)
        sd=(p.xyz['silicon'][:,0]<c.ls/2)|(p.xyz['silicon'][:,0]>c.length-c.ld/2)
        neutral=float(np.dot(p.volume[sd],nn[sd]-p.donors[sd])/max(np.dot(p.volume[sd],p.donors[sd]),1e-30))
        row=dict(key=key,iteration=start+counter,vg_v=vg,vd_v=vd,id_a=cur,
                 phi_residual_v=float(max(abs(res))),charge_residual=float(nr),current_change=float(ir),
                 electron_number=float(np.dot(qnew['n'],mapping.sitevol)*1e-21),
                 sd_neutrality_relative_error=neutral,
                 occupation_change=float(max(abs(qnew['occupation_tracking']-previous_occup))),
                 eigenvalue_change_ev=float(max(abs(qnew['eigenvalue_tracking']-previous_energy))),
                 mixing='newton_krylov',gain=1.,elapsed_s=time.perf_counter()-tic)
        append_csv(out/'scf_history.csv',row);previous_current=cur
        previous_occup=qnew['occupation_tracking'];previous_energy=qnew['eigenvalue_tracking']
        last=(new,mapping.potential(new),qnew,row,n)
        atomic_npz(out/'checkpoints'/f'{key}.npz',n=n,iteration=start+counter,vg=vg,vd=vd)
        print(f'[{key}] Newton {counter:02d} I={cur:.4e} phiR={row["phi_residual_v"]:.3e} nR={nr:.3e}',flush=True)
    remaining=max(1,min(30,c.max_iterations-start))
    try:
        root=newton_krylov(fun,last[0],rdiff=1e-5,inner_maxiter=15,outer_k=8,
                           maxiter=remaining,f_tol=c.phi_tol*.2,line_search='armijo',callback=callback)
        # Callback also runs for final accepted iterate; perform once if initial guess met tolerance.
        if counter==0: callback(root,fun(root))
    except NoConvergence as exc:
        if counter==0: callback(exc.args[0],fun(exc.args[0]))
    except Exception as exc:
        (out/f'{key}_newton_failure.txt').write_text(repr(exc))
        print('Newton fallback failed:',repr(exc),flush=True)
    row=last[3]
    accepted=bool(row['phi_residual_v']<c.phi_tol and row['charge_residual']<c.charge_tol and row['current_change']<c.current_tol)
    return last,accepted


def run_scf(c,p,mapping,builders,grid,mu,vg,vd,out,key,initial=None,resume=False):
    out=Path(out);checkpoint=out/'checkpoints'/f'{key}.npz';checkpoint.parent.mkdir(exist_ok=True)
    n=initial.copy() if initial is not None else p.donors.copy()
    start_iteration=0
    if resume and checkpoint.exists():
        cp=np.load(checkpoint)
        if float(cp['vg'])!=vg or float(cp['vd'])!=vd: raise ValueError('Checkpoint bias mismatch')
        # An exhausted nonconverged checkpoint supplies a restart guess, not
        # zero remaining iterations. Keep the previous failure in CSV/logs.
        n=cp['n'];start_iteration=int(cp['iteration'])
        if start_iteration>=c.max_iterations: start_iteration=0
    # One-shot ML input is independent of the converged solution (avoid target leakage).
    p.set_charge(p.donors)
    seed_phi=mapping.potential(p.solve(vg,vd))
    mixer=Anderson(p.volume,c.mix,c.history);tic=time.perf_counter();previous_current=None
    previous_occup=None;previous_energy=None;converged=False;last=None
    sd=(p.xyz['silicon'][:,0]<c.ls/2)|(p.xyz['silicon'][:,0]>c.length-c.ld/2)
    for it in range(start_iteration+1,c.max_iterations+1):
        p.set_charge(n);phi_node=p.solve(vg,vd);phi=mapping.potential(phi_node)
        q=electronic(c,builders,grid,phi,mu,vd)
        raw=mapping.charge(q['n']);r=raw-n
        charge_res=np.sqrt(np.dot(p.volume,r*r)/max(np.dot(p.volume,raw*raw),1e-30))
        # UNMIXED fixed-point Poisson residual (small damping cannot fake convergence).
        p.set_charge(raw);raw_phi_node=p.solve(vg,vd)
        phi_res=float(max(abs(raw_phi_node-phi_node)))
        cur=q['diag']['id_a'];current_res=abs(cur-previous_current)/max(abs(cur),1e-15) if previous_current is not None else np.inf
        number=float(np.dot(q['n'],mapping.sitevol)*1e-21)
        neutrality=float(np.dot(p.volume[sd],raw[sd]-p.donors[sd])/max(np.dot(p.volume[sd],p.donors[sd]),1e-30))
        # Population and spectrum tracking are basis-invariant traces/spectral lists.
        populations=q['occupation_tracking'];energies=q['eigenvalue_tracking']
        occ_change=float(max(abs(populations-previous_occup))) if previous_occup is not None else np.nan
        energy_change=float(max(abs(energies-previous_energy))) if previous_energy is not None else np.nan
        target=raw
        if c.preconditioner=='screened':
            # Approximate positive local quantum capacitance, not an extra
            # physical carrier model. Used ONLY to precondition the fixed point.
            chi=raw/(KB_EV*c.temperature)
            screened=p.solve(vg,vd,screening=(phi_node,chi))
            direction=r+chi*(screened-phi_node)
            target=np.maximum(n+direction,0)
        n_next,method,gain=mixer.update(n/c.donor_sd,target/c.donor_sd);n_next*=c.donor_sd
        if c.preconditioner=='screened': method='screened_'+method
        row=dict(key=key,iteration=it,vg_v=vg,vd_v=vd,id_a=cur,phi_residual_v=phi_res,
                 charge_residual=float(charge_res),current_change=float(current_res),electron_number=number,
                 sd_neutrality_relative_error=neutrality,occupation_change=occ_change,eigenvalue_change_ev=energy_change,
                 mixing=method,gain=gain,elapsed_s=time.perf_counter()-tic)
        append_csv(out/'scf_history.csv',row)
        print(f'[{key}] {it:03d} I={cur:.4e} phiR={phi_res:.3e} nR={charge_res:.3e} {method}',flush=True)
        converged=bool(it>=3 and charge_res<c.charge_tol and phi_res<c.phi_tol and current_res<c.current_tol)
        # Saved charge is always the input to the saved H/phi, never a stale output.
        atomic_npz(checkpoint,n=n if converged else n_next,iteration=it,vg=vg,vd=vd)
        last=(phi_node,phi,q,row,n.copy())
        if converged: break
        if it>=40 and it%20==0:
            # Switch algorithm on sustained residual, not on tiny mixed steps.
            # History remains on disk. The current actual residual triggers a
            # nonlinear solve early enough to retain the requested iteration cap.
            if phi_res>10*c.phi_tol and charge_res>c.charge_tol:
                print(f'[{key}] switching from Anderson to Newton-Krylov',flush=True)
                break
        n=n_next;previous_current=cur;previous_occup=populations;previous_energy=energies
    if last is None: raise RuntimeError('Checkpoint exhausted max_iterations; use a new run with a larger limit')
    if not converged and last[3]['iteration']<c.max_iterations:
        last,converged=newton_refine(c,p,mapping,builders,grid,mu,vg,vd,out,key,last,tic)
    phi_node,phi,q,row,n=last
    # Re-evaluate density/current on the exact saved potential with detailed diagnostics.
    final=electronic(c,builders,grid,phi,mu,vd,detail=True)
    raw=mapping.charge(final['n']);res=np.sqrt(np.dot(p.volume,(raw-n)**2)/max(np.dot(p.volume,raw*raw),1e-30))
    converged=bool(converged and res<c.charge_tol)
    current_ok=final['diag']['current_conservation_abs_a']<max(1e-15,abs(final['diag']['id_a'])*1e-6)
    result=dict(key=key,vg_v=vg,vd_v=vd,id_a=final['diag']['id_a'],
                status='SCF_CONVERGED_PROTOTYPE' if converged and current_ok else 'NOT_CONVERGED',
                iterations=row['iteration'],phi_residual_v=row['phi_residual_v'],charge_residual=float(res),
                sd_neutrality_relative_error=row['sd_neutrality_relative_error'],
                continuum_hierarchy_pass=final['diag']['continuum_hierarchy_pass'],
                current_conservation_relative=final['diag']['current_conservation_relative'],
                quantitative_validation_pass=False,backend=c.backend,runtime_s=time.perf_counter()-tic)
    bandwidth=min(4*KIN/(v.mx*c.dx**2) for v in VALLEYS)
    result['minimum_longitudinal_bandwidth_ev']=bandwidth
    result['bias_inside_minimum_bandwidth']=bool(abs(vd)+8*KB_EV*c.temperature<bandwidth)
    state=out/'states'/f'{key}.npz';state.parent.mkdir(exist_ok=True)
    data=dict(xyz_nm=grid.xyz,x_nm=grid.x,y_nm=grid.y,z_nm=grid.z,silicon_mask=grid.mask,
              phi_v=phi,initial_neutral_poisson_phi_v=seed_phi,
              electron_density_cm3=final['n'],poisson_density_cm3=n,
              node_xyz_nm=p.xyz['silicon'],node_phi_v=phi_node,donors_cm3=p.donors,
              node_volume_nm3=p.volume,cell_volume_nm3=grid.cell_volume,bond_current_a=final['bonds'],
              vg_v=vg,vd_v=vd,mu_l_ev=mu,mu_r_ev=mu-vd,
              converged=(result['status']=='SCF_CONVERGED_PROTOTYPE'),
              quantitative_validation_pass=False,backend=c.backend)
    if final['parts']:
        data['density_from_l_cm3']=sum(v[0] for v in final['parts'])
        data['density_from_r_cm3']=sum(v[1] for v in final['parts'])
    for j,(sys,s) in enumerate(final['states']):
        prefix=f'valley{j}_';data.update({prefix+k:s[k] for k in ['e','weights','origin','occupation','ambiguous','driven_eigen_occupation']})
        centre=len(grid.x)//2
        data[prefix+'subbands_ev']=np.stack([v[:min(c.modes,len(v))] for v in sys.subbands])
        data[prefix+'centre_modes']=sys.bases[centre]
    atomic_npz(state,**data)
    (out/'states'/f'{key}.json').write_text(json.dumps({**result,**final['diag']},indent=2))
    p.set_charge(n)
    p.solve(vg,vd)  # restore the exact unscreened physical field for VTK/export
    return result,n,final
