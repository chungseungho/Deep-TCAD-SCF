"""Frozen STEP08 electronic evaluator for same-runtime refactoring regression.

Copied unchanged from restored STEP08. Not imported by production code.
The historical NPZ golden is BLAS-sensitive for a truncated degenerate MS basis.
"""
import numpy as np
from device.regions import mode_regions
from modespace.hamiltonian import density_from_matrix
from multispace.occupation import solve_states
from multispace.constraints import solve_number

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

