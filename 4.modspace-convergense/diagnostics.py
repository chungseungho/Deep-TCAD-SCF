#!/usr/bin/env python3
"""Offline manufactured-field transfer + exact dense NEGF lead reference.

Works without DEVSIM and Gmsh; cannot assert real SCF convergence.
"""
import itertools
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from transfer import Transfer
from quantum import lead_matrix_self_energy, spectral_batched


def manufactured_cube():
    points=np.array(list(itertools.product([0.,1.],repeat=3)))
    by={tuple(p):i for i,p in enumerate(points)}
    tets=[]
    for perm in itertools.permutations(range(3)):
        v=np.zeros(3);q=[by[(0.,0.,0.)]]
        for idx in perm:
            v[idx]=1.;q.append(by[tuple(v)])
        tets.append(q)
    axes=(np.linspace(0,1,17),np.linspace(.025,.975,20),np.linspace(.025,.975,20))
    return Transfer.from_tetra(points,np.array(tets),points,axes)


def direct_reference(onsite,t,lead_l,lead_r,E,eta):
    """Full dense inverse; validates retarded G blocks and spectral density."""
    n,m,_=onsite.shape
    blocks=np.zeros((n*m,n*m),complex)
    for i in range(n):
        blocks[i*m:(i+1)*m,i*m:(i+1)*m]=onsite[i]
        if i+1<n:
            blocks[i*m:(i+1)*m,(i+1)*m:(i+2)*m]=-t*np.eye(m)
            blocks[(i+1)*m:(i+2)*m,i*m:(i+1)*m]=-t*np.eye(m)
    out=[];als=[];ars=[]
    for en in E:
        sl,gl=lead_matrix_self_energy([en],lead_l,t,eta)
        sr,gr=lead_matrix_self_energy([en],lead_r,t,eta)
        A=(en+1j*eta)*np.eye(n*m)-blocks
        A[:m,:m]-=sl[0];A[-m:,-m:]-=sr[0]
        G=np.linalg.inv(A)
        G0N=G[:m,-m:]
        out.append(np.trace(gl[0]@G0N@gr[0]@G0N.conj().T).real)
        l=[];r=[]
        for j in range(n):
            gi0=G[j*m:(j+1)*m,:m];giN=G[j*m:(j+1)*m,-m:]
            l.append(gi0@gl[0]@gi0.conj().T)
            r.append(giN@gr[0]@giN.conj().T)
        als.append(l);ars.append(r)
    return np.array(out),np.array(als),np.array(ars)


def run_diagnostics(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    tr=manufactured_cube()
    diag=tr.diagnostics()
    x,y,z=tr.axes; X,Y,Z=np.meshgrid(x,y,z,indexing='ij')
    p=1+.3*tr.nodes_nm[:,0]-.1*tr.nodes_nm[:,1]+.12*tr.nodes_nm[:,2]
    exact=1+.3*X-.1*Y+.12*Z
    got=tr.to_regular(p)
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
    ax[0].imshow(got[:, :, len(z)//2].T,origin='lower',aspect='auto',extent=[0,1,0,1]);ax[0].set_title('Exact tetra affine interpolation (synthetic field)')
    ax[1].imshow((got-exact)[:,:,len(z)//2].T,origin='lower',aspect='auto',extent=[0,1,0,1]);ax[1].set_title('Potential error (V)')
    for a in ax:a.set(xlabel='x',ylabel='y')
    fig.savefig(out/'01_tetra_to_fd.png',dpi=160);plt.close(fig)
    # Nonuniform (positive) charge tests the P.T volume-weighted transfer.
    n=(1+X+Y*Y+Z)*1e17
    nn=tr.to_nodes(n)
    Q=tr.totals(n,nn)
    diag.update(Q)
    fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    axs[0].bar(['FD charge integral','Tetra FE-lumped integral'],[Q['grid_particles'],Q['tet_particles']])
    axs[0].set(title='Geometric total-charge conservation',ylabel='Carrier number')
    axs[1].scatter(tr.nodes_nm[:,0],nn,c=tr.nodes_nm[:,1],s=40)
    axs[1].set(xlabel='x',ylabel='Mapped node density (cm$^{-3}$)',title='Adjoint density transfer (synthetic)')
    fig.savefig(out/'02_charge_conservation.png',dpi=160);plt.close(fig)
    nsite=9;m=3;t=.4;E=np.linspace(.15,1.2,32);eta=2e-5
    H0=np.diag([.25,.45,.7]).astype(complex)+2*t*np.eye(m)
    H0[0,1]=H0[1,0]=.075
    H0[1,2]=H0[2,1]=-.045
    onsite=np.tile(H0,(nsite,1,1))
    onsite[nsite//2,0,0]+=.09
    onsite[nsite//2,0,2]=onsite[nsite//2,2,0]=.03
    left=H0.copy();right=H0.copy()
    T,AL,AR=spectral_batched(onsite,t,left,right,E,eta,True)
    td,ald,ard=direct_reference(onsite,t,left,right,E,eta)
    tmax=float(np.max(abs(T-td)))
    amax=float(max(np.max(abs(AL-ald)),np.max(abs(AR-ard))))
    # Visualize real/imag parts and positive broadening of coupled leads.
    el=np.linspace(-.2,2.0,300)
    Sigma,Gamma=lead_matrix_self_energy(el,H0,t,eta)
    ev,U=np.linalg.eigh(H0)
    Sd=np.einsum('am,eab,bn->emn',U.conj(),Sigma,U,optimize=True)
    fig,axs=plt.subplots(1,2,figsize=(11,3.9),layout='constrained')
    for m in range(3):
        axs[0].plot(el,Sd[:,m,m].real,label=f'lead eigenmode {m+1}')
        axs[1].plot(el,(-2*Sd[:,m,m].imag).real,label=f'lead eigenmode {m+1}')
    axs[0].set(xlabel='Energy (eV)',ylabel='Re Sigma (eV)',title='Retarded lead energy shift')
    axs[1].set(xlabel='Energy (eV)',ylabel='Gamma (eV)',title='Contact broadening (nonnegative)')
    for a in axs:a.legend(fontsize=8);a.grid(alpha=.25)
    fig.savefig(out/'04_lead_selfenergy.png',dpi=160);plt.close(fig)
    diag.update({'retarded_selfenergy_dense_T_max_error':tmax,
                 'dense_injectivity_spectral_max_error':amax})
    fig,ax=plt.subplots(figsize=(7,4.2))
    ax.plot(E,T,'-',lw=2,label='Recursive GF, full-matrix leads')
    ax.plot(E,td,'o',fillstyle='none',ms=4,label='Dense full inverse reference')
    ax.set(xlabel='Energy (eV)',ylabel='Transmission',title='Contact self-energy / RGF validation')
    ax.legend();ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(out/'03_selfenergy_dense_reference.png',dpi=160);plt.close(fig)
    diag['synthetic_tests_pass']=bool(diag['linear_phi_error_V']<1e-10 and
      diag['relative_charge_mismatch']<1e-12 and tmax<1e-8 and amax<1e-6)
    (out/'diagnostics.json').write_text(json.dumps(diag,indent=2))
    (out/'STATUS.md').write_text('# STEP04 — Numerical diagnostics\n\n'
      '**Only manufactured field and small exact NEGF tests. Not a real-device SCF pass.**\n\n'
      +f"- PASS: {diag['synthetic_tests_pass']}\n- tetra interpolation max error: {diag['linear_phi_error_V']:.3e} V\n"
      +f"- Geometric discrete charge mismatch: {diag['relative_charge_mismatch']:.3e}\n"
      +f"- RGF vs full inverse max T error: {tmax:.3e}\n"
      +f"- RGF spectral density matrix max error: {amax:.3e}\n\n"
      + '\n'.join(f'![{f.name}]({f.name})' for f in out.glob('*.png'))+'\n')
    print(json.dumps(diag,indent=2))
    if not diag['synthetic_tests_pass']:
        raise SystemExit('FAILED numerical diagnostics')
    return diag


if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='figures')
    a=ap.parse_args();run_diagnostics(a.out)
