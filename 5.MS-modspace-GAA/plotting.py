"""Figures from saved CAD/mesh/solver data only. No invented I-V points."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

COLORS={'Silicon':'#5d9fd1','GateOxide':'#e6b254','Spacer':'#75bfb3'}


def parse_mesh(path):
    lines=Path(path).read_text().splitlines();i=lines.index('$PhysicalNames');n=int(lines[i+1]);names={}
    for row in lines[i+2:i+2+n]:
        dim,tag,name=row.split(maxsplit=2);names[(int(dim),int(tag))]=name.strip('"')
    i=lines.index('$Nodes');n=int(lines[i+1]);raw=[r.split() for r in lines[i+2:i+2+n]]
    byid={int(r[0]):j for j,r in enumerate(raw)};xyz=np.array([[float(v) for v in r[1:]] for r in raw])
    faces={};tets={};i=lines.index('$Elements');n=int(lines[i+1])
    for row in lines[i+2:i+2+n]:
        a=list(map(int,row.split()));typ,nt=a[1:3];dim={2:2,4:3}.get(typ)
        if dim is None: continue
        name=names.get((dim,a[3]),'unlabelled');ind=[byid[v] for v in a[3+nt:]]
        (faces if dim==2 else tets).setdefault(name,[]).append(ind)
    return xyz,faces,tets


def tetra_cut(vertices,tets,axis,value):
    polygons=[];sizes=[]
    for tet in tets:
        p=vertices[tet];d=p[:,axis]-value
        if d.min()>0 or d.max()<0: continue
        hits=[]
        for a,b in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
            if d[a]*d[b]<0: hits.append(p[a]+(p[b]-p[a])*d[a]/(d[a]-d[b]))
            if abs(d[a])<1e-9: hits.append(p[a])
            if abs(d[b])<1e-9: hits.append(p[b])
        if len(hits)<3: continue
        q=np.unique(np.round(np.array(hits),8),axis=0)
        q=np.delete(q,axis,axis=1)
        if len(q)<3: continue
        centre=q.mean(axis=0);angle=np.arctan2(q[:,1]-centre[1],q[:,0]-centre[0]);q=q[np.argsort(angle)]
        polygons.append(q)
        sizes.append(max(np.linalg.norm(p[a]-p[b]) for a in range(4) for b in range(a)))
    return polygons,sizes


def plot_geometry(c,meshdir,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    xyz,faces,tets=parse_mesh(Path(meshdir)/'gaa_nm.msh')
    fig=plt.figure(figsize=(13,5));ax=fig.add_subplot(121,projection='3d')
    # Real mesh contact/interface faces, not artistic blocks.
    for group,color,alpha in [('SiExterior','#5d9fd1',.5),('SpacerExterior','#75bfb3',.35),
                              ('Source','#437ca3',.9),('Drain','#437ca3',.9),('SiOx','#5d9fd1',.6),
                              ('SiSpacer','#5d9fd1',.6),('OxSpacer','#75bfb3',.4),('Gate','#ce8f27',.18)]:
        tri=np.asarray(faces.get(group,[]),int)
        if len(tri): ax.add_collection3d(Poly3DCollection(xyz[tri],facecolors=color,edgecolors='none',alpha=alpha))
    ax.set(xlim=(0,c.length),ylim=(-c.sd_width/2,c.sd_width/2),zlim=(-c.sd_thickness/2,c.sd_thickness/2),
           xlabel='x (nm)',ylabel='y (nm)',zlabel='z (nm)',title='Actual conformal GAA mesh surfaces')
    ax.set_box_aspect((c.length,c.sd_width,c.sd_thickness));ax.view_init(22,-60)
    ax.set_yticks([-c.sd_width/2,0,c.sd_width/2]);ax.set_zticks([-c.sd_thickness/2,0,c.sd_thickness/2])
    ax=fig.add_subplot(122)
    for name,tt in tets.items():
        polygons,_=tetra_cut(xyz,tt,1,.001)
        ax.add_collection(PolyCollection(polygons,facecolors=COLORS[name],edgecolors='none',label=name))
    ax.plot([c.gate_start,c.gate_end],[c.thickness/2+c.tox]*2,color='#806000',lw=3,label='GAA gate boundary')
    ax.plot([c.gate_start,c.gate_end],[-c.thickness/2-c.tox]*2,color='#806000',lw=3)
    ax.autoscale();ax.set_aspect('equal');ax.set(xlabel='x (nm)',ylabel='z (nm)',title='x-z section of original tetrahedra');ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'01_geometry.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(14,4));mesh_collections=[];mesh_sizes=[]
    for name,tt in tets.items():
        polys,sizes=tetra_cut(xyz,tt,0,(c.gate_start+c.gate_end)/2+.001)
        axs[0].add_collection(PolyCollection(polys,facecolors=COLORS[name],edgecolors='none'))
        polys,sizes=tetra_cut(xyz,tt,1,.001)
        coll=PolyCollection(polys,array=np.array(sizes),cmap='viridis',edgecolors='none')
        axs[2].add_collection(coll);mesh_collections.append(coll);mesh_sizes.extend(sizes)
        if name=='Silicon':
            centres=np.array([p.mean(axis=0) for p in polys]);points=np.column_stack([centres[:,0],np.zeros(len(centres)),centres[:,1]])
            coll2=PolyCollection(polys,array=np.log10(np.maximum(c.doping(points),1e10)),cmap='plasma',edgecolors='none')
            axs[1].add_collection(coll2)
    from matplotlib.colors import Normalize
    norm=Normalize(min(mesh_sizes),max(mesh_sizes))
    for collection in mesh_collections: collection.set_norm(norm)
    for ax in axs: ax.autoscale();ax.set_aspect('equal')
    axs[0].set(title='Gate-centre y-z materials',xlabel='y (nm)',ylabel='z (nm)')
    axs[1].set(title='Donors: log10(cm^-3)',xlabel='x (nm)',ylabel='z (nm)');fig.colorbar(coll2,ax=axs[1],shrink=.6)
    axs[2].set(title='Local tetra edge size (nm)',xlabel='x (nm)',ylabel='z (nm)');fig.colorbar(coll,ax=axs[2],shrink=.6)
    fig.tight_layout();fig.savefig(out/'02_material_doping_mesh.png',dpi=160);plt.close(fig)


def rows(path):
    return list(csv.DictReader(path.open())) if path.exists() else []


def plot_run(out):
    out=Path(out);figdir=out/'figures';figdir.mkdir(exist_ok=True)
    iv=rows(out/'idvgs.csv');od=rows(out/'idvds.csv')
    if iv or od:
        fig,axs=plt.subplots(1,2,figsize=(11,4))
        good=[r for r in iv if r['status']=='SCF_CONVERGED_PROTOTYPE']
        if good:
            good.sort(key=lambda r:float(r['vg_v']));v=np.array([float(r['vg_v']) for r in good]);i=np.array([float(r['id_a']) for r in good])
            axs[0].semilogy(v,np.maximum(abs(i),1e-22),'o-',label='SCF converged')
            axs[0].legend()
        failed_label_added=False
        for vg in sorted({float(r['vg_v']) for r in od}):
            group=[r for r in od if float(r['vg_v'])==vg]
            group.sort(key=lambda r:float(r['vd_v']))
            if group:
                values=[float(r['id_a'])*1e6 if r['status']=='SCF_CONVERGED_PROTOTYPE' else np.nan for r in group]
                axs[1].plot([float(r['vd_v']) for r in group],values,'o-',label=f'Vgs={vg:g} V')
                failed=[float(r['vd_v']) for r in group if r['status']!='SCF_CONVERGED_PROTOTYPE']
                if failed:
                    axs[1].scatter(failed,np.full(len(failed),.02),transform=axs[1].get_xaxis_transform(),
                                   marker='x',color='red',label=None if failed_label_added else 'SCF failed (bottom marks)')
                    failed_label_added=True
        if od: axs[1].legend()
        axs[0].set(xlabel='Vgs (V)',ylabel='|Id| (A)',title='Finite-reservoir transfer')
        axs[1].set(xlabel='Vds (V)',ylabel='Id (microampere)',title='Finite-reservoir output')
        if any(r.get('bias_inside_minimum_bandwidth')=='False' for r in od):
            axs[1].set_title('Finite-reservoir output\nCoarse FD: high-Vds fall is not physical NDR',fontsize=10,color='#9b382d')
        for ax in axs: ax.grid(alpha=.25)
        fig.suptitle('Calculated prototype: numerical/continuum convergence not yet certified',fontsize=10)
        fig.tight_layout(rect=(0,0,1,.94));fig.savefig(figdir/'03_iv.png',dpi=160);plt.close(fig)
    hist=rows(out/'scf_history.csv')
    if hist:
        fig,axs=plt.subplots(1,2,figsize=(11,4))
        for key in dict.fromkeys(r['key'] for r in hist):
            rs=[r for r in hist if r['key']==key]
            for ax,col in zip(axs,['phi_residual_v','charge_residual']):
                ax.semilogy([int(r['iteration']) for r in rs],[max(float(r[col]),1e-16) for r in rs],alpha=.75)
        axs[0].set(ylabel='Unmixed potential residual (V)');axs[1].set(ylabel='Unmixed charge residual')
        for ax in axs: ax.set_xlabel('SCF iteration');ax.grid(alpha=.25)
        fig.tight_layout();fig.savefig(figdir/'04_scf.png',dpi=160);plt.close(fig)
    files=sorted((out/'states').glob('vg*.npz')) if (out/'states').exists() else []
    files=[f for f in files if not f.name.endswith(('_spectra.npz','.tmp.npz'))]
    if not files: return
    # Select the last converged snapshot; never present a failed snapshot as accepted.
    valid=[f for f in files if bool(np.load(f)['converged'])]
    if not valid: return
    s=np.load(valid[-1]);mask=s['silicon_mask'];shape=mask.shape
    phi=np.full(shape,np.nan);n=phi.copy();phi[mask]=s['phi_v'];n[mask]=s['electron_density_cm3']
    x,y,z=s['x_nm'],s['y_nm'],s['z_nm'];jy=np.argmin(abs(y));ix=len(x)//2
    fig,axs=plt.subplots(2,2,figsize=(11,7))
    for ax,data,title in [(axs[0,0],phi,'Potential (V)'),(axs[0,1],np.log10(np.maximum(n,1e10)),'log10 n (cm^-3)')]:
        im=ax.pcolormesh(x,z,data[:,jy,:].T,shading='nearest');fig.colorbar(im,ax=ax);ax.set(title=title,xlabel='x (nm)',ylabel='z (nm)')
    for j in range(3): axs[1,0].plot(x,s[f'valley{j}_subbands_ev'][:,0],label=f'valley {j}')
    axs[1,0].set(title='Lowest local subband',xlabel='x (nm)',ylabel='Energy (eV)');axs[1,0].legend()
    u=s['valley0_centre_modes']
    im=axs[1,1].pcolormesh(y,z,abs(u[:,:,0].T)**2,shading='nearest');fig.colorbar(im,ax=axs[1,1])
    axs[1,1].set(title='Gate-centre mode probability',xlabel='y (nm)',ylabel='z (nm)')
    fig.suptitle(f'Vgs={float(s["vg_v"]):g} V, Vds={float(s["vd_v"]):g} V, actual SCF snapshot')
    fig.tight_layout();fig.savefig(figdir/'05_fields_modes.png',dpi=160);plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for j in range(3):
        e=s[f'valley{j}_e'];o=s[f'valley{j}_driven_eigen_occupation'];w=s[f'valley{j}_weights']
        use=(e<float(s['mu_l_ev'])+.5)&(e>float(s['mu_r_ev'])-.3)
        axs[0].scatter(e[use],o[use],s=10,c=w[use,0]-w[use,2],cmap='coolwarm',vmin=-1,vmax=1)
    for mu in [float(s['mu_l_ev']),float(s['mu_r_ev'])]: axs[0].axvline(mu,ls='--',c='gray')
    axs[0].set(xlabel='Finite-system eigenenergy (eV)',ylabel='Driven eigenbasis population',title='Colour: L minus R spatial weight')
    axs[1].plot((x[1:]+x[:-1])/2,s['bond_current_a']*1e6)
    axs[1].set(xlabel='x (nm)',ylabel='Bond electron flux (microampere)',title='Undriven channel continuity')
    fig.tight_layout();fig.savefig(figdir/'06_occupation_current.png',dpi=160);plt.close(fig)
    spectra=valid[-1].with_name(valid[-1].stem+'_spectra.npz')
    if spectra.exists():
        spec=np.load(spectra);fig,ax=plt.subplots(figsize=(7,4))
        ax.plot(spec['energy_ev'],spec['transmission'])
        for mu in [float(s['mu_l_ev']),float(s['mu_r_ev'])]: ax.axvline(mu,ls='--',color='gray')
        ax.set(xlabel='Energy (eV)',ylabel='Spin/valley summed T(E)',title='Post-SCF finite-bath transmission (sampled)')
        ax.grid(alpha=.25);fig.tight_layout();fig.savefig(figdir/'08_transmission.png',dpi=160);plt.close(fig)
    if 'density_from_l_cm3' in s:
        l=np.full(shape,np.nan);r=l.copy();l[mask]=s['density_from_l_cm3'];r[mask]=s['density_from_r_cm3']
        fig,axs=plt.subplots(1,2,figsize=(11,4))
        for ax,d,title in zip(axs,[l,r],['L injected density','R injected density']):
            im=ax.pcolormesh(x,z,np.log10(np.maximum(d[:,jy,:].T,1e10)),shading='nearest');fig.colorbar(im,ax=ax)
            ax.set(title=title+' log10(cm^-3)',xlabel='x (nm)',ylabel='z (nm)')
        fig.tight_layout();fig.savefig(figdir/'07_resolved_density.png',dpi=160);plt.close(fig)
