from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[2];out=Path(__file__).parent
s=np.load(root/'results/final_raised_equilibrium/states/vg0.400000_vd0.000000.npz');xyz=s['xyz_nm'];axes=[np.unique(xyz[:,i]) for i in range(3)]
indices=tuple(np.searchsorted(a,xyz[:,i]) for i,a in enumerate(axes));shape=tuple(map(len,axes));phi=np.full(shape,np.nan);density=phi.copy();phi[indices]=s['phi_v'];density[indices]=s['electron_density_cm3']
y=int(np.argmin(abs(axes[1])));x,z=axes[0],axes[2];extent=(x[0]-.25,x[-1]+.25,z[0]-.5,z[-1]+.5)
fig,ax=plt.subplots(2,1,figsize=(10,5),sharex=True)
a=ax[0].imshow(phi[:,y,:].T,origin='lower',extent=extent,aspect='auto',cmap='coolwarm');fig.colorbar(a,ax=ax[0],label='Cell-average potential (V)')
b=ax[1].imshow(np.log10(np.maximum(density[:,y,:].T,1.)),origin='lower',extent=extent,aspect='auto',cmap='viridis');fig.colorbar(b,ax=ax[1],label='log10 density (cm^-3)')
for a in ax:a.set_ylabel('z (nm)')
ax[1].set_xlabel('x (nm)');fig.suptitle('Actual Raised S/D equilibrium SCF, central-y slice\nNumerical/physical device convergence not certified');fig.tight_layout();fig.savefig(out/'raised_equilibrium_slices.png',dpi=180)
