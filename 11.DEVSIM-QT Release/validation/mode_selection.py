"""Conservative starting basis plus mandatory observable-driven refinement.

An energy/occupation criterion is a proposal, never a convergence certificate.
Degenerate multiplets are retained. Junction evanescent guard is explicit.
"""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from config import Config,KIN,KB_EV
from modespace.hamiltonian import Grid,transverse_matrix
from device.definition import DeviceDefinition
from transport.contact_model import solve_reference

def select(config,out,evanescent_guard_ev=.1):
    c=Config(**json.loads(Path(config).read_text()));g=Grid.build(c);valleys=DeviceDefinition(transport_axis=c.transport_axis).frame.valleys();r=solve_reference(c,g,valleys)
    report={'scope':'starting proposal only; run refine at off/threshold/on and low/high drain','contact_full_basis':len(r.phi),'evanescent_guard_ev':evanescent_guard_ev,'valleys':[],'numerical_convergence_verified':False};nm=1;ns=1
    for v in valleys:
        sd=r.spectra[v.name][0];i=np.flatnonzero((g.x>=c.gate_start)&(g.x<=c.gate_end))[0]
        ch=eigh(transverse_matrix(g.mask[i],KIN/(v.my*c.dy**2),KIN/(v.mz*c.dz**2)).toarray(),eigvals_only=True)
        # Scan configured gate window: a positive gate can lower channel bands.
        channel_ceiling=r.mu+1.2+12*KB_EV*c.temperature+evanescent_guard_ev
        contact_ceiling=r.mu+12*KB_EV*c.temperature+evanescent_guard_ev
        def count(e,cut):
            n=max(1,int(np.sum(e<=cut)))
            while n<len(e) and abs(e[n]-e[n-1])<1e-6:n+=1
            return min(n,len(e))
        a,b=count(ch,channel_ceiling),count(sd,contact_ceiling);nm=max(nm,a);ns=max(ns,b)
        report['valleys'].append({'valley':v.name,'channel_modes':a,'contact_modes':b,'contact_cutoff_ev':contact_ceiling,'channel_cutoff_ev':channel_ceiling})
    values=c.to_dict();values.update(modes=nm,sd_modes=ns);Path(out).write_text(json.dumps(values,indent=2))
    report.update(requested_modes=nm,requested_sd_modes=ns,warning='A full gate window can require full channel basis; nested enrichment adds further contact states. Resource and fully SCF convergence checks remain mandatory.')
    Path(out).with_suffix('.selection.json').write_text(json.dumps(report,indent=2));return report
