"""Quantify the retained-subspace mismatch at an abrupt raised-S/D step.

Dx refinement at fixed transverse truncation is not a uniform refinement of
the same continuum domain. A missing narrow-to-wide expansion component is
multiplied by tx~dx^-2. This diagnostic does not correct/fitted-rescale current.
"""
import json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scipy.linalg import orth
from config import Config,VALLEYS,KIN
from modespace.hamiltonian import Grid,ModeBuilder


def run(out):
    rows=[]
    for channel,sd in [(2,4),(3,6),(5,10),(7,14),(9,18)]:
        c=Config(dx=.5,modes=channel,sd_modes=sd);g=Grid.build(c)
        for v in VALLEYS:
            s=ModeBuilder(c,g,v).build(np.zeros(len(g.xyz)))
            i=np.where(g.x<c.ls)[0][-1]
            wide=s.bases[i].reshape(-1,s.bases[i].shape[-1]);narrow=s.bases[i+1].reshape(-1,s.bases[i+1].shape[-1])
            overlap=wide.T@narrow;singular=np.linalg.svd(overlap,compute_uv=False)
            deficit=max(0.,float(np.linalg.eigvalsh(np.eye(narrow.shape[1])-overlap.T@overlap).max()))
            enriched=orth(np.column_stack([wide,narrow]));exact=enriched.T@narrow
            row={'valley':v.name,'channel_modes':channel,'sd_modes':sd,
                 'overlap_singular_values':singular.tolist(),'maximum_missing_norm_squared':deficit,
                 'unresolved_interface_scale_ev':{str(dx):KIN/(v.mx*dx**2)*deficit for dx in [.5,.25,.125]},
                 'span_enrichment_rank':enriched.shape[1],
                 'span_enrichment_missing_norm_squared':float(np.linalg.norm(np.eye(narrow.shape[1])-exact.T@exact))}
            rows.append(row)
    result={'results':rows,'scope':'zero-potential subspace diagnostic, not SCF or an implemented production enrichment',
            'implication':'Couple longitudinal and transverse/mode refinement; severe truncated overlap can create an artificial abrupt-junction barrier',
            'next_implementation':'A junction-compatible enriched wide basis can contain the zero-extended narrow subspace; validate full projected kinetic/potential and contacts before using it'}
    Path(out).write_text(json.dumps(result,indent=2));return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--out',default='results/junction_overlap.json');a=p.parse_args();run(a.out)
