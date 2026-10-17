"""Transparent timing wrapper around the real driver; no numerical substitutions."""
import sys,time,json,functools
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);args=sys.argv[2:]
if args and args[0]=='--':args=args[1:]
counts={};last=time.monotonic();start=last

def wrap(module,name,label,announce=False):
    original=getattr(module,name)
    @functools.wraps(original)
    def measured(*a,**kw):
        global last
        tic=time.perf_counter()
        if announce:print('PROFILE START',label,flush=True)
        try:return original(*a,**kw)
        finally:
            row=counts.setdefault(label,{'calls':0,'inclusive_seconds':0.,'energy_nodes':0});row['calls']+=1;row['inclusive_seconds']+=time.perf_counter()-tic
            if name in ('rgf','lead_matrices'):row['energy_nodes']+=len(a[0])
            if announce or time.monotonic()-last>30:
                (out/'stage_progress.json').write_text(json.dumps({'elapsed_s':time.monotonic()-start,'stages':counts},indent=2));last=time.monotonic()
                print('PROFILE',label,json.dumps(row),flush=True)
    setattr(module,name,measured)

import transport.negf_reference as nr
import transport.resonances as rs
import transport.quadrature as quad
import transport.bound_states as bs
import transport.contact_model as cm
from modespace.hamiltonian import ModeBuilder
from modespace.junction_basis import JunctionCompatibleBuilder
from poisson.device import Poisson
from coupling.mapping import Mapping
for module,name,label,announce in [(nr,'rgf','RGF contour/injectivity',False),(nr,'lead_matrices','contact self-energy',False),(rs,'resonance_breaks','banded resonance search',True),(quad,'integrate_vector','real-energy integration',True),(bs,'audit','bound audit',False),(cm,'solve_reference','lead transverse SP',True),(ModeBuilder,'build','transverse basis/base blocks',False),(JunctionCompatibleBuilder,'build','nested blocks',False),(Poisson,'solve','Poisson solve',False)]:wrap(module,name,label,announce)
import run
sys.argv=['run.py',*args]
try:code=run.main()
finally:(out/'stage_totals.json').write_text(json.dumps({'elapsed_s':time.monotonic()-start,'stages':counts,'timings_are_inclusive_do_not_sum_nested_categories':True},indent=2))
raise SystemExit(code)
