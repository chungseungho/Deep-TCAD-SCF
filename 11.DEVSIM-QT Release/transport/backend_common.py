"""Shared grid/Hamiltonian setup for built-in transport implementations."""
import numpy as np
from modespace.hamiltonian import density_from_matrix
from transport.results import TransportResult


def initialize(context):
    c=context['config']; builders=context['builders']; grid=context['grid']; phi=context['potential']
    systems=[b.build(phi) for b in builders]
    return {'c':c,'grid':grid,'bias':context.get('bias',0.),'systems':systems,
        'deg':[2*s.valley.multiplicity for s in systems],
        'density':np.zeros(len(grid.xyz)),'left':np.zeros(len(grid.xyz)),
        'right':np.zeros(len(grid.xyz)),'unresolved':np.zeros(len(grid.xyz)),
        'states':[],'diagnostics':[],'matrices':[],
        'bonds':np.zeros(len(grid.x)-1),'tracking':[],
        'mu':float(context['fixed_mu']),'constraint':{},'detail':context.get('detail',False)}


def finish(d):
    c=d['c'];grid=d['grid'];bonds=d['bonds']
    current=float('nan') if d.get('current_is_undefined',False) else float(
        bonds[(grid.x[:-1]>=c.ls)&(grid.x[1:]<=c.drain_start)].mean())
    optional={'config':d['c'],'systems':d['systems'],'states':d['states'],'matrices':d['matrices'],
        'tracking':np.asarray(d['tracking']),'constraint':d['constraint'],
        'density_l':d['left'],'density_r':d['right'],
        'density_unresolved':d['unresolved'],'valleys':d['diagnostics']}
    return TransportResult(d['density'],current,bonds,d['mu'],d['mu']-d['bias'],
        diagnostics={'valleys':d['diagnostics']},optional=optional)

