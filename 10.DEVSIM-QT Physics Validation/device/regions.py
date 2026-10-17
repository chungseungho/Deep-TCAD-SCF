"""Five geometric regions; reservoir placement is diagnosed, never fitted."""
import numpy as np

NAMES=('left_reservoir','left_extension','channel','right_extension','right_reservoir')


def five_regions(coordinate,c,shift=0.):
    x=np.asarray(coordinate)
    left=c.reservoir_left_nm+shift;right=c.length-c.reservoir_right_nm-shift
    if not 0<left<=c.ls or not c.drain_start<=right<c.length:
        raise ValueError('Shifted reservoir lies outside electrode')
    return np.select([x<left,x<c.gate_start,x<c.gate_end,x<right],[0,1,2,3],default=4)


def three_regions(coordinate,c,shift=0.):
    five=five_regions(coordinate,c,shift)
    return np.where(five==0,0,np.where(five==4,2,1))


def mode_regions(system,c,shift=0.):
    labels=three_regions(system.grid.x,c,shift)
    return np.repeat(labels,np.diff(system.offsets))


def diagnose(grid,density,potential,c):
    labels=five_regions(grid.xyz[:,0],c);dv=grid.cell_volume*1e-21
    donor=c.doping(grid.xyz);report={}
    for i,name in enumerate(NAMES):
        sel=labels==i;ne=float(density[sel].sum()*dv);nd=float(donor[sel].sum()*dv)
        ph=potential[sel]
        xs=grid.xyz[sel,0];unique=np.unique(xs)
        means=np.array([ph[xs==x].mean() for x in unique])
        gradient=float(max(abs(np.diff(means)/np.diff(unique)))) if len(unique)>1 else 0.
        report[name]={'electrons':ne,'donor_electrons':nd,'net_electrons':ne-nd,
            'neutrality_relative':(ne-nd)/nd if nd else None,
            'potential_spread_v':float(np.ptp(ph)),'max_mean_gradient_v_nm':gradient}
    report['reservoir_screening_diagnostic_pass']=bool(all(
        abs(report[n]['neutrality_relative'])<.05 and report[n]['max_mean_gradient_v_nm']<.001
        for n in (NAMES[0],NAMES[-1])))
    report['screening_thresholds']={'relative_neutrality':.05,'mean_gradient_v_nm':.001,
        'meaning':'diagnostic screening targets, not universal semiconductor validity bounds'}
    return report
