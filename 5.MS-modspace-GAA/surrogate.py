#!/usr/bin/env python3
"""Small RBF-ridge log-density baseline. NOT the lab's DeepSCF neural network.

Prediction is a warm start and is always checked by the full reference SCF.
This module supplies data plumbing and a reproducible whole-bias holdout, not
evidence that a trained neural model can replace the quantum solver.
"""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from scipy.linalg import solve


def train(dataset,destination):
    d=np.load(dataset);test=d['test_mask'];b=d['bias_v']
    if len(b)<8 or test.sum()<1: raise ValueError('Need >=8 complete bias points with an independent holdout')
    x=b[~test];y=np.log10(np.maximum(d['target_channel_density_cm3'][~test],1e10))
    scale=np.maximum(np.ptp(x,axis=0),.05);centre=x.mean(axis=0);xn=(x-centre)/scale
    distances=cdist(xn,xn,'sqeuclidean');nonzero=distances[distances>0]
    length=float(np.sqrt(np.median(nonzero))) if len(nonzero) else 1.
    K=np.exp(-distances/(2*length**2));mean=y.mean(axis=0)
    U,s,Vt=np.linalg.svd(y-mean,full_matrices=False);rank=min(12,len(s))
    components=Vt[:rank];codes=(y-mean)@components.T
    weights=solve(K+1e-7*np.eye(len(K)),codes,assume_a='pos')
    xt=(b[test]-centre)/scale;kt=np.exp(-cdist(xt,xn,'sqeuclidean')/(2*length**2))
    prediction=10**np.clip(mean+(kt@weights)@components,10,23)
    actual=d['target_channel_density_cm3'][test]
    report={'model':'RBF-ridge reduced log-density baseline; not DeepSCF',
            'heldout_entire_bias_points':b[test].tolist(),
            'heldout_log10_rmse':float(np.sqrt(np.mean((np.log10(np.maximum(actual,1e10))-np.log10(prediction))**2))),
            'heldout_relative_l2':float(np.linalg.norm(prediction-actual)/np.linalg.norm(actual)),
            'training_physics_quantitatively_validated':bool(np.all(d['quantitative_validation_pass'])),
            'inference_role':'warm start; reference SCF mandatory'}
    destination=Path(destination);destination.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(destination,train_x=xn,scale=scale,centre=centre,length=length,weights=weights,
                        mean=mean,components=components,bias_min=x.min(axis=0),bias_max=x.max(axis=0),
                        fingerprint=d['fingerprint'],channel_xyz_nm=d['channel_xyz_nm'])
    destination.with_suffix('.json').write_text(json.dumps(report,indent=2));return report


def predict(model,vg,vd,fingerprint,channel_xyz):
    m=np.load(model);bias=np.array([vg,vd])
    if str(m['fingerprint'])!=fingerprint or not np.array_equal(m['channel_xyz_nm'],channel_xyz):
        raise ValueError('Surrogate geometry/physics/grid mismatch')
    if np.any(bias<m['bias_min']-1e-9) or np.any(bias>m['bias_max']+1e-9):
        raise ValueError('Surrogate refuses extrapolation beyond training bias box')
    x=(bias-m['centre'])/m['scale'];K=np.exp(-np.sum((m['train_x']-x)**2,axis=1)/(2*float(m['length'])**2))
    return 10**np.clip(m['mean']+(K@m['weights'])@m['components'],10,23)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('dataset');p.add_argument('destination');a=p.parse_args()
    print(json.dumps(train(a.dataset,a.destination),indent=2))
