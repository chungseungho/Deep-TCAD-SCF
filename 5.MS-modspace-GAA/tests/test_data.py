import sys,json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pytest
from data import export_dataset
from surrogate import train,predict


def test_prototype_export_is_explicit_and_split_is_by_bias(tmp_path):
    run=tmp_path/'run';(run/'states').mkdir(parents=True)
    config={'ls':2,'lg':4,'spacer':1}
    (run/'settings.json').write_text(json.dumps({'config':config,'fingerprint':'test-geometry'}))
    xyz=np.array([[.5,0,0],[2.5,0,0],[3.5,0,0],[6.5,0,0],[8.5,0,0]])
    for i in range(10):
        np.savez_compressed(run/'states'/f'vg{i:02d}_vd.npz',converged=True,quantitative_validation_pass=False,
             xyz_nm=xyz,initial_neutral_poisson_phi_v=np.full(5,i*.01),phi_v=np.full(5,99.),
             electron_density_cm3=np.full(5,1e17*(1+i*.2)),vg_v=i*.1,vd_v=.1)
    # An interrupted temporary snapshot must not be read as a real result.
    (run/'states'/'vg99_vd.tmp.npz').write_text('interrupted')
    with pytest.raises(ValueError,match='No eligible'): export_dataset(run,tmp_path/'blocked.npz')
    file=tmp_path/'dataset.npz';report=export_dataset(run,file,True);d=np.load(file)
    assert report['snapshots']==10 and report['heldout_bias_points']==2
    assert d['input_poisson_phi_v'].max()<1  # cannot leak the final phi=99 target
    assert d['test_mask'].shape==(10,)  # entire bias points, no random spatial split
    model=tmp_path/'baseline.npz';r=train(file,model)
    assert r['training_physics_quantitatively_validated'] is False
    pred=predict(model,.4,.1,'test-geometry',xyz[d['channel_mask']])
    assert pred.min()>0
    with pytest.raises(ValueError,match='mismatch'): predict(model,.4,.1,'another-grid',xyz[d['channel_mask']])
    with pytest.raises(ValueError,match='extrapolation'): predict(model,2.,.1,'test-geometry',xyz[d['channel_mask']])
