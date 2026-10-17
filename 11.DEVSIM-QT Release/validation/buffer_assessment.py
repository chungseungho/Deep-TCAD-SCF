import json
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1]
a=np.load(ROOT/'results/raised_eq_final/states/vg0.400000_vd0.000000.npz');b=np.load(ROOT/'results/buffer_eq/states/vg0.400000_vd0.000000.npz')
xyz=b['xyz_nm'].copy();xyz[:,0]-=6;distance,index=cKDTree(xyz).query(a['xyz_nm'])
if max(distance)>1e-9:raise ValueError('Original target domain not preserved')
x=a['xyz_nm'][:,0];na=a['electron_density_cm3'];nb=b['electron_density_cm3'][index]
rows={}
for name,mask in [('channel',(x>=15)&(x<=39)),('junction_extension',((x>=12)&(x<15))|((x>39)&(x<=42))),('original_reservoir',(x<12)|(x>42))]:
 rows[name]={'density_l2':float(np.linalg.norm((na-nb)[mask])/np.linalg.norm(nb[mask])),'phi_max_difference_v':float(max(abs(a['phi_v'][mask]-b['phi_v'][index][mask])))}
report={'category':'B: computational lead extension with original gate/channel/junction coordinates restored','outer_extension_each_side_nm':6,'original_target_domain_matched':True,'potential_reference_and_neutral_mu_unchanged':True,'bias_v':0.,'regions':rows,'status':'FAIL' if any(r['density_l2']>.01 for r in rows.values()) else 'INCONCLUSIVE','reason':'One extension at equilibrium does not establish a plateau or biased transport convergence; remeshing is a secondary numerical contributor.'}
(ROOT/'results/assessment/buffer_comparison.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
