#!/usr/bin/env python3
"""Generate reviewable figures/status from actual saved values only."""
import csv,json,re,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from validation.convergence import summarize
from provenance import source_digest


def rows(path):return list(csv.DictReader(path.open())) if path.exists() else []
def read(path):return json.loads(path.read_text()) if path.exists() else {}
def accepted(r):return r.get('status')=='SCF_CONVERGED_PROTOTYPE'
CAMPAIGNS=['convergence','extended_convergence','additional_convergence']


def parameter_label(group,c):
    if group=='gamma':return f"{c['gamma']*1000:g}"
    if group=='dx':return f"{c['dx']:g}"
    if group=='transverse':return f"{c['dy']:g}/{c['dz']:g}"
    if group=='modes':return f"{c['modes']}/{c['sd_modes']}"
    if group=='length':return f"{c['ls']:g}/{c['ld']:g}"
    if group=='mesh':return f"{c['mesh_h']:g}/{c['channel_h']:g}/{c['junction_h']:g}"
    if group=='basis':return f"{c['basis']} {c['modes']}/{c['sd_modes']}"
    if group=='scf':return f"{c['phi_tol']:.2g}/{c['charge_tol']:.2g}/{c['current_tol']:.2g}"
    return group


PARAMETER_UNITS={'gamma':'Gamma (meV)','dx':'dx (nm)','transverse':'dy/dz (nm)',
    'modes':'Channel / S-D modes','length':'Ls/Ld (nm)','mesh':'h/channel/junction (nm)',
    'basis':'Basis and channel/S-D modes','scf':'Phi (V) / charge / current limits'}


def figures(root):
    out=root/'results';dest=out/'figures';dest.mkdir(exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120})
    transfer=rows(out/'iv/idvgs.csv');output=rows(out/'output/idvds.csv')
    if transfer:
        transfer.sort(key=lambda r:float(r['vg_v']))
        v=np.array([float(r['vg_v']) for r in transfer]);i=np.array([float(r['id_a']) if accepted(r) else np.nan for r in transfer])
        gm=np.full(len(v),np.nan)
        for j in range(1,len(v)-1):
            if np.isfinite(i[j-1:j+2]).all():gm[j]=(i[j+1]-i[j-1])/(v[j+1]-v[j-1])
        fig,ax=plt.subplots(1,3,figsize=(13,3.8))
        for a,d,title,unit in zip(ax,[i*1e6,abs(i),gm*1e6],['Transfer','Transfer, log','Transconductance'],['Id (uA)','|Id| (A)','gm (uA/V)']):
            a.plot(v,d,'o-',ms=3,color='#236e90');a.set(xlabel='Vgs (V)',ylabel=unit,title=title);a.grid(alpha=.2)
        ax[1].set_yscale('log')
        fig.suptitle('Actual SCF, Vds=0.1 V | provisional numerical/contact configuration',fontsize=11)
        fig.tight_layout();fig.savefig(dest/'transfer_gm.png',dpi=170);plt.close(fig)
    if output:
        fig,ax=plt.subplots(figsize=(6.8,4.2))
        for vg in sorted({float(r['vg_v']) for r in output}):
            group=sorted([r for r in output if float(r['vg_v'])==vg],key=lambda r:float(r['vd_v']))
            ax.plot([float(r['vd_v']) for r in group],[float(r['id_a'])*1e6 if accepted(r) else np.nan for r in group],'o-',ms=3,label=f'Vgs={vg:g} V')
        ax.set(xlabel='Vds (V)',ylabel='Id (uA)',title='Actual SCF output; reference convergence not certified')
        ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(dest/'output.png',dpi=170);plt.close(fig)
    history_cases=[('device_lowbias','point','vg0.800000_vd0.100000','Finite, modes 3/6'),
                   ('output','idvds','vg1.200000_vd0.050000','Output difficult bias, all attempts'),
                   ('nested_negf_refined','point','vg0.800000_vd0.100000','Biased NEGF, incomplete')]
    fig,axes=plt.subplots(1,2,figsize=(11.5,4))
    for directory,phase,key,label in history_cases:
        history=[r for r in rows(out/directory/'scf_history.csv') if r['key']==key]
        if not history:continue
        final=next((r for r in rows(out/directory/(phase+'.csv')) if r['key']==key and accepted(r)),None)
        for ax,field in zip(axes,['phi_residual_v','charge_residual']):
            line=ax.semilogy(np.arange(1,len(history)+1),[max(float(r[field]),1e-16) for r in history],
                             label=label,linestyle='-' if final else '--')[0]
            if final:ax.plot(len(history)+1,max(float(final[field]),1e-16),'*',color=line.get_color(),ms=10)
    for ax,tolerance,ylabel in zip(axes,[1e-4,2e-3],['Potential residual (V)','Relative charge residual']):
        ax.axhline(tolerance,color='black',ls=':',lw=1,label='Acceptance limit')
        ax.set(xlabel='Recorded update index',ylabel=ylabel);ax.grid(alpha=.2)
    axes[1].legend(fontsize=8)
    fig.suptitle('SCF history, restarts retained; stars are residuals on the exact saved state',fontsize=11)
    fig.tight_layout();fig.savefig(dest/'scf_history.png',dpi=170);plt.close(fig)
    files=sorted((out/'iv/states').glob('vg*.npz'))
    valid=[p for p in files if bool(np.load(p)['converged'])]
    if valid:
        s=np.load(valid[-1]);xyz=s['xyz_nm'];mask=s['silicon_mask'];shape=mask.shape
        field=np.full(shape,np.nan);density=field.copy();field[mask]=s['phi_v'];density[mask]=s['electron_density_cm3']
        x,y,z=s['x_nm'],s['y_nm'],s['z_nm'];iy=np.argmin(abs(y));iz=np.argmin(abs(z))
        # Show exterior Si cells. A stride through flattened regular cells
        # aliases the z coordinate and can misleadingly look like stacked sheets.
        padded=np.pad(mask,1);interior=mask.copy()
        for axis in range(3):
            for shift in [-1,1]:
                slices=[slice(1,-1)]*3
                slices[axis]=slice(0,-2) if shift<0 else slice(2,None)
                interior &= padded[tuple(slices)]
        keep=np.flatnonzero((mask & ~interior)[mask])
        if len(keep)>18000:keep=np.sort(np.random.default_rng(6).choice(keep,18000,replace=False))
        fig=plt.figure(figsize=(13.5,5.2))
        for j,(values,label,cmap) in enumerate([(np.log10(np.maximum(s['electron_density_cm3'],1e9)),'log10 n (cm^-3)','viridis'),(s['phi_v'],'Potential (V)','coolwarm')]):
            ax=fig.add_subplot(1,2,j+1,projection='3d')
            image=ax.scatter(*xyz[keep].T,c=values[keep],s=5,marker='s',depthshade=False,cmap=cmap,rasterized=True)
            ax.set(xlabel='x (nm)',ylabel='y (nm)',zlabel='z',title=label);ax.set_box_aspect((54,20,15));ax.view_init(22,-65)
            from matplotlib.ticker import MaxNLocator
            ax.yaxis.set_major_locator(MaxNLocator(4));ax.zaxis.set_major_locator(MaxNLocator(4));ax.tick_params(labelsize=8,pad=1)
            ax.zaxis.labelpad=0
            fig.colorbar(image,ax=ax,shrink=.65,pad=.13)
        fig.suptitle(f'Si exterior cells, coordinates in nm: Vgs={float(s["vg_v"]):g} V, Vds={float(s["vd_v"]):g} V')
        fig.tight_layout();fig.savefig(dest/'density_potential_3d.png',dpi=170);plt.close(fig)
        fig,ax=plt.subplots(1,3,figsize=(13,3.8))
        ax[0].plot(x,field[:,iy,iz],label='centre-line phi');ax[0].set(xlabel='x (nm)',ylabel='Potential (V)',title='Source / channel / drain')
        ax[1].plot(x,-field[:,iy,iz],label='Ec, Si reference=0')
        for v in range(3):ax[1].plot(x,s[f'valley{v}_subbands_ev'][:,0],label=f'lowest subband {v}')
        ax[1].set(xlabel='x (nm)',ylabel='Energy (eV)',title='Conduction edge and confinement');ax[1].legend(fontsize=7)
        c=read(out/'iv/settings.json')['config'];line=np.nansum(density,axis=(1,2))*c['dy']*c['dz']*1e-21
        ax[2].plot(x,line);ax[2].set(xlabel='x (nm)',ylabel='Electrons / nm',title='Cross-section integrated density')
        for a in ax:a.grid(alpha=.2)
        fig.tight_layout();fig.savefig(dest/'device_profiles.png',dpi=170);plt.close(fig)
    spectra=out/'post_scf_current/transmission.npz'
    if spectra.exists():
        s=np.load(spectra);fig,ax=plt.subplots(figsize=(6.8,4))
        ax.plot(s['energy_ev'],s['spin_valley_summed_transmission']);ax.set(xlabel='Energy (eV)',ylabel='Spin/valley summed T(E)',title='Actual GAA finite-bath independent post-SCF spectrum')
        ax.grid(alpha=.2);fig.tight_layout();fig.savefig(dest/'transmission.png',dpi=170);plt.close(fig)
    spectra=out/'frozen_contact_comparison/transmission.npz'
    if spectra.exists():
        s=np.load(spectra);fig,ax=plt.subplots(figsize=(7.4,4))
        ax.plot(s['energy_ev'],s['finite_bath_spin_valley_summed_transmission'],label='Finite distributed bath')
        ax.plot(s['energy_ev'],s['semiinfinite_spin_valley_summed_transmission'],label='Semi-infinite wire leads',alpha=.8)
        ax.set(xlabel='Energy (eV)',ylabel='Spin/valley summed T(E)',title='Same frozen GAA field; contact comparison only')
        ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(dest/'frozen_contact_comparison.png',dpi=170);plt.close(fig)
    allcases=[];configurations={}
    for directory in CAMPAIGNS:
        p=out/directory
        if p.exists():
            allcases+=summarize(p)['cases']
            for file in p.glob('*.json'):
                c=read(file)
                if isinstance(c,dict) and 'backend' in c and 'dx' in c:configurations[file.stem]=c
    groups=sorted({name.rsplit('_',1)[0] for name in configurations})
    if groups:
        fig,axes=plt.subplots((len(groups)+2)//3,3,figsize=(12,3.5*((len(groups)+2)//3)),squeeze=False)
        for ax,group in zip(axes.ravel(),groups):
            cases=sorted([name for name in configurations if name.rsplit('_',1)[0]==group])
            saved={r['case']:r for r in allcases}
            ax.plot(range(len(cases)),[float(saved[name]['id_a'])*1e6 if name in saved and accepted(saved[name]) else np.nan for name in cases],'o-')
            ax.set_xticks(range(len(cases)),[parameter_label(group,configurations[name]) for name in cases],rotation=25 if group in ['mesh','basis','scf'] else 0)
            ax.tick_params(axis='x',labelsize=8)
            missing=[i for i,name in enumerate(cases) if name not in saved or not accepted(saved[name])]
            if missing:
                ax.scatter(missing,np.full(len(missing),.035),transform=ax.get_xaxis_transform(),color='#a23b32',marker='x',label='No accepted result')
                ax.legend(fontsize=7)
            ax.set_xlim(-.15,len(cases)-.85)
            ax.set(title=group,xlabel=PARAMETER_UNITS.get(group,group),ylabel='Id (uA)');ax.grid(alpha=.2)
        for ax in axes.ravel()[len(groups):]:ax.axis('off')
        fig.suptitle('Fully coupled 3D SCF sensitivity at Vgs=0.8 V, Vds=0.1 V; not a global certificate')
        fig.tight_layout();fig.savefig(dest/'numerical_convergence.png',dpi=170);plt.close(fig)


def write_reports(root):
    result=root/'results';checks={}
    for label,file in [('transfer','iv/idvgs.csv'),('output','output/idvds.csv')]:
        rs=rows(result/file);checks[label]={'requested':25 if label=='transfer' else 75,'computed':len(rs),'accepted':sum(accepted(r) for r in rs),
             'failed_biases':[[float(r['vg_v']),float(r['vd_v']),r['status']] for r in rs if not accepted(r)]}
    physics=read(result/'benchmarks/summary.json');current=read(result/'post_scf_current/summary.json');old=read(result/'step05_interface_audit.json')
    test=(result/'tests.txt').read_text() if (result/'tests.txt').exists() else ''
    excluded={'development_point','negf_pilot','spectral_pilot'}
    mappings=[read(path) for path in result.rglob('mapping.json') if not set(path.relative_to(result).parts)&excluded]
    states=[read(path) for path in result.rglob('vg*.json') if path.parent.name=='states' and not set(path.relative_to(result).parts)&excluded]
    good=[s for s in states if accepted(s)]
    campaigns={name:summarize(result/name) for name in CAMPAIGNS if (result/name).exists()}
    nested={name:rows(result/name/'point.csv') for name in ['nested_spectral_dx0p5','nested_spectral_dx0p25']}
    nested_change=None
    if all(nested.values()) and all(accepted(value[0]) for value in nested.values()):
        a,b=[float(value[0]['id_a']) for value in nested.values()]
        nested_change=abs(a/b-1)
    metrics=read(result/'iv/transfer_metrics.json')
    frozen=read(result/'frozen_contact_comparison/summary.json')
    cross_sweep=read(result/'cross_sweep_consistency.json')
    biased_references={}
    for name in ['negf_lowbias','nested_negf','nested_negf_refined']:
        table=rows(result/name/'point.csv');history=rows(result/name/'scf_history.csv')
        biased_references[name]={'completed_bias_rows':table,
            'execution_status':table[-1]['status'] if table else 'NO_ACCEPTED_FINAL_STATE',
            'last_saved_scf_iteration':history[-1] if history else None}
    level3=bool(good and all(float(s.get('gauss_relative_error',1))<1e-7 for s in good))
    tests_pass=bool(re.search(r'\b[0-9]+ passed',test)) and not bool(re.search(r'\b[1-9][0-9]* (?:failed|errors?)\b',test))
    summary={'sweeps':checks,'levels':{'1_code':True,'2_tests':tests_pass,'3_device_scf':level3,'4_numerical_reference':False,'5_ai_teacher':False},
             'delivered_core_source_sha256':source_digest(),
             'saved_state_core_sources':sorted({s['source_sha256'] for s in states if 'source_sha256' in s}),
             'source_snapshots':read(result/'source_snapshots/index.json'),
             'all_results_quantitatively_validated':False,'simulation_model':'electron-only ballistic hard-wall Hartree EMA, finite wide-band baths',
             'semiinfinite_equilibrium_scf':read(result/'negf_equilibrium/states/vg0.800000_vd0.000000.json'),
             'semiinfinite_biased_status':'unresolved narrow-resonance integration; ERROR retained',
             'semiinfinite_biased_attempts':biased_references,
             'maximum_saved_gauss_relative_error':max((s.get('gauss_relative_error',0) for s in good),default=None),
             'maximum_saved_current_conservation_absolute_a':max((s.get('current_conservation_abs_a',0) for s in good),default=None),
             'maximum_saved_current_conservation_relative':max((s.get('current_conservation_relative',0) for s in good),default=None),
             'maximum_mapping_constant_error':max((s.get('constant_density_relative_error',0) for s in mappings),default=None),
             'maximum_mapping_work_defect':max((s.get('mapping_work_relative_defect',0) for s in good),default=None),
             'current_comparison':current,'transfer_metrics_provisional':metrics,
             'frozen_contact_current_only':frozen,
             'cross_sweep_consistency':cross_sweep,
             'field_export_integrity':read(result/'field_export_audit.json'),
             'resume_export_regression':read(result/'resume_export_regression/resume_artifact_check.json'),
             'numerical_campaigns':campaigns,'nested_trace_comparison':nested,
             'nested_two_grid_current_relative_change':nested_change,
             'nested_comparison_plateau_certified':False}
    (result/'step06_summary.json').write_text(json.dumps(summary,indent=2))
    t,o=checks['transfer'],checks['output']
    status=f'''# STEP06 — actual implementation and validation status

**Level 3: 실제 3D GAA fully SCF를 검증했다. Level 4/5 정량 reference 또는 AI teacher는 아직 아니다.**

| 판정 | 실제 근거 |
|---|---|
| Level 1 — Code implemented | 독립 Gmsh–DEVSIM–EMA solver, finite spectral / NEGF comparison |
| Level 2 — Numerical tests passed | {test.strip().splitlines()[-1] if test.strip() else 'pending'} |
| Level 3 — Device SCF converged | 실제 GAA에서 unmixed phi/n/current, Gauss 검사 |
| Level 4 — Numerical reference | 미달: contact/grid/mode/basis 등 전 영역 plateau 미확립 |
| Level 5 — AI teacher ready | 미달: 모든 state의 quantitative_validation_pass=False |

| 실제 sweep | 계산 / 목표 | accepted SCF | 실패 |
|---|---:|---:|---:|
| Id–Vgs, Vds=.1 V | {t['computed']}/25 | {t['accepted']} | {t['computed']-t['accepted']} |
| Id–Vds, Vgs=.4,.8,1.2 V | {o['computed']}/75 | {o['accepted']} | {o['computed']-o['accepted']} |

현재 sweep 설정은 dx=.25 nm, dy=dz=1 nm, modes=2/4, gamma=10 meV, coarse Poisson mesh다.
이는 **실제로 계산된 provisional I–V**이며 mode/contact/continuum convergence를 인증한 소자 사양이 아니다.
별도 campaign에서 dx, gamma, mode, length, basis, Poisson/transverse mesh, SCF tolerance를 실제 coupled SCF로 비교했다. 현재 SCF tolerance 외의 축은 지정한 plateau 기준을 통과하지 못했다. Fine dx의 세 번째 coupled 계산은 미완료다.

## 검출하고 수정한 오류

- STEP05 원본을 saved bias Vgs=.4/Vds=.1 V에서 다시 읽어 실행했을 때 Si/oxide 공유 전위 점프가 {old.get('silicon_oxide',{}).get('maximum_potential_jump_v',float('nan')):.6g} V, oxide/spacer는 {old.get('oxide_spacer',{}).get('maximum_potential_jump_v',float('nan')):.6g} V였다. 삼중 접합의 겹치는 interface 조건을 공유 P1 자유도로 교체했다. 물질별 epsilon은 유지한다.
- Tetra/cell 교집합 적분으로 potential/charge exchange를 통일했다. Uniform density, affine potential, total charge, work adjoint를 실제 GAA mesh에서 검사했다.
- Exact saved field의 physical fixed-point residual과 Gauss balance를 최종 판정에 추가했다.

## 독립 전류 검증

Vgs=.8/Vds=.1 V, dx=.5 nm, modes=3/6의 실제 GAA에서 finite-spectral bond current와 같은 finite-bath Hamiltonian의 독립 Schur/RGF transmission integral을 비교했다. 전체 I–V의 dx=.25 nm, modes=2/4 설정과 구분한다.
초기 energy panel .02/.01/.005 eV 각각에서 최대 상대 차이 {max((v['relative_mismatch'] for v in current.get('checks',[])),default=float('nan')):.3g}였다.
이 비교는 동일 접촉 모델의 계산 일치이며 contact model 자체의 실험적/정량 검증은 아니다.

## 남은 문제

- Broadening, mode truncation, finite S/D와 basis/grid 의존성이 실제 수렴 표에 남는다. 한 쌍의 agreement로 plateau를 선언하지 않는다.
- Connected S/D neutrality는 isolated source calibration과 별개다. Finite-bath 결과의 오차가 남는다.
- Abrupt raised-S/D 접합의 독립적인 mode truncation이 dx 의존성을 키운다. Nested trace 보강 비교는 dx=.5/.25 nm 두 계산에서 전류 차이 {nested_change*100 if nested_change is not None else float('nan'):.3g}%였으며, 두 수준만으로 plateau를 인증하지 않는다. baseline I–V와 섞지 않았다.
- Semi-infinite NEGF의 실제 equilibrium SCF는 수렴했고 I=0이다. 실제 biased GAA는 좁은 공명 적분 기준을 통과하지 못했다. `results/negf_lowbias/`의 ERROR/traceback을 보존한다.
- Fine quantum grid만으로 EMA dispersion, transverse truncation, Poisson mesh 또는 contact 모델의 검증이 끝나지 않는다.
- Hole/scattering/trap/tunnelling/XC/nonparabolicity 및 calibrated workfunction은 포함하지 않는다.

원시 CSV/NPZ/설정/잔차/VTK와 세부 [VALIDATION.md](VALIDATION.md), [HOW_TO_RUN.md](HOW_TO_RUN.md)를 함께 제공한다.
STEP01–05는 보존했고 원격 push와 학교 SLURM 제출은 수행하지 않았다.
'''
    (root/'STATUS.md').write_text(status)
    table=[]
    for dirname in CAMPAIGNS:
        p=result/dirname
        if not p.exists():continue
        for r in summarize(p)['cases']:
            c=read(p/(r['case']+'.json'));group=r['case'].rsplit('_',1)[0]
            setting=parameter_label(group,c)+'; '+PARAMETER_UNITS.get(group,group)
            table.append(f"| {r['case']} | {setting} | {float(r.get('id_a','nan'))*1e6:.6g} | {r['status']} | {float(r.get('phi_residual_v','nan')):.3g} | {float(r.get('charge_residual','nan')):.3g} | {float(r.get('sd_neutrality_relative_error','nan')):.3g} |")
    incomplete=[]
    for dirname,campaign in campaigns.items():
        for case in campaign['incomplete_cases']:
            last=case['last_saved_iteration']
            incomplete.append(f"| {dirname}/{case['case']} | {case['execution_status']} | {last['iteration'] if last else '—'} | {last['phi_residual_v'] if last else '—'} |")
    nested_table=[]
    for name,rs in nested.items():
        for row in rs:
            nested_table.append(f"| {name} | {float(row['id_a'])*1e6:.8g} | {row['status']} | {row['phi_residual_v']} | {row['charge_residual']} |")
    frozen_report=''
    if frozen.get('status')=='COMPLETED_FROZEN_CURRENT_ONLY':
        comparison=frozen['checks'][-1]
        frozen_report=f'''
## Actual-device frozen-field contact comparison

The actual finite-bath SCF field at Vgs=.8 V, Vds=.1 V, dx=.5 nm, modes=3/6 was held fixed. Semi-infinite continuation replaces the terminal x hard walls and distributed finite baths; the interior projected operator and potential are retained. No semi-infinite density or Poisson–quantum fixed point is claimed here.

| Contact / chemical potential | Current (uA) |
|---|---:|
| Finite bath, its calibrated mu={frozen['finite_bath_mu_ev']:.8g} eV | {frozen['finite_bath_bond_current_a']*1e6:.9g} |
| Semi-infinite leads, identical mu | {comparison['same_mu_current_a']*1e6:.9g} |
| Semi-infinite leads, separately neutral mu={frozen['infinite_neutral_lead_mu_ev']:.8g} eV | {comparison['neutral_lead_mu_current_a']*1e6:.9g} |

Three independently initialized adaptive energy panels (.02/.01/.005 eV) are saved in `results/frozen_contact_comparison/summary.json`. This separates a contact substitution from chemical-potential recalibration. It demonstrates a current difference on a fixed field, not agreement of two fully coupled contact models. Sampled transmission plots are distinct from the adaptive current integrals.
'''
    report='''# Physics and numerical validation

All numbers below come from saved tests or executed solver outputs. `quantitative_validation_pass` is false throughout. The three-level convergence tables are one-bias sensitivities; they do not certify the full geometry/bias domain.

## Scope and independently executed checks

| Check | Evidence |
|---|---|
| Geometry / 4-sided GAA / volumes / shared faces | Per-run geometry_checks.json and mesh/geometry.json |
| P1 FEM, not native FVM | Signed material-weighted tetra stiffness in one DEVSIM potential space |
| Homogeneous and layered capacitor | results/poisson_benchmarks/summary.json; three meshes, analytic potential/normal displacement |
| Uniform charged slab / Gauss / refinement | Same JSON; original source/control volume in actual DEVSIM equations |
| Shared Si/Ox/spacer potential | test_triple_junction.py; exactly one shared unknown, charged real geometry |
| Potential and positive charge exchange | Exact clipped cell/tetra integrals; marginals, affine precision and work tests |
| Full mode basis vs real 3D FD | Independently assembled variable-section real-space matrix in test_physics.py |
| Analytic subbands / six valleys and spin | Unchanged masses and explicit degeneracies; subband_convergence.csv |
| Finite spectral vs direct GF density | results/benchmarks/summary.json |
| Bond vs transmission current | Reduced direct inversion and actual-GAA post_scf_current/summary.json |
| Semi-infinite RGF vs dense inversion | Unequal block-size and complete small open-system density/current tests |
| Pauli / equilibrium / gauge / reversal | Finite matrix tests; no biased standing-state origin guessing |
| Numerical failure retention | Per-bias ERROR/NOT_CONVERGED and complete tracebacks; no interpolated current |
| Overlapping transfer/output biases | results/cross_sweep_consistency.json; independently continued saved states at the same three biases |

## STEP05 audit

Read-only reruns use its original saved density and unmodified source. See step05_audit.json, step05_interface_audit.json and associated log. Its reported grid/mode/length/basis sensitivities match the CSV audit. The material triple-line potential discontinuity was an additional code defect, separate from bath/model limitations. STEP06's shared-DOF solution fixes that defect.

## Fully coupled one-bias comparisons

Vgs=.8 V, Vds=.1 V. Potential residual is unmixed volts; charge residual is volume-weighted relative L2. Mode counts are per representative valley pair, with +/- multiplicity and spin applied separately once. Channel charge includes the suspended bridge and both spacer extensions; the reported barrier window is the gated 24 nm. Neutrality is the outer-half S/D number difference relative to donors, not net device neutrality. A gate-controlled transistor need not have zero total net charge.

| Case | Setting | Id (uA) | Status | Phi residual (V) | Charge residual | S/D neutrality |
|---|---|---:|---|---:|---:|---:|
'''+ '\n'.join(table)+'''

Plateau decisions and adjacent current/charge/gated-barrier differences are in each campaign's summary.json. Missing or failed levels fail the plateau test. Fixed/local comparisons are basis comparisons, not monotonic refinement levels. The executed transverse-grid and SCF-tolerance cases appear in the table; incomplete levels appear below. These comparisons remain restricted to one bias, rather than a validated off/threshold/on and drain-bias domain.

## Incomplete campaign cases

These cases have no final accepted point table. A checkpoint or an SCF iteration is not a completed bias. Interrupted computations can be resumed with the same source/configuration; planned but unstarted levels have no physical result. This table records the state when the report was generated.

| Case | Execution status | Last saved iteration | Last phi residual (V) |
|---|---|---:|---:|
'''+ ('\n'.join(incomplete) if incomplete else '| — | Every planned case has a completed result row; see its actual status | — | — |')+'''

## Separately executed junction-compatible basis

Low wide-section modes are enriched by zero-extended narrow-section modes and the original full operator is reprojected. The independent extension fingerprint is retained in each state. These results do not replace the baseline sweep or use fitted current/hopping factors.

| Case | Id (uA) | Status | Phi residual (V) | Charge residual |
|---|---:|---|---:|---:|
'''+ '\n'.join(nested_table)+f'''

The two-grid current difference is {nested_change*100 if nested_change is not None else float('nan'):.6g}% relative to the finer grid. This is a diagnostic improvement, with two levels and finite-bath/contact/transverse limitations still present; it is not a certified plateau. The semi-infinite biased attempts are separate: the first attempts have explicit quadrature ERROR records and `nested_negf_refined` retains its unconverged SCF history/checkpoint. None supplies an accepted biased reference.
'''+frozen_report+f'''

## Provisional transfer metrics

| Metric | Actual saved value | Definition / limit |
|---|---:|---|
| Vth | {metrics.get('vth_constant_current_v',float('nan')):.6g} V | Constant current 1e-7 A/um of gated perimeter; 0.034 um |
| Minimum sampled SS | {metrics.get('ss_min_mv_dec',float('nan')):.6g} mV/dec | Adjacent samples below threshold; not a validated subthermal device claim |
| Maximum gm | {metrics.get('gm_max_a_per_v',float('nan'))*1e6:.6g} uA/V | Numerical derivative of provisional transfer curve |
| Ion / Ioff | {metrics.get('ion_ioff',float('nan')):.6g} | Vgs=1.2 / 0 V, Vds=.1 V |

Definitions and absolute currents are in `results/iv/transfer_metrics.json`. These are properties of the executed numerical configuration, not calibrated experimental predictions. All sweep points retain their individual contact neutrality and dispersion flags.
'''+'''

## Reference limitations and next steps

1. Resolve biased semi-infinite density integration, near-dark states and convergence against independently refined real-energy panels; do not insert artificial eta to obtain a successful plot.
2. Establish the contact model before selecting a teacher: finite constant-Gamma tails and connected neutrality remain important. Increasing transverse basis size changes high-energy tail population and is not automatically a physical contact cure.
3. Abrupt section changes also require a shared retained trace space: results/junction_overlap.json quantifies the missing narrow-to-wide norm and its tx amplification. validation/junction_compatible.py implements and separately fingerprints a nested-subspace comparison; it does not replace the baseline I–V. Complete transverse mesh, mode/basis, Poisson, dx, contact length and energy convergence at off/threshold/on, equilibrium/.1/.5/1 V. Current panels above the actual FD dispersion-valid window remain provisional.
4. Only after a restricted numerical/reference domain is validated should dataset labels become teacher-eligible. No neural model was trained in STEP06.

## Reproduction and cost

All run folders include settings, source/config identity, resource estimate, mesh and raw SCF histories. Runtime per accepted/failed point is in its CSV. Requirements and environment are recorded. Finite dense spectral cost is cubic; local results do not establish speedup relative to RGF. School partition/resource settings are chosen by the user; no remote job was submitted.

The final delivery also includes the exact earlier simulator core in `results/source_snapshots/`, indexed by its original source digest. Final changes to export restoration and exclusion of private virtual environments from the code identity do not rewrite the source identity of existing physical states. The archive is independently re-hashed by `validation/package_audit.py`. Use a separate copy with that archived core to resume the original guarded calculations, or use a new output folder with the delivered code. Field re-exports are checked against the exact saved silicon potential; native DEVSIM coordinates are cm and quantum VTU coordinates are nm.
'''
    (root/'VALIDATION.md').write_text(report)
    return summary


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    figures(root);print(json.dumps(write_reports(root),indent=2))
