"""Build reports/figures exclusively from saved executions; no success interpolation."""
import sys,json,csv
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from provenance import source_digest

ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';F=R/'figures';F.mkdir(exist_ok=True)
def read(folder):return json.loads((R/folder/'summary.json').read_text())
def write(name,text): (ROOT/name).write_text(text.strip()+'\n')
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])
def fmt(x):return f'{x:.6g}'
def pct(x):return f'{100*x:.4g}%'
def save(name):plt.savefig(F/name,dpi=180,bbox_inches='tight');plt.close()

def main():
    u=read('uniform');lo=read('gaa_modes');hi=read('gaa_modes_batched');matched=read('matched_mu_modes');chain=read('open_chain');small=read('constriction');origin=read('origin_sensitivity');tests=read('tests_batched')
    # The high-mode run starts again at M20. Preserve the earlier M14->M20
    # refinement fields while replacing absolute quantities with the batched run.
    merged={}
    for r in lo['cases']+hi['cases']:
        key=(r['modes'],r['prescription'])
        merged[key]={**merged.get(key,{}),**r}
    cases=list(merged.values())
    cases.sort(key=lambda x:(x['modes'],x['prescription']))
    fixed=[r for r in cases if r['prescription']=='fixed_mu'];fixedn=[r for r in cases if r['prescription']=='fixed_N']
    # Export all mode quantities, preserving the complete nested JSON alongside.
    with (R/'mode_summary.csv').open('w',newline='') as f:
        keys=['modes','sd_modes','prescription','mu_ev','electrons','open_electrons','closed_open_relative_l2','linear_response_conductance_fixed_contact_mu_s']
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:r[k] for k in keys} for r in cases)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2})
    fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    ax[0,0].plot([r['modes'] for r in fixed],[r['electrons'] for r in fixed],'o-',label='Fixed contact chemical potential')
    ax[0,0].plot([r['modes'] for r in fixedn],[r['electrons'] for r in fixedn],'s--',label='Fixed N = 72.6');ax[0,0].set_ylabel('Total electrons');ax[0,0].legend(fontsize=8)
    ax[0,1].plot([r['modes'] for r in fixedn],[r['mu_ev'] for r in fixedn],'o-',color='#cf6b24');ax[0,1].set_ylabel('Fixed-N chemical potential (eV)')
    for subset,label in [(fixed,'Fixed mu'),(fixedn,'Fixed N')]:
        rr=[r for r in subset if 'density_change_from_previous' in r]
        ax[1,0].semilogy([r['modes'] for r in rr],[100*r['density_change_from_previous'] for r in rr],'o-',label=label)
    ax[1,0].axhline(2,color='grey',ls='--',label='2% exploration target');ax[1,0].set_ylabel('Successive density L2 change (%)');ax[1,0].legend(fontsize=8)
    ax[1,1].plot([r['modes'] for r in fixed],[1e6*r['linear_response_conductance_fixed_contact_mu_s'] for r in fixed],'o-');ax[1,1].set_ylabel('Sampled linear conductance (uS)')
    for a in ax.flat:a.set_xlabel('Requested channel modes (S/D = 2x)')
    fig.suptitle('Frozen GAA mode refinement — no complete density/current plateau');save('mode_refinement.png')
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for L in (12,24,48,96):
        d=np.load(R/f'uniform/L{L}_M300_fixed_mu.npz');ax[0].plot(d['x_nm'],d['line_density_per_nm']/3,label=f'L={L} nm')
    ax[0].axhline(1,color='black',ls='--');ax[0].set(xlabel='Position (nm)',ylabel='Closed density / neutral open density',xlim=(0,12));ax[0].legend(fontsize=8)
    rr=[r for r in u['cases'] if r['modes']==300 and r['prescription']=='fixed_mu']
    ax[1].semilogy([r['length_nm'] for r in rr],[r['center_relative_error_to_open'] for r in rr],'o-',label='Central 6 nm')
    ax[1].semilogy([r['length_nm'] for r in rr],[abs(r['end_region_neutrality']) for r in rr],'s-',label='Outer 3 nm');ax[1].set(xlabel='Finite wire length (nm)',ylabel='Relative error / depletion');ax[1].legend()
    fig.suptitle('Uniform 20 x 15 nm wire — same mu, full transverse basis');save('reservoir_boundary.png')
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for m,style in [(2,'--'),(20,'-')]:
        d=np.load(R/f'matched_mu_modes/M{m}_fixed_mu.npz');x=d['xyz_nm'][:,0];xx=np.unique(x);dx=xx[1]-xx[0]
        for name,col in [('density_cm3','#315ba8'),('open_density_cm3','#d27531')]:
            line=np.array([d[name][x==z].sum()*.5e-21/dx for z in xx]);ax[0].plot(xx,line,ls=style,color=col,label=f'{"Closed" if name=="density_cm3" else "Open"}, M={m}')
    ax[0].set(xlabel='Position (nm)',ylabel='Electron line density (1/nm)');ax[0].legend(fontsize=8)
    rr=[r for r in matched['cases'] if r['prescription']=='fixed_mu']
    for region in ('channel','extension','left','right'):
        ax[1].semilogy([r['modes'] for r in rr],[100*r['regional'][region]['closed_open_relative_l2'] for r in rr],'o-',label=region)
    ax[1].set(xlabel='Requested channel modes',ylabel='Closed/open density difference (%)');ax[1].legend(fontsize=8)
    fig.suptitle('Same frozen GAA field and chemical potential — boundaries remain different');save('matched_mu_density.png')
    f=R/'repaired_ms/failures/vg0.800000_vd0.050000_valley0_occupation.npz';d=np.load(f);p=np.load(R/'state_audit/vg0.800000_vd0.050000_valley0_occupation_profiles.npz')
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained');ids=p['state_indices'];order=np.argsort(d['occupation_uncertainty_per_state'][ids])[-4:]
    for j in order:ax[0].plot(p['x_nm'],p['longitudinal_probability'][j],label=f'State {ids[j]}')
    ax[0].set(xlabel='Position (nm)',ylabel='Probability per longitudinal slice');ax[0].legend(fontsize=8)
    a=d['active'];sc=ax[1].scatter(d['weights'][a,0],d['weights'][a,2],c=d['weights'][a,1],s=18,cmap='viridis');ax[1].plot([0,.5],[0,.5],'k--',lw=1);ax[1].set(xlabel='L projector weight',ylabel='R projector weight');fig.colorbar(sc,ax=ax[1],label='Channel weight')
    fig.suptitle('Rejected GAA iterate at Vg=0.8, Vd=0.05 V — no assigned injection');save('ambiguous_states.png')
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
    for folder,label in [('repaired_negf','Vg=0.8, Vd=0.05 V'),('repaired_negf_vg04','Vg=0.4, Vd=0.01 V')]:
        for file in (R/folder).glob('*trials.jsonl'):
            events=[json.loads(z) for z in file.read_text().splitlines()];accepted=[e for e in events if e.get('accepted')]
            if accepted:ax.semilogy([e['iteration']+1 for e in accepted],[e['residual_inf'] for e in accepted],'o-',label=label)
    ax.axhline(1e-4,color='grey',ls='--',label='Potential residual threshold');ax.set(xlabel='Accepted domain Newton step',ylabel='Original map residual (V)',title='Actual GAA trial recovery — runs interrupted, not converged');ax.legend(fontsize=8);save('domain_newton_trials.png')
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for gamma in (.0025,.005,.01,.02):
        rr=[r for r in chain['finite_bath_contact_sweep'] if r['gamma_ev']==gamma]
        ax[0].semilogy([r['electrode_cells'] for r in rr],[100*r['current_relative_to_open'] for r in rr],'o-',label=f'Gamma={gamma:g} eV')
        ax[1].semilogy([r['electrode_cells'] for r in rr],[100*r['channel_density_relative_to_open'] for r in rr],'o-')
    ax[0].set_ylabel('Current difference from open lead (%)');ax[1].set_ylabel('Channel density difference (%)')
    for a in ax:a.set_xlabel('Cells in each finite bath')
    ax[0].legend(fontsize=8);fig.suptitle('Finite-bath model: contact parameters affect density and current');save('finite_bath_convergence.png')
    fig,ax=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for panel in (.02,.01,.005):
        rr=[r for r in chain['open_chain'] if r['energy_panel_ev']==panel]
        ax[0].plot([r['vd_v'] for r in rr],[1e6*r['current_a_per_degeneracy'] for r in rr],'o-',label=f'Panel {panel:g} eV')
        ax[1].semilogy([r['vd_v'] for r in rr],[r['density_max_abs_error'] for r in rr],'o-')
    ax[0].set(xlabel='Applied chain voltage (V)',ylabel='Current per degeneracy (uA)');ax[0].legend(fontsize=8)
    ax[1].set(xlabel='Applied chain voltage (V)',ylabel='Density max error vs dense integration')
    fig.suptitle('Independent open-chain reference — frozen ramp, not GAA I-V');save('open_chain_reference.png')
    mt=table(['Channel / S-D','N (fixed mu)','mu (fixed N, eV)','density refinement (fixed mu)','channel charge refinement (fixed mu)'],
        [[f"{a['modes']}/{a['sd_modes']}",fmt(a['electrons']),fmt(b['mu_ev']),pct(a['density_change_from_previous']) if 'density_change_from_previous' in a else '—',pct(a['channel_charge_change_from_previous']) if 'channel_charge_change_from_previous' in a else '—'] for a,b in zip(fixed,fixedn)])
    write('MODE_CONVERGENCE.md',f'''# STEP08 — Mode convergence

**GAA 전체 수렴은 미달성이다.** M1–M5를 넘어 14/28, 20/40, 28/56, 36/72까지 실제 frozen-field 계산을 수행했다. 마지막 단계에서도 전체 density L2 변화는 {pct(fixed[-1]['density_change_from_previous'])}, transmission sampled L2 변화는 {pct(fixed[-1]['transmission_relative_l2_change'])}다. 채널 charge만 안정되는 것을 전체 검증으로 취급하지 않았다.

## Finite-N STEP07 평형장을 고정한 결과

모든 비교는 같은 xyz 격자, 동일한 frozen phi, 같은 온도다. fixed-mu 열의 mu={fixed[0]['mu_ev']:.12g} eV는 해당 finite-N 평형장의 저장값이며 neutral open-lead mu와 다르다. fixed-N은 72.6 electrons를 유지하도록 mu를 다시 푼다.

{mt}

실제 S/D rank는 요청한 sd_modes에 narrow trace를 보강하므로 더 클 수 있다. 정확한 valley별 dimension은 각 summary.json에 있다. 고정 mu에서 N이 증가하는 것 자체는 오류가 아니다. 같은 H의 더 넓은 공간을 채우는 결과이고, 별도로 mu를 움직이는 fixed-N 비교와 구분해야 한다.

`linear_response_conductance_fixed_contact_mu_s`는 저장한 221개 T(E) sample의 thermal integral 진단값이다. 좁은 resonance에 대한 독립 adaptive 검증을 통과한 비평형 전류가 아니며, fixed-N 행에도 비교용 같은 contact mu의 값이 들어 있다. 이를 실제 GAA I–V로 사용하지 않는다. Lowest local subbands, occupied eigenstate counts, regional charges, T(E), closed/open densities도 저장했다. 공간 전체 LDOS(E,r)와 충분한 GAA occupied-subband plateau는 아직 인증하지 않았다.

## Uniform contact DOS

같은 full-basis neutral mu={u['full_basis_neutral_mu_ev']:.12g} eV에서 4 modes는 {pct(u['bulk_mode_sweep'][0]['relative_missing_density'])}, 20 modes는 {pct(next(r for r in u['bulk_mode_sweep'] if r['modes']==20)['relative_missing_density'])}의 bulk line charge를 놓친다. 64/100/300 modes는 각각 {', '.join(pct(r['relative_missing_density']) for r in u['bulk_mode_sweep'] if r['modes'] in (64,100,300))}. 이것은 rectangular 20x15 nm uniform lead의 결과이며 GAA의 contact-mode 기준을 자동 결정하지 않는다.

## Independent real-space limit

`results/constriction`은 dx=1,0.5,0.25 nm × electrode length=2,4,8 nm × 5 bases = {len(small['cases'])} cases이다. 독립 3D stencil로 만든 H와 `U†HU`, 고유값, closed density, reconstruction, dense real-space open inversion의 T/current를 비교했다. Full-basis identity 통과={small['full_basis_identity_pass']}. Truncated-basis 오차는 삭제하지 않았다. Nested trace 보강은 projection completeness와 occupied-density convergence를 같은 것으로 만들지 않는다.

## 메모리 수정

Unbatched M28 실험은 8 GiB cgroup에서 OOM(exit -9)로 중단됐다. 결과를 성공으로 처리하지 않았다. Contour pole 8개, transmission energy 16개씩 처리하도록 수정한 뒤 M28과 M36 계산을 완료했다. Pole 수, self-energy, tolerance는 그대로다. Batch 1과 batch 96/200 비교를 테스트했다. 고차 GAA 전체 고유벡터와 real-axis adaptive panel 배열은 여전히 비용이 크다.

![mode](results/figures/mode_refinement.png)

Raw: `results/gaa_modes`, `results/gaa_modes_batched`, `results/matched_mu_modes`, `results/mode_summary.csv`. `gaa_modes_high`는 실패한 unbatched 실행이다. GAA dx/transverse/Poisson joint convergence는 미실행이며 `HPC_PLAN.md`를 따른다.
''')
    uniform_table=table(['L (nm)','center relative error','outer 3 nm neutrality error'],[[r['length_nm'],fmt(r['center_relative_error_to_open']),pct(r['end_region_neutrality'])] for r in u['cases'] if r['modes']==300 and r['prescription']=='fixed_mu'])
    write('CONTACT_CONVERGENCE.md',f'''# STEP08 — Contact and reservoir physics

**전극 길이를 늘리는 것과 hard wall 옆 reservoir가 bulk-like해지는 것은 다르다.** 같은 uniform wire의 중앙은 열린 전극 밀도로 수렴하지만, 경계 바로 옆의 고정된 구간은 depletion을 유지한다.

{uniform_table}

이 실험은 fixed uniform potential이며 gradient=0은 입력 조건이다. SCF screening을 증명한 값이 아니다. L=12/24/48/96 nm와 mode=4/8/12/16/20/32/48/64/100/300, fixed mu/fixed N을 모두 계산했다. LDOS는 `-df/dE`로 thermal convolution한 delta spectrum이다. GF에 인위적인 eta를 넣지 않았다. Transverse spacing 1/0.5/0.25 nm에서 full-basis bulk density와 neutrality mu도 비교했다.

STEP07의 약 -48% reservoir charge 오차를 단일 원인으로 단정할 수 없다. 불충분한 occupied modes, 끝 hard wall, finite-N mu와 neutral-contact mu 차이, 충분하지 않은 screening이 함께 작용한다. 전체 N을 맞춰도 지역별 전하가 맞는 것은 아니다. 또한 nodal donor quadrature 72.840918878과 quantum-cell donor 72.6의 차이는 서로 다른 source discretization의 오차이며, 이를 fitting으로 renormalize하지 않았다.

Finite WBL path C는 4 electrode lengths × 4 gamma 값을 독립 small chain에서 비교했다. L=64의 current error는 gamma .0025/.005/.01/.02 eV에서 각각 {', '.join(pct(r['current_relative_to_open']) for r in chain['finite_bath_contact_sweep'] if r['electrode_cells']==64)}. 특정 gamma의 일치를 열린 전극의 물리적 인증으로 삼지 않는다. 채널 density 오차도 동시에 기록했다.

Energy-dependent embedding 대안은 positive-residue finite-chain pole expansion을 평가했다. 8/16/32/64/128 poles의 complex-axis self-energy 오차는 `results/uniform/summary.json`에 있다. 이 결과는 equilibrium contour representation의 가능성만 보여 준다. Biased GF-free density, initial preparation, long-time/continuum limit을 해결한 구현은 아직 없다.

![boundary](results/figures/reservoir_boundary.png)
![finite bath](results/figures/finite_bath_convergence.png)

실제 GAA electrode length × mode × mu의 fully coupled plateau는 미달성이다. Uniform/reduced-wire 결과로 이를 대신 주장하지 않는다.
''')
    region_table=table(['Region','M2 closed/open difference','M20 closed/open difference'],[[name,pct(matched['cases'][0]['regional'][name]['closed_open_relative_l2']),pct([r for r in matched['cases'] if r['prescription']=='fixed_mu'][-1]['regional'][name]['closed_open_relative_l2'])] for name in ('channel','extension','left','right')])
    write('MS_NEGF_COMPARISON.md',f'''# STEP08 — MS versus open NEGF

**동일 mu의 실제 GAA density 차이는 해결되지 않았다.** 원래 STEP07의 independently coupled equilibrium 차이 24.0657%와 frozen-field 차이 24.7878%는 서로 다른 비교다. 이번에는 저장된 fixed-mu MS 평형장에 두 density solver를 적용해 frozen 24.7878%를 재현했다. M20/40에서도 {pct([r for r in matched['cases'] if r['prescription']=='fixed_mu'][-1]['closed_open_relative_l2'])}가 남는다.

{region_table}

1. Uniform wire: 같은 H dispersion, full modes, T, mu에서 중앙 밀도의 finite-length limit을 확인했다. 경계부 전체를 같다고 주장하지 않는다.
2. Small constriction: 같은 potential, full real-space reference에서 projection/density/T/current identities를 확인했다. Length와 dx를 분리한 45 case를 계산했다.
3. Actual frozen GAA: fixed-N field와 fixed-contact-mu field를 별도로 사용했다. mode sweep에서 각각 fixed mu와 fixed N을 분리했다.
4. Fully coupled GAA: STEP07의 낮은 모드 평형들을 warm-start 재현했다. 고차 모드의 독립 SCF 비교와 비평형 인증은 미완료다.

## Occupation and tracking

실패 파일의 모든 valley를 저장하도록 바꿨고, 0.01/0.05/0.1 V에서 active ambiguous state 186건의 energies, weights, profile, localization, nearest gap, reach, contrast, boundary stability를 추출했다. `state_audit/active_states.csv`와 NPZ가 근거다. 이것은 실패한 iterate의 상태이지 SCF 해가 아니다.

추가로 `phi(vd)=phi_eq+vd*x/L`의 **명시적 frozen family**에서 M2/M4/M10과 7 bias를 계산했다. 이 controlled family는 actual Poisson bias continuation이 아니다. overlap tracking과 degenerate projector 비교를 수행해도 origin은 별도 미결정으로 남긴다. 

대칭 uniform chain은 active {chain['standing_state_counterexample']['active_states']}개가 모두 ambiguous인데 open T={chain['standing_state_counterexample']['open_transmission_at_mu']:.12g}이다. 이는 closed standing-state probability만으로 일반적인 injection origin을 결정할 수 없다는 반례다. 원래 MS-DFT 전체 이론이 틀렸다는 주장은 아니다.

![comparison](results/figures/matched_mu_density.png)
![states](results/figures/ambiguous_states.png)
''')
    write('ROOT_CAUSE_REPORT.md',f'''# STEP08 — Root-cause report

프로젝트: **DEVSIM-QT**. 이 패키지는 원인 분석과 검증 가능한 수정을 포함하지만, **정량 검증된 비평형 GAA I–V / Fully-SCF Teacher는 아직 완성되지 않았다.**

## 원인과 실제 수정

| 문제 | 분류 | 근거 | 이번 조치 |
|---|---|---|---|
| unknown origin을 number equation에서 L로 계산 | Algorithmic error | STEP07 `multispace/constraints.py`의 `origin==2` 외 모두 fL | N_low/high, mu interval을 먼저 계산하고 material ambiguity는 density 생성 전에 거부 |
| index 기반 상태 비교 | Algorithmic weakness | eigenvalue crossing/degenerate rotation 가능 | Hungarian overlap matching, degenerate-subspace singular values; origin 결정과 분리 |
| invalid Newton trial이 전체 SCF를 즉시 abort | Algorithmic error | STEP07 failure trace, 실제 GAA bound audit | 유효한 domain에서만 Armijo 감소 스텝 허용; half-step 성공, invalid derivative probe와 stalled search 기록 |
| occupied transverse modes 부족 | Numerical approximation | uniform 4-mode missing charge {pct(u['bulk_mode_sweep'][0]['relative_missing_density'])}; actual GAA 36-mode density도 변화 | M1–M5 및 M14/M20/M28/M36, fixed-mu/fixed-N 분리 실험 |
| finite wall와 열린 전극 차이 | Boundary/model + finite-size approximation | 중앙은 길이에 수렴하나 outer 3 nm depletion 유지 | uniform mode/length/mu 결합 실험; region별 density comparison |
| standing eigenstate의 injection origin 부재 | Physical applicability limit | symmetric wire: T=1인데 active 19개 모두 ambiguous | threshold 완화나 임의 filling 없이 거부; open path 별도 유지 |
| isolated bound occupation | Physical model incompleteness | Vg=.8,Vd=.001 initial trial의 pole E=0.0052550675 eV, fL−fR=0.00961798 | audit 보존; final solution에서도 유효해야 함. 준비/relaxation 모델 임의 추가하지 않음 |
| GF pole/energy 배열의 과도한 메모리 | Implementation resource problem | 8 GiB에서 M28 unbatched OOM | energy/pole batching; 동일 수치 확인, M36까지 완료 |
| current tolerance가 acceptance에 미사용 | Algorithmic omission | STEP07 acceptance 조건 | A의 후처리는 별도, B/C는 연속 current change도 기존 tolerance로 확인 |

## 수정 전후 재현

Original STEP07 source digest `c0fbe6dc...`를 별도 복원해 실행했다. 이전 density를 초기값으로 사용한 **warm-start 재현**이며 independent cold-start 수렴이라고 주장하지 않는다. Vg=.8 equilibrium MS/NEGF가 다시 수렴했고, MS Vd=.01/.05/.1은 다시 거부됐다. .05의 첫 valley에서 기존 11개 ambiguous state도 재현했다. Finite spectral Vg=.8,Vd=.1은 약 4.49181 uA로 재현됐다.

전자수 수정 후 .05 V의 feasible mu interval은 [-0.04742291141,-0.04686515073] eV이고 audit mu에서 전체 N 범위는 [72.04560806,73.16232108] electrons이다. 1.11671302 electrons의 불확실성을 먼저 드러낸다. 이를 제거했다고 주장하거나 arbitrary L로 맞추지 않는다.

NEGF domain Newton의 실제 GAA 로그에서는 invalid full step 뒤 valid half-step이 채택됐다. 두 audit run은 계산 비용 때문에 중단했고 `INTERRUPTED.json`으로 표시했다. 이전의 즉시-abort 구현 오류는 수정됐지만, 그 실행들을 최종 fixed point 성공/실패로 판정하지 않았다. .001/.01의 초기 trial에서 pole이 발견된 경우 역시 모든 해의 비존재 증명이 아니다.

## 최종 연구 질문에 대한 답

1. **실패는 복합 원인이다.** 구현 오류, 부족한 mode/contact convergence, closed-origin approximation의 물리적 한계를 각각 확인했다.
2. **Uniform wire 내부는 수렴한다. 실제 GAA 전체는 미인증이다.** finite boundary와 open reservoir를 무조건 동일시하면 안 된다.
3. **검사한 GAA 비평형 MS 점유는 안정적으로 정의되지 않았다.** 모드 증가와 tracking만으로 해결되지 않는다.
4. **정당한 이식 범위는 projector/number/gauge/density bookkeeping이다.** SIESTA의 orbital averaging/tie-break/XC/variational constrained search까지 EMA가 구현했다고 할 수 없다.
5. **Open NEGF는 small-system 독립 reference를 통과한다. GAA reference 완성은 남았다.** bound preparation, embedded dark states, joint convergence가 아직 필요하다. C는 별개의 finite bath model이다.
6. **검증된 정상 GAA I–V는 이번에도 제출할 수 없다.** 전제인 density/mode/contact plateau가 통과하지 않았으므로 100개 bias sweep으로 이를 숨기지 않았다. 지원 코드와 조건부 실행 명령은 제공한다.

## 남은 구체적 한계

Bound audit는 모든 propagating lead band 아래의 poles를 검사한다. Embedded dark states의 전역적인 배제, above-band occupied poles, 전체 GAA density-matrix PSD는 완전히 검증하지 않았다. Finite real-space quadrature의 매우 좁은 resonance는 별도 integration budget 문제다. 원래 donor/charge source discretization 차이와 GAA mesh refinement도 미해결이다. 기존 좁은 mode configuration에서 얻은 수렴값은 정량적인 reference가 아니다.

원본 SIESTA 새 첨부는 전달 실패했다. 이전 버전만 분석했으며 소스 비교 범위는 `MSDFT_SOURCE_COMPARISON.md`에 명시했다. 모든 단계가 완료됐다는 표현이나 AI training reference 인증은 하지 않는다.

![domain](results/figures/domain_newton_trials.png)
''')
    write('VALIDATION.md',f'''# STEP08 — Validation status

**Unit/native tests {tests['total_tests']} passed; GAA physical certification NOT passed.**

| 검증 | 실행 범위 | 판정 |
|---|---|---|
| Quantum/numerical tests | {tests['groups'][0]['tests']} tests | pass |
| 실제 DEVSIM/Gmsh | {tests['groups'][1]['tests']} tests, 별도 process | pass |
| Number bounds / degenerate tracking / invalid trial recovery | 신규 tests | pass |
| Batch-size invariance | 1 vs 96 poles, 1 vs 200 energies | pass |
| Uniform wire | 4 lengths × 10 mode levels × 2 number prescriptions | completed; 내부 수렴 확인 |
| Small constriction | 3 dx × 3 lengths × 5 bases | full-basis identity pass |
| Actual GAA frozen field | M2–M36, two number prescriptions | 전체 plateau 미달성 |
| Same-mu frozen GAA | M2–M20 | closed/open density 차이 남음 |
| Actual MS ambiguous state | 186 records + 3 modes × 7 frozen biases | origin unresolved |
| Open chain density/current | Vd=0/.001/.01/.05/.1 × 3 energy refinements | independent dense reference pass |
| GAA domain Newton | valid half-step recovery 기록 | interrupted, biased solution 미인증 |
| GAA equilibrium / C reference | warm-start 재현 | numerical convergence only |
| Grid/contact coupled GAA convergence | 계획 생성 | 미실행 |
| 정상 GAA I–V / AI teacher | prerequisites failed | 미제출 / 사용 불가 |

Open chain independent density max absolute error={max(r['density_max_abs_error'] for r in chain['open_chain']):.6g}, current absolute error={max(r['current_abs_error_a'] for r in chain['open_chain']):.6g} A per degeneracy. Current conservation과 full dense rho PSD도 확인했다. 이 작은 시스템 검증을 GAA certification으로 확장하지 않는다.

Test logs는 `results/tests_batched/`, 이전 core version 검증은 `results/test_report/`에 있다. Native/quantum 그룹은 분리 실행하며 JUnit completion XML을 요구한다. 예전 monolithic native 종료 문제를 성공으로 세지 않았다.

Final production digest: `{source_digest()}`. Snapshot은 `results/source_snapshots/`. 원본 source snapshot, 초기 seed, source/config identities와 실패 trace를 포함했다. 검증 script 자체도 최종 package manifest로 고정한다.

![chain](results/figures/open_chain_reference.png)
''')
    status={'step':8,'project_name':'DEVSIM-QT','source_sha256':source_digest(),'tests_passed':tests['total_tests'],
        'mode_levels_completed':[r['modes'] for r in fixed],'biased_gaa_certified_points':0,
        'quantitative_validation_pass':False,'ready_for_training':False,'full_task_physics_goal_complete':False,
        'siesta_latest_attachment_read':False,'previous_siesta_archive_analyzed':True}
    (ROOT/'STATUS.json').write_text(json.dumps(status,indent=2))
    write('STATUS.md','# STEP08 — status\n\n실행 상태는 `STATUS.json`, 과학적 판단은 `ROOT_CAUSE_REPORT.md`에 기록했다.\n\n64 tests passed. Mode/contact/GAA I–V certification remains incomplete. All recorded GAA quantitative_validation_pass flags are false.')
    print('reports and figures complete',flush=True)

if __name__=='__main__':main()
