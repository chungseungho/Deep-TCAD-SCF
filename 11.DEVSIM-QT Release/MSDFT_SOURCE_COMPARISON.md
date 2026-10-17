> Historical baseline document retained for research provenance. Its results and paths are not STEP11 validation evidence. See RESULTS.md and RELEASE_CERTIFICATE.json for this candidate.

# STEP08 — 연구실 MS-DFT 소스와 EMA 구현의 구분

이번 단계의 이름은 **DEVSIM-QT**, 폴더는 `8.DEVSIM-QT/`이다. STEP07은 수정 전 비교 기준이다.

## 확인한 자료와 한계

이번 `siesta-5.4.2.tar(1).gz` 첨부는 전달 실패로 읽지 못했다. 대신 이전 STEP07에서 제공받았던 `siesta-5.4.2.tar.gz`를 복원해 실제 Fortran을 다시 읽었다. SHA256은 `43f05274cd94a7a966a0164e11ad9b4b5576d987925844e00289b580a5b33fc0`. 새 첨부와 동일한지는 확인하지 않았다. 이 분석은 **복원된 버전**에만 해당한다. 파일별 해시와 줄 수는 `research/siesta_source_provenance.json`에 있다. 원본 Fortran, 압축파일, 실행 바이너리는 결과물에 포함하지 않았다. SIESTA를 컴파일하거나 실제 입력 FDF로 실행하지 않았다.

## 실제 호출 경로

| 단계 | 실제 소스 근거 | 확인된 의미 |
|---|---|---|
| Hamiltonian / Hartree / XC | `Src/dhscf.F`, `dhscf` 534행부터; Hartree 1526행부터, generalized Poisson 1600–1669행 | 밀도에서 Hartree와 XC 항을 구성. EMA/Hartree에는 원자궤도, pseudopotential, XC가 없다. |
| Density solver 선택 | `Src/compute_dm.F`, `compute_dm`, 206행 및 271행 `diagon`, 284행 `transiesta` | 대각화와 TranSIESTA는 별도 분기이다. |
| MS 대각화 선택 | `Src/diagon.F` 243행 `diagg_msdft`, 303행 `diagk_msdft` | Gamma / k점 분리. 178–183행 collinear 및 캐시/병렬 조건도 있다. |
| 일반화 고유문제 | `Src/diagg_msdft.F` 대각화 부분 및 `Src/diagk_msdft.F` 223행 `cdiag(Haux,Saux,...)` | `Hc = ESc`, overlap metric 사용. EMA 격자는 직교 기저다. |
| 분할 진단 | `diagg_msdft.F` 215–223행; `diagk_msdft.F` 230–239행 | 선택적 진단/Mulliken 분할. |
| 점유 설정 | `diagg_msdft.F` 266행, `diagk_msdft.F` 266행 `msdft`; `msdft.F` 507–568행 dispatch | MS 전용 점유와 기존 ground-state 점유를 구분한다. |
| 밀도행렬 | `diagg_msdft.F` 294–355행; `diagk_msdft.F` 414–503행 | 점유 가중치와 고유벡터 외적으로 D를 구성한다. |
| Mixing | `Src/m_mixing_scf.F90`, `mixers_scf_init` 46행부터 | DM/SCF mixer 설정, weight/history/method. 점유의 물리적 유일성을 보장하는 루틴이 아니다. |
| 수렴 | `Src/scfconvergence_test.F` 110–139행 | 사용자가 선택한 DM, H, 에너지 조건의 논리곱. 실제 tolerance는 FDF를 읽어야 알 수 있다. Harris 특수경로 198행 이후는 별도다. |
| 열린 전극 | `Src/m_ts_electrode.F90`, `SSR_sGreen_DOS` 46행, `SSR_sGreen_NoDOS` 428행, `create_Green` 709행 | 전극 surface Green function. |
| Self-energy | `Src/m_ts_elec_se.F90`, `UC_expansion_Sigma` 209행, `UC_expansion_Sigma_GammaT` 279행 | 에너지 의존 Σ와 Γ. |
| 열린 밀도 | `Src/m_ts_dm_update.F90`, `update_DM` 457–663행 | TranSIESTA 경로의 밀도 갱신. |
| 전류 후처리 | `Util/TS/TBtrans/m_tbt_tri_scat.F90`, `A_Gamma`, `orb_current` 1511–1623행; `m_tbt_save.F90` | 별도 TBtrans 계산. `msfermid`가 자동으로 SCF 뒤에 self-energy와 전류를 호출하는 것은 아니다. |

## `msfermid`의 실제 규칙

`Src/msdft.F` 571–1209행이 핵심이다.

* **legacy weight (740–839행):** 영역의 계수 절댓값 제곱을 더한 뒤 그 영역의 orbital 수로 나눈다. 따라서 `sum_R w_R = 1`인 양의 공간 projector 확률과 다르다. 영역의 궤도 수를 바꾸면 비교 자체가 바뀐다.
* **분류 (812–839행):** 가장 큰 영역이 L이면 L, R이면 R; 채널 우세 상태는 L/R weight를 다시 비교해 C-left 또는 C-right로 나눈다. R 우세 분기의 `>=`, 채널 내부의 `L >= R`는 동률을 결정하는 코드 선택이다. 거의 축퇴한 부분공간에 대해 회전 불변인 주입 원점을 추론하는 처리는 이 점유 경로에서 확인되지 않았다.
* **Mulliken 옵션:** `m_msdft_mulliken_partition.F90` 246–289행은 `Re[c_i* (Sc)_i]`를 영역별로 합친다. 전체 합은 S-normalization과 연결되지만 개별 영역값이 음수가 될 수 있다. 311–333행의 비교 규칙 역시 분류 규칙이지 scattering injectivity가 아니다.
* **global μ:** 895–920행에서 전체 전자수 `qtot`을 먼저 맞춘다. `tol=1e-10`, `nitmax=200`은 629–632행. 이는 총 전자수 오차이지 국소 reservoir neutrality 오차가 아니다.
* **비평형 offset:** 924–1075행에서 voltage window 내부의 점유를 고치고 공통 `ef_mod`를 다시 조정한다. `BiasReference=left`이면 `μL=ef_mod`, `μR=ef_mod−VSD*(eV→Ry)`이다. `ry2ev`라는 변수 이름과 달리 662행 값은 `0.073498618`이다. 대칭 기준 옵션도 별도로 존재한다.
* **window:** equilibrium μ 중심으로 bias의 2.5배 범위를 정하고, 바깥 점유는 이전 equilibrium 계산값을 유지한다. EMA 구현의 모든 상태 Fermi filling과 정확히 동일하지 않다.
* **추가 옵션:** weighted/half/eta-half/fixed-state 설정이 존재한다. 470–487행 weighted 루틴은 음의 Mulliken weight를 clip하고 분모가 작으면 반반 점유한다. 옵션의 존재를 물리적 정당화나 논문 검증으로 취급하지 않았다. STEP08에는 이러한 우회 filling을 도입하지 않았다.
* **전자수와 neutrality:** `qtot` 만족은 전체 조건 하나다. L/R 각각의 donor charge와 electron charge를 같게 만드는 두 독립 제약은 아니다.

## EMA에 적용한 것과 적용하지 않은 것

적용: 서로 다른 L/R chemical potential, 양의 직교 공간 projector, 전체 전자수 조건과 고정 contact μ의 구분, `rho = U f U†`, 일관된 electrostatic gauge, SCF 이후 열린 전극 전류 비교.

수정: STEP07의 미분류 상태를 왼쪽으로 간주하던 전자수 방정식. 이제 각 미분류 상태의 `[min(fL,fR),max(fL,fR)]`로 `N_low(μ),N_high(μ)`를 만들고, `N_high(μ_lo)=Ntarget`, `N_low(μ_hi)=Ntarget`의 허용 μ 구간을 먼저 계산한다. 그 구간의 중간값은 **불확실성 평가 지점**이며 물질적으로 미분류된 점유나 밀도를 반환하지 않는다. 두 reservoir μ를 평균해 점유하는 알고리즘이 아니다. 허용한 작은 tail 차이는 전자수 오차 경계로 공개한다.

추가: wavefunction overlap의 최대 대응 및 축퇴 부분공간 principal angles. 상태 추적은 분류와 별도이며 주입 원점을 만들어 주지 않는다. 고유값 index만으로 이전 상태의 점유를 복사하지 않는다.

적용하지 않음: per-orbital 평균을 EMA 확률로 복사, 임의 tie-break, 임의 회전/반반 점유, DFT constrained-search variational energy를 구현했다고 주장하는 것. 현재 A 경로는 **Multi-Space-inspired EMA/Hartree approximation**이다.

## 문헌과 구현의 경계

Lee et al., Advanced Science (2020), DOI [10.1002/advs.202001038](https://doi.org/10.1002/advs.202001038)는 실제 금속 electrode의 voltage-drop interface를 reservoir와 구분한다. PNAS QFL 연구 DOI [10.1073/pnas.1921273117](https://doi.org/10.1073/pnas.1921273117)는 채널에만 국소화된 bias-window state가 없다는 조건을 사용한다. 해당 보충자료를 이번에 새로 확보했다고 주장하지 않는다.

Bound-state의 준비 이력 문제는 Stefanucci, PRB 75, 195115 (2007), [원 논문](https://arxiv.org/abs/cond-mat/0608401)을 참고했다. 이 일반적인 사실을 근거로 어느 GAA trial 하나가 실패했다고 모든 fixed point가 없다고 결론내리지 않았다. 에너지 적분의 좁은 resonance 문제와 물리적 isolated-state population 문제는 별도로 기록한다.
