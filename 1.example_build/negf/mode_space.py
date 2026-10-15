"""Si effective-mass, slice-local transverse mode basis and ballistic RGF NEGF.

x transport, yz confinement, transverse hard-wall (Dirichlet), 3 Delta
valley orientations, valley multiplicities (2), spin degeneracy (2).
Energy units eV; mesh units nm; density output cm^-3. No EMA fitting.

This is an *energy-discretized demonstration solver*, not yet a validated
production-grade nanosheet model. Inelastic scattering is not included.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import diags, eye, kron

Q = 1.602176634e-19
H_PLANCK = 6.62607015e-34
KB_EV = 8.617333262145e-5
HB2_OVER_2M0 = 0.0380998212  # eV nm^2


@dataclass
class QuantumResults:
    current_a: float
    electron_density_cm3: np.ndarray
    spectrum_ev: np.ndarray
    transmission: np.ndarray
    valley_currents_a: dict


def fermi(e, mu, kt):
    return 1.0 / (1.0 + np.exp(np.clip((e-mu) / kt, -100, 100)))


def surface_self_energy(e, mode_energy, t, eta):
    """Semi-infinite 1D lead self-energy per mode, retarded branch."""
    z = complex(e, eta)
    a = z - (mode_energy + 2*t)
    root = np.sqrt(a*a - 4*t*t)
    s1, s2 = (a-root)/2, (a+root)/2
    return np.where(np.abs(s1) <= np.abs(s2), s1, s2)


def transverse_hamiltonian(cfg, valley, phi_yz):
    ny, nz = cfg.ny, cfg.nz
    dy = cfg.width_nm / (ny+1)
    dz = cfg.thickness_nm / (nz+1)
    ty = HB2_OVER_2M0 / (valley.my * dy * dy)
    tz = HB2_OVER_2M0 / (valley.mz * dz * dz)
    ly = diags([-np.ones(ny-1), 2*np.ones(ny), -np.ones(ny-1)],
               [-1,0,1]).tocsr() * ty
    lz = diags([-np.ones(nz-1), 2*np.ones(nz), -np.ones(nz-1)],
               [-1,0,1]).tocsr() * tz
    h = (kron(ly, eye(nz)) + kron(eye(ny), lz)).toarray()
    # electron potential energy E_c - q*phi, numerically -phi in eV for V.
    h += np.diag(cfg.ec0_ev - phi_yz.ravel())
    return h


def build_mode_slices(cfg, valley, phi_3d):
    nx, nm = cfg.nx, cfg.modes
    e = np.empty((nx, nm))
    u = np.empty((nx, cfg.ny*cfg.nz, nm))
    for i in range(nx):
        evals, evecs = eigh(transverse_hamiltonian(cfg, valley, phi_3d[i]),
                            subset_by_index=(0,nm-1), check_finite=False)
        e[i], u[i] = evals, evecs
        if i:
            for m in range(nm):
                if np.dot(u[i-1,:,m], u[i,:,m]) < 0:
                    u[i,:,m] *= -1
    dx = cfg.length_nm / (nx-1)
    t = HB2_OVER_2M0 / (valley.mx * dx * dx)
    onsite = np.array([np.diag(ei + 2*t) for ei in e], complex)
    hop = np.array([-t*(u[i].T @ u[i+1]) for i in range(nx-1)], complex)
    return e, u, t, onsite, hop


def spectral_rgf(onsite, hop, t, contact_e0, contact_eN, energy, eta,
                 with_density=True):
    """Return transmission and optional diagonal-slice lesser density matrix.

    gR is right-connected; gL is left-connected. We calculate source/drain
    injected propagators G_{i,0} and G_{i,N-1}; therefore no full inverse.
    """
    ns, nm = onsite.shape[:2]
    z = complex(energy, eta)
    I = np.eye(nm, dtype=complex)
    sig_l = surface_self_energy(energy, contact_e0, t, eta)
    sig_r = surface_self_energy(energy, contact_eN, t, eta)
    gamma_l = np.maximum(0, -2*np.imag(sig_l))
    gamma_r = np.maximum(0, -2*np.imag(sig_r))
    left = np.diag(sig_l)
    right = np.diag(sig_r)

    gr = np.empty((ns,nm,nm), complex)
    gr[-1] = np.linalg.inv(z*I-onsite[-1]-right)
    for i in range(ns-2,-1,-1):
        gr[i] = np.linalg.inv(z*I-onsite[i]-hop[i] @ gr[i+1] @ hop[i].conj().T)
    g00 = np.linalg.inv(z*I-onsite[0]-left-hop[0] @ gr[1] @ hop[0].conj().T)
    Gi0 = np.empty_like(gr)
    Gi0[0] = g00
    for i in range(1,ns):
        Gi0[i] = gr[i] @ hop[i-1].conj().T @ Gi0[i-1]

    # For G_{0,N}, start at G_{N,N} and step backwards with left-connected gr.
    gl = np.empty_like(gr)
    gl[0] = np.linalg.inv(z*I-onsite[0]-left)
    for i in range(1,ns):
        gl[i] = np.linalg.inv(z*I-onsite[i]-hop[i-1].conj().T @ gl[i-1] @ hop[i-1])
    GiN = np.empty_like(gr)
    GiN[-1] = np.linalg.inv(z*I-onsite[-1]-right-hop[-1].conj().T @ gl[-2] @ hop[-1])
    for i in range(ns-2,-1,-1):
        GiN[i] = gl[i] @ hop[i] @ GiN[i+1]

    g0N = GiN[0]
    T = float(np.sum((np.abs(g0N)**2) * gamma_l[:,None] * gamma_r[None,:]))
    if not with_density:
        return T, None, None
    aL = (Gi0 * gamma_l[None,None,:]) @ np.swapaxes(Gi0.conj(),1,2)
    aR = (GiN * gamma_r[None,None,:]) @ np.swapaxes(GiN.conj(),1,2)
    return T, aL, aR


def evaluate_negf(cfg, phi_3d, vg, vd, density=True):
    """Ballistic, mode-coupled NEGF per valley with equilibrium reservoirs.

    mu_source is input relative to reference E_c, mu_drain=mu_source-Vd.
    NEGF G^< includes injection from each contact, not only |G^R|.
    """
    phi_3d = np.asarray(phi_3d, float)
    if phi_3d.shape != (cfg.nx,cfg.ny,cfg.nz):
        raise ValueError('potential grid shape must be (nx,ny,nz)')
    energies = np.linspace(cfg.energy_min_ev,cfg.energy_max_ev,cfg.energy_points)
    de = energies[1]-energies[0]
    weights = np.full(len(energies),de)
    weights[[0,-1]] *= .5
    kt = KB_EV*cfg.temperature_k
    muL = cfg.efermi_source_ev
    muR = muL-vd
    fL = fermi(energies,muL,kt)
    fR = fermi(energies,muR,kt)
    total_current = 0.0
    total_n = np.zeros((cfg.nx,cfg.ny,cfg.nz),float)
    total_T = np.zeros(len(energies),float)
    valley_currents = {}
    dx = cfg.length_nm/(cfg.nx-1)
    dy = cfg.width_nm/(cfg.ny+1)
    dz = cfg.thickness_nm/(cfg.nz+1)
    for valley in cfg.valleys:
        e, u, t, onsite, hop = build_mode_slices(cfg,valley,phi_3d)
        rho = np.zeros((cfg.nx,cfg.modes,cfg.modes),complex)
        valley_T = np.zeros(len(energies),float)
        for k, energy in enumerate(energies):
            T,aL,aR = spectral_rgf(onsite,hop,t,e[0],e[-1],energy,cfg.eta_ev,
                                  with_density=density)
            valley_T[k] = T
            if density:
                rho += (weights[k]/(2*np.pi)) * (fL[k]*aL + fR[k]*aR)
        # spin factor 2, explicit valley multiplicity
        cur = 2 * valley.multiplicity * Q*Q/H_PLANCK * np.sum(weights*valley_T*(fL-fR))
        valley_currents[valley.name] = cur
        total_current += cur
        total_T += valley_T*valley.multiplicity
        if density:
            # density at each discrete (x,y,z) site from local density matrices.
            cell_nm3 = dx*dy*dz
            for i in range(cfg.nx):
                ni = np.real(np.einsum('am,mn,an->a',u[i],rho[i],u[i].conj()))
                total_n[i] += (2*valley.multiplicity*ni.reshape(cfg.ny,cfg.nz)
                               / cell_nm3 * 1e21)  # nm^-3 -> cm^-3
    return QuantumResults(float(total_current),total_n,energies,total_T,valley_currents)
