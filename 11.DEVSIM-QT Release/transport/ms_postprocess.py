"""Semi-infinite scattering current after SCF; no density returned to MS SCF."""
import numpy as np
from config import Q,HBAR_EV_S,KB_EV
from multispace.occupation import fermi
from transport.negf_reference import device_blocks,lead_matrices,rgf
from transport.quadrature import integrate_vector
from transport.resonances import resonance_breaks


def transmission(system,c,energy,vd):
    ref=c.potential_reference_v
    ons,hop,tx,leads=device_blocks(system,c,(ref,ref+vd));e=np.atleast_1d(energy).astype(complex)
    values=[]
    for start in range(0,len(e),c.transmission_batch_size):
        ee=e[start:start+c.transmission_batch_size]
        sl,sr=lead_matrices(ee,tx,leads);_,amp=rgf(ee,ons,hop,sl,sr,False,0)
        gam=1j*(sr-sr.conj().transpose(0,2,1));a=amp[-1]
        values.append(np.einsum('eij,ejk,eik->e',gam,a,a.conj(),optimize=True).real)
    return np.concatenate(values)


def postprocess(system,c,mu,vd):
    ref=c.potential_reference_v
    ons,hop,tx,leads=device_blocks(system,c,(ref,ref+vd));size=len(hop)
    if abs(vd)<1e-14:return {'current_a':0.,'bonds_a':np.zeros(size),'nodes':0,'conservation_relative':0.}
    lo=min(mu,mu-vd)-c.fermi_tail*KB_EV*c.temperature
    hi=max(mu,mu-vd)+c.fermi_tail*KB_EV*c.temperature
    breaks,info=resonance_breaks(ons,hop,tx,leads,lo,hi,c.energy_step,c.energy_order)
    def integrand(e):
        sl,sr=lead_matrices(e.astype(complex),tx,leads);_,a=rgf(e,ons,hop,sl,sr,False,0)
        gamma=1j*(sr-sr.conj().transpose(0,2,1));v=a[-1]
        T=np.einsum('eij,ejk,eik->e',gamma,v,v.conj(),optimize=True).real
        flux=[-2*np.einsum('ab,ebk,eak->e',t,a[i+1],a[i].conj(),optimize=True).imag for i,t in enumerate(hop)]
        df=(fermi(e,mu,c.temperature)-fermi(e,mu-vd,c.temperature))/(2*np.pi)
        return np.column_stack([T,*flux])*df[:,None]
    integral,error=integrate_vector(integrand,breaks,c.energy_order,c.energy_rtol,2e-15,c.energy_max_nodes)
    g=2*system.valley.multiplicity;curr=integral.real*Q/HBAR_EV_S*g
    return {'current_a':float(curr[0]),'bonds_a':curr[1:],'nodes':error['gf_nodes'],
            'conservation_relative':float(max(abs(curr[1:]-curr[0]))/max(abs(curr[0]),1e-15)),
            'quadrature':error,'pole_guidance':info,'gf_inside_electronic_scf':False}


def direct_transmission(system,c,energy,vd):
    """Independent dense inversion, for small models and sparse spectrum checks."""
    ref=c.potential_reference_v
    ons,hop,tx,leads=device_blocks(system,c,(ref,ref+vd));off=np.r_[0,np.cumsum([len(x) for x in ons])]
    H=np.zeros((off[-1],off[-1]),complex)
    for i,h in enumerate(ons):H[off[i]:off[i+1],off[i]:off[i+1]]=h
    for i,t in enumerate(hop):
        H[off[i]:off[i+1],off[i+1]:off[i+2]]=t;H[off[i+1]:off[i+2],off[i]:off[i+1]]=t.conj().T
    ans=[]
    for e in np.atleast_1d(energy):
        sl,sr=lead_matrices(np.array([complex(e)]),tx,leads);s=np.zeros_like(H)
        s[:len(ons[0]),:len(ons[0])]+=sl[0];s[-len(ons[-1]):,-len(ons[-1]):]+=sr[0]
        G=np.linalg.solve(e*np.eye(len(H))-H-s,np.eye(len(H)))[:len(ons[0]),-len(ons[-1]):]
        gl=1j*(sl[0]-sl[0].conj().T);gr=1j*(sr[0]-sr[0].conj().T)
        ans.append(np.trace(gl@G@gr@G.conj().T).real)
    return np.asarray(ans)
