"""DEVSIM 3D linear Poisson with user-supplied NEGF electron concentration.

All Poisson lengths in cm, densities cm^-3, electrostatic potential in volts.
No DD equations are used. Source/drain/gate are Dirichlet; remaining exterior
boundaries have natural zero normal displacement flux.
"""
import numpy as np
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator

Q = 1.602176634e-19
EPS0_F_PER_CM = 8.8541878128e-14


class DevsimPoisson:
    def __init__(self, cfg, mesh_path):
        import devsim as d
        self.d = d
        self.cfg = cfg
        self.device = 'mono_gaa'
        self.si = 'silicon'
        self.ox = 'oxide'
        mesh = 'mono_gaa_mesh'
        d.create_gmsh_mesh(mesh=mesh, file=str(mesh_path))
        d.add_gmsh_region(mesh=mesh, gmsh_name='Silicon', region=self.si, material='Silicon')
        d.add_gmsh_region(mesh=mesh, gmsh_name='GateOxide', region=self.ox, material='Oxide')
        d.add_gmsh_contact(mesh=mesh, gmsh_name='Source', region=self.si, name='source', material='Metal')
        d.add_gmsh_contact(mesh=mesh, gmsh_name='Drain', region=self.si, name='drain', material='Metal')
        d.add_gmsh_contact(mesh=mesh, gmsh_name='Gate', region=self.ox, name='gate', material='Metal')
        d.add_gmsh_interface(mesh=mesh, gmsh_name='SiOx_Interface', region0=self.si,
                             region1=self.ox, name='si_ox')
        d.finalize_mesh(mesh=mesh)
        d.create_device(mesh=mesh, device=self.device)
        print('DEVSIM regions/contacts/interface:', d.get_region_list(device=self.device),
              d.get_contact_list(device=self.device), d.get_interface_list(device=self.device))
        for reg, er in ((self.si, cfg.eps_si), (self.ox, cfg.eps_ox)):
            d.set_parameter(device=self.device, region=reg, name='Permittivity',
                            value=er * EPS0_F_PER_CM)
            d.node_solution(device=self.device, region=reg, name='Potential')
            d.edge_from_node_model(device=self.device, region=reg, node_model='Potential')
            d.edge_model(device=self.device, region=reg, name='ElectricField',
                         equation='(Potential@n0-Potential@n1)*EdgeInverseLength')
            d.edge_model(device=self.device, region=reg, name='ElectricField:Potential@n0',
                         equation='EdgeInverseLength')
            d.edge_model(device=self.device, region=reg, name='ElectricField:Potential@n1',
                         equation='-EdgeInverseLength')
            d.edge_model(device=self.device, region=reg, name='DField',
                         equation='Permittivity*ElectricField')
            d.edge_model(device=self.device, region=reg, name='DField:Potential@n0',
                         equation='Permittivity*EdgeInverseLength')
            d.edge_model(device=self.device, region=reg, name='DField:Potential@n1',
                         equation='-Permittivity*EdgeInverseLength')
            if reg == self.si:
                # Input is set externally, NOT an additional unknown to solve.
                d.node_solution(device=self.device, region=reg, name='InjectedElectrons')
                d.node_model(device=self.device, region=reg, name='SpaceCharge',
                             equation='-ElectronCharge*(NetDoping-InjectedElectrons)')
                d.node_model(device=self.device, region=reg,
                             name='SpaceCharge:Potential', equation='0')
                d.set_parameter(device=self.device, region=reg, name='ElectronCharge', value=Q)
                d.set_parameter(device=self.device, region=reg, name='NdSourceDrain',
                                value=cfg.donor_sd_cm3)
                d.set_parameter(device=self.device, region=reg, name='NdChannel',
                                value=cfg.donor_channel_cm3)
                d.set_parameter(device=self.device, region=reg, name='LSource',
                                value=cfg.ls_nm*1e-7)
                d.set_parameter(device=self.device, region=reg, name='LChannelEnd',
                                value=(cfg.ls_nm + cfg.lg_nm)*1e-7)
                d.node_model(device=self.device, region=reg, name='NetDoping',
                             equation='ifelse(x < LSource, NdSourceDrain, '
                                      'ifelse(x > LChannelEnd, NdSourceDrain, NdChannel))')
                d.equation(device=self.device, region=reg, name='PotentialEquation',
                           variable_name='Potential', edge_model='DField',
                           node_model='SpaceCharge', variable_update='default')
            else:
                d.equation(device=self.device, region=reg, name='PotentialEquation',
                           variable_name='Potential', edge_model='DField',
                           variable_update='default')
        d.interface_model(device=self.device, interface='si_ox', name='PhiContinuity',
                          equation='Potential@r0-Potential@r1')
        d.interface_model(device=self.device, interface='si_ox', name='PhiContinuity:Potential@r0',
                          equation='1')
        d.interface_model(device=self.device, interface='si_ox', name='PhiContinuity:Potential@r1',
                          equation='-1')
        d.interface_equation(device=self.device, interface='si_ox', name='PotentialEquation',
                             interface_model='PhiContinuity', type='continuous')
        for contact in ('source', 'drain', 'gate'):
            d.set_parameter(device=self.device, name=f'{contact}_bias', value=0.0)
            d.contact_node_model(device=self.device, contact=contact,
                                 name=f'{contact}_bc', equation=f'Potential-{contact}_bias')
            d.contact_node_model(device=self.device, contact=contact,
                                 name=f'{contact}_bc:Potential', equation='1')
            d.contact_equation(device=self.device, contact=contact, name='PotentialEquation',
                               node_model=f'{contact}_bc', edge_charge_model='DField')
        self.si_nodes_nm = self.get_coords(self.si)
        self.ox_nodes_nm = self.get_coords(self.ox)
        self.set_charge(np.zeros(len(self.si_nodes_nm)))
        self._si_interp = None

    def get_coords(self, region):
        d = self.d
        return np.stack([np.asarray(d.get_node_model_values(device=self.device,
                         region=region, name=axis), dtype=float) * 1e7
                         for axis in ('x','y','z')], axis=1)

    def set_charge(self, n_cm3):
        n_cm3 = np.asarray(n_cm3, float)
        if n_cm3.shape != (len(self.si_nodes_nm),):
            raise ValueError('Silicon node electron density length mismatch')
        if not np.all(np.isfinite(n_cm3)) or n_cm3.min() < -1e-6:
            raise ValueError('Invalid electron density from NEGF')
        self.d.set_node_values(device=self.device, region=self.si,
                               name='InjectedElectrons', values=np.maximum(0, n_cm3).tolist())

    def set_charge_from_regular_grid(self, x, y, z, n):
        interp = RegularGridInterpolator((x, y, z), n, bounds_error=False,
                                         fill_value=0, method='linear')
        self.set_charge(interp(self.si_nodes_nm))

    def solve(self, vg, vd):
        d = self.d
        for c, bias in (('source',0), ('drain',vd), ('gate',vg)):
            d.set_parameter(device=self.device, name=f'{c}_bias', value=float(bias))
        d.solve(type='dc', absolute_error=1e-14, relative_error=1e-9,
                maximum_iterations=20, solver_type='direct')
        return np.asarray(d.get_node_model_values(device=self.device, region=self.si,
                                                  name='Potential'), float)

    def sample_si_potential(self, x, y, z):
        phi = np.asarray(self.d.get_node_model_values(device=self.device, region=self.si,
                                                      name='Potential'), dtype=float)
        if self._si_interp is None:
            from scipy.spatial import Delaunay
            # triangulation used repeatedly over iterations; avoid rebuilding it
            self._si_interp = Delaunay(self.si_nodes_nm)
        interpolator = LinearNDInterpolator(self._si_interp, phi, fill_value=np.nan)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        vals = interpolator(np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])).reshape(X.shape)
        if np.any(~np.isfinite(vals)):
            raise RuntimeError('Poisson sampling outside silicon mesh: check mesh and transport grid')
        return vals

    def export_vtk(self, path):
        self.d.write_devices(file=str(path), type='vtk')
