import unittest
import os
import numpy as np
from collections import OrderedDict

from pydiag.utils.ParIO import Parameters


class TestParameters(unittest.TestCase):
    def setUp(self):
        self.par = Parameters()
        self.par.Read_Pars("parameters.dat")
        self.pnt = self.par.asnamedtuple()

    def testRead_Pars(self):
        reference = OrderedDict(
                [('n_procs_s', 1), ('n_procs_v', 1), ('n_procs_w', 1), ('n_procs_x', 1),
                 ('n_procs_y', 1), ('n_procs_z', 2), ('n_procs_sim', 2), ('n_spec', 1), ('nx0', 24),
                 ('nky0', 8), ('nz0', 16), ('nv0', 32), ('nw0', 8), ('kymin', 0.4707), ('lv', 3.0),
                 ('lw', 9.0), ('lx', 160.0), ('x0', 0.5), ('n0_global', 30), ('diagdir', "'.//'"),
                 ('read_checkpoint', False), ('write_checkpoint', True), ('istep_field', 20),
                 ('istep_mom', 20), ('istep_nrg', 20), ('istep_vsp', 50), ('istep_schpt', 5000),
                 ('istep_energy', 50), ('istep_prof', 20), ('istep_srcmom', 40),
                 ('istep_fsa_moments', 40), ('write_std', True), ('nonlinear', True),
                 ('x_local', False), ('comp_type', "'IV'"), ('perf_vec', '1 1 2 2 1 1 1 1 1'),
                 ('nblocks', 64), ('timescheme', "'RK4'"), ('dt_max', 0.00538),
                 ('dt_vlasov', 0.00538), ('ev_coll', 0.021793), ('courant', 1.25),
                 ('timelim', 1000), ('ntimesteps', 100), ('underflow_limit', -1.0), ('beta', 0.0),
                 ('debye2', 0.0), ('collision_op', "'landau'"), ('coll', 0.0001),
                 ('coll_cons_model', "'xu_rosenbluth'"), ('init_cond', "'db'"), ('hyp_z', 0.2),
                 ('hyp_v', 0.2), ('l_buffer_size', 0.025), ('lcoef_krook', 1.0),
                 ('u_buffer_size', 0.025), ('ucoef_krook', 1.0), ('explicit_buffer', True),
                 ('rad_bc_type', 1), ('magn_geometry', "'circular'"), ('trpeps', 0.18),
                 ('major_R', 1.0), ('minor_r', 0.36),
                 ('q_coeffs', '0.85000000,   0.0000000,   2.2000000'), ('mag_prof', True),
                 ('rhostar', 0.0056), ('dpdx_term', "'gradB_eq_curv'"), ('dpdx_pm', 0.0),
                 ('norm_flux_projection', False), ('name1', "'ions'"), ('prof_type1', 2),
                 ('kappa_T1', 6.96), ('LT_center1', 0.5), ('LT_width1', 0.3), ('kappa_n1', 2.23),
                 ('Ln_center1', 0.5), ('Ln_width1', 0.3), ('mass1', 1.0), ('temp1', 1.0),
                 ('dens1', 1.0), ('charge1', 1), ('step_time', 0.4095),
                 ('number of computed time steps', 100), ('time for initial value solver', 40.95),
                 ('calc_dt', True), ('nltdt_off', False), ('ev_coll_est', 0.033719517),
                 ('init_time', 24.7879), ('n_fields', 1), ('n_moms', 6), ('nrgcols', 10),
                 ('ly', 13.3492), ('PRECISION', 'DOUBLE'), ('ENDIANNESS', 'LITTLE'),
                 ('OMP_NUM_THREADS', 1), ('GIT_BRANCH', '4db53734ac0718005d04494da24a986923f8d76c'),
                 ('GIT_MASTER', '6982ecdba20201e3443742016b7e9fafd4474e37'),
                 ('RELEASE', '1.8 - alpha 0'), ('nu_ei', 0.02424), ('nustar_i', 0.00276),
                 ('Tref', 1), ('ky0_ind',0), ('kx_center',0.0),
                 ('omega_prec', 1E-3),
                 ('Omega0_tor', 0.0), ('ExBrate', 0.0),
                 ('with_coriolis', False), ('with_centrifugal', False),
                 ('with_comoving_other', False), ('with_bxphi0', False),
                 ('sign_Ip_CW', 1),('sign_Bt_CW', 1),('n_pol',1),
                 ('nref', 1.0), ('Bref', 1.0), ('mref', 1.999),
                 ("Lref", 1.0)])

        nmlref = OrderedDict([('n_procs_s', 'parallelization'), ('n_procs_v', 'parallelization'),
                              ('n_procs_w', 'parallelization'), ('n_procs_x', 'parallelization'),
                              ('n_procs_y', 'parallelization'), ('n_procs_z', 'parallelization'),
                              ('n_procs_sim', 'parallelization'), ('n_spec', 'box'), ('nx0', 'box'),
                              ('nky0', 'box'), ('nz0', 'box'), ('nv0', 'box'), ('nw0', 'box'),
                              ('kymin', 'box'), ('lv', 'box'), ('lw', 'box'), ('lx', 'info'),
                              ('x0', 'box'), ('n0_global', 'box'), ('diagdir', 'in_out'),
                              ('read_checkpoint', 'in_out'), ('write_checkpoint', 'in_out'),
                              ('istep_field', 'in_out'), ('istep_mom', 'in_out'),
                              ('istep_nrg', 'in_out'), ('istep_vsp', 'in_out'),
                              ('istep_schpt', 'in_out'), ('istep_energy', 'in_out'),
                              ('istep_prof', 'in_out'), ('istep_srcmom', 'in_out'),
                              ('istep_fsa_moments', 'in_out'), ('write_std', 'in_out'),
                              ('nonlinear', 'general'), ('x_local', 'general'),
                              ('comp_type', 'general'), ('perf_vec', 'general'),
                              ('nblocks', 'general'), ('timescheme', 'general'),
                              ('dt_max', 'general'), ('dt_vlasov', 'general'),
                              ('ev_coll', 'general'), ('courant', 'general'),
                              ('timelim', 'general'), ('ntimesteps', 'general'),
                              ('underflow_limit', 'general'), ('beta', 'general'),
                              ('debye2', 'general'), ('collision_op', 'general'),
                              ('coll', 'general'), ('coll_cons_model', 'general'),
                              ('init_cond', 'general'), ('hyp_z', 'general'), ('hyp_v', 'general'),
                              ('l_buffer_size', 'nonlocal_x'), ('lcoef_krook', 'nonlocal_x'),
                              ('u_buffer_size', 'nonlocal_x'), ('ucoef_krook', 'nonlocal_x'),
                              ('explicit_buffer', 'nonlocal_x'), ('rad_bc_type', 'nonlocal_x'),
                              ('magn_geometry', 'geometry'), ('trpeps', 'geometry'),
                              ('major_R', 'geometry'), ('minor_r', 'geometry'),
                              ('q_coeffs', 'geometry'), ('mag_prof', 'geometry'),
                              ('rhostar', 'geometry'), ('dpdx_term', 'geometry'),
                              ('dpdx_pm', 'geometry'), ('norm_flux_projection', 'geometry'),
                              ('name1', 'species1'), ('prof_type1', 'species1'),
                              ('kappa_T1', 'species1'), ('LT_center1', 'species1'),
                              ('LT_width1', 'species1'), ('kappa_n1', 'species1'),
                              ('Ln_center1', 'species1'), ('Ln_width1', 'species1'),
                              ('mass1', 'species1'), ('temp1', 'species1'), ('dens1', 'species1'),
                              ('charge1', 'species1'), ('step_time', 'info'),
                              ('number of computed time steps', 'info'),
                              ('time for initial value solver', 'info'), ('calc_dt', 'info'),
                              ('nltdt_off', 'info'), ('ev_coll_est', 'info'), ('init_time', 'info'),
                              ('n_fields', 'info'), ('n_moms', 'info'), ('nrgcols', 'info'),
                              ('ly', 'info'), ('PRECISION', 'info'), ('ENDIANNESS', 'info'),
                              ('OMP_NUM_THREADS', 'info'), ('GIT_BRANCH', 'info'),
                              ('GIT_MASTER', 'info'), ('RELEASE', 'info'), ('nu_ei', 'info'),
                              ('nustar_i', 'info'),
                              ('Tref', 'units'), ('ky0_ind', 'box'), ('kx_center', 'box'),
                              ('omega_prec', 'general'),
                              ('Omega0_tor', 'external_contr'),('ExBrate', 'external_contr'),
                              ('with_coriolis', 'external_contr'), ('with_centrifugal', 'external_contr'),
                              ('with_comoving_other', 'external_contr'), ('with_bxphi0', 'external_contr'),
                              ('sign_Ip_CW', 'geometry'), ('sign_Bt_CW', 'geometry'), ('n_pol', 'geometry'),
                              ('nref', 'units'), ('Bref', 'units'), ('mref', 'units'),
                              ("Lref", 'units')])

        self.assertEqual(self.par.pardict, reference)
        self.assertEqual(self.par.nmldict, nmlref)

    def testWrite_Pars(self):
        self.par.Write_Pars("parameters_testwrite")
        testpar = Parameters()
        testpar.Read_Pars("parameters_testwrite")
        self.assertTrue(dict(self.par.pardict) == dict(testpar.pardict))
        self.assertEqual(dict(self.par.nmldict), dict(testpar.nmldict))
        os.remove("parameters_testwrite")

    def test_nokey(self):
        with np.testing.assert_raises(KeyError):
            self.par.pardict["thisisatest"]

    def test_notuplekey(self):
        with np.testing.assert_raises(AttributeError):
            self.pnt.thisisatest

    def test_tuplekey(self):
        self.assertEquals(self.pnt.timescheme, "'RK4'")
