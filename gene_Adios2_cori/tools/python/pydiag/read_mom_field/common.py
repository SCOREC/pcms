import numpy as np
from os.path import join


class common_data:
    def __init__(self, p):
        global pars
        pars = p
        self.nx = int(pars.pardict['nx0'])
        self.ny = int(pars.pardict['nky0'])
        self.nz = int(pars.pardict['nz0'])
        self.kymin = float(pars.pardict['kymin'])
        self.kygrid = np.array([self.kymin*j for j in range(self.ny)], dtype=np.float64)
        try:
            lx = float(pars.pardict['lx'])
        except KeyError:
            lx = 1
        self.zgrid = np.array([-np.pi + k*2*np.pi/self.nz for k in range(self.nz)],
                              dtype=np.float64)
        try:
            self.x_local = pars.pardict["x_local"]
        except KeyError:
            self.x_local = True
        if self.x_local:
            self.kxmin = 2*np.pi/lx
            self.kxgridp = [self.kxmin*i for i in range(int(self.nx/2 + 1))]
            self.kxgridn = [self.kxmin*i for i in range(int(-self.nx/2 + 1), 0)]
            self.kxgrid = np.array(self.kxgridp + self.kxgridn, dtype=np.float64)
            # ordered kx array
            self.kxgrid_ord = np.array(self.kxgridn + self.kxgridp[0:-1], dtype=np.float64)
        else:
            self.xgrid = np.array([-lx/2. + lx/(self.nx - 1)*i for i in range(0, self.nx)])
            try:
                self.minor_r = float(pars.pardict['minor_r'])
            except KeyError:
                pass
            self.rhostar = float(pars.pardict['rhostar'])
            self.x0 = float(pars.pardict['x0'])
            self.x_agrid = self.xgrid*self.rhostar + self.x0
        # initialize number of timesteps for different diags
        self.tlen = 0
        self.nt_en = 0
        self.nt_mom = 0

    def set_profiles(self, spec):
        try:
            dat = np.genfromtxt(join(self.rundir, 'profiles_' + spec + self.fileext))
            xgrid = dat[:, 0]
            self.temp_prof = dat[:, 2]
            self.dens_prof = dat[:, 3]
        except OSError:
            self.temp_prof = 1.
            self.dens_prof = 1.
