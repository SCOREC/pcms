import numpy as np


class averages:
    def __init__(self, cm, gm):
        global common, geom, nx, ny, nz
        common = cm
        geom = gm
        nx = common.nx
        ny = common.ny
        nz = common.nz
        if not common.x_local:
            self.jaco3d = np.repeat(geom.jacobian[:, :], ny, axis=1).reshape(nz, ny, nx)

    def init_weights(self, times):
        # init weights for time averages
        tlen = len(times)
        weights = np.empty(tlen, dtype=np.float64)
        if tlen > 1:
            for i in range(1, tlen - 1):
                weights[i] = 0.5*(times[i + 1] - times[i - 1])
            weights[0] = times[1] - times[0]
            weights[-1] = times[-1] - times[-2]
        else:
            weights[0] = 1
        return weights

    def t(self, var):
        if len(var) == common.tlen:
            return np.average(var, weights=common.tweights, axis=0)
        elif len(var) == common.nt_en:
            return np.average(var, weights=common.tweights_en, axis=0)
        elif len(var) == common.nt_mom:
            return np.average(var, weights=common.tweights_mom, axis=0)

    def z(self, var):
        if common.x_local:
            return np.average(var, weights=geom.jacobian, axis=0)
        else:
            return np.average(var, weights=self.jaco3d, axis=0)
