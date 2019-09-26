#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# from data.base_file import GENE_file
from utils.loader import Loader
from diagnostics.diag_ballamp import diag_ballamp
from diagnostics.diag_contours import diag_contours
from diagnostics.diag_flux_spectra import diag_flux_spectra
from diagnostics.diag_amplitude_spectra import diag_amplitude_spectra
from diagnostics.diag_shearing_rate import diag_shearing_rate
from utils.loader import Loader
import utils.derivatives as der
from utils.ParIO import Parameters as Par
from data.data import Data
from src.run import Simulation

import numpy as np
import matplotlib.pyplot as plt

# import h5py
# import numpy as np
# file='./tests/test_local_std/3spec/field_0001.h5'
# fid = h5py.File(file, 'r')
# data=fid.get('/field/phi/0000000000').value

folder='./tests/test_local_std/3spec/';
folder='./tests/test_gene3d/';
runextensions  = ['_1'];
t_start = 10.0
t_end = 10.0


sim = Simulation(folder, folder, runextensions)
# act_run = Run(folder, runextensions[0])
# sim.run.append(act_run)
sim.Prepare()

sim.data = Data()
sim.data.load_in(sim.run[0], sim.extensions)

selected_diags = []

#selected_diags.append(diag_ballamp(sim.data.avail_vars,sim.run[0].specnames))
#selected_diags.append(diag_flux_spectra(sim.data.avail_vars,sim.run[0].specnames))
selected_diags.append(diag_contours(sim.data.avail_vars,sim.run[0].specnames))
#selected_diags.append(diag_amplitude_spectra(sim.data.avail_vars,sim.run[0].specnames))
#selected_diags.append(diag_shearing_rate(sim.data.avail_vars,sim.run[0].specnames))


loader = Loader(selected_diags, sim.run[0], sim.data)
loader.set_interval(sim.data, t_start, t_end)

for it, time in enumerate(loader.processable_times):
    print(" time {}".format(time))
    for i_d, diag in enumerate(selected_diags):
        diag.execute(sim.data, loader.processable[it])

for i_d, diag in enumerate(selected_diags):
    diag.plot(loader.processable_times)
