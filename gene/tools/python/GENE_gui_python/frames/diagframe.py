#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import Frame
from tkinter.ttk import Notebook as Notebook

from tkinter import Label, Button, Checkbutton, END, BooleanVar
import os

""" class containing the notebook with all tabs for the different diags

    Naming of the frames  inside the notebook MUST be consistent
"""


class DiagFrame:
    def __init__(self, parent, grandparent):
        self.notebook = Notebook(parent)
        self.notebook.grid(row=0, column=0, rowspan=3, sticky="nsew", padx=10, pady=10)

        self.Basicframe = Frame(self.notebook)
        self.notebook.add(self.Basicframe, text="Basic")
        self.Basic_diags = 0
        self.Basic_fld = "field"

        self.Localframe = Frame(self.notebook)
        self.notebook.add(self.Localframe, text="Local")
        self.Local_diags = 0
        self.Local_fld = "field"

        self.Globalframe = Frame(self.notebook)
        self.notebook.add(self.Globalframe, text="Global")
        self.Global_diags = 0
        self.Global_fld = "nrg"

        self.GENE3Dframe = Frame(self.notebook)
        self.notebook.add(self.GENE3Dframe, text="GENE3D")
        self.GENE3D_diags = 0
        self.GENE3D_fld = "nrg"

        self.nrgframe = Frame(self.notebook)
        self.notebook.add(self.nrgframe, text="nrg")
        self.nrg_diags = 0
        self.nrg_fld = "nrg"

        self.vspframe = Frame(self.notebook)
        self.notebook.add(self.vspframe, text="vsp")
        self.vsp_diags = 0
        self.vsp_fld = "vsp"

        self.diag_list = []
        self.autoadd_diagnostic(grandparent)

    """ Parses the diagnostic folder and adds each of the present files to the GUI"""

    def autoadd_diagnostic(self, parent):
        """ so, we need to know if the run has been set.
            If we know the run we are rpocessing we know the list of variables
            that are available. if not we dont know it so the add will not be able to show things
        """

        for filename in os.listdir('diagnostics'):
            if filename[-3:] == '.py' and filename[0:-3] != 'diagnostic' and \
                    filename[0:-3] != 'baseplot' and filename[0:-3] != 'diagspace':
                exec('from diagnostics.' + filename[0:-3] + ' import ' + filename[0:-3])
                exec('parent.diag_list.append(Diag_Frame_Entry(' + filename[0:-3] + '(parent.sim.data.avail_vars,parent.sim.specnames), self, parent))')


class Diag_Frame_Entry(object):
    def __init__(self, fname, diagframe, GUI):

        self.diag = fname
        self.my_tab = self.diag.tab
        self.active = BooleanVar()
        self.active.initialize(False)
        self.my_row = eval('diagframe.' + self.my_tab + '_diags')
        self.my_frame = eval('diagframe.' + self.my_tab + 'frame')

        self.check = Checkbutton(self.my_frame, variable=self.active,
                                 command=lambda var_in=self.active, GUI=GUI: self.selector(var_in,
                                                                                           GUI))
        self.check.grid(row=self.my_row, column=0, sticky="nsew", padx=10, pady=10)

        self.label = Label(self.my_frame, text=self.diag.name)
        self.label.grid(row=self.my_row, column=1, sticky="nswe", padx=10, pady=10)

        self.button = Button(self.my_frame, text="Options",
                             command=lambda: self.diag.options_gui.startup(self.diag))
        self.button.grid(row=self.my_row, column=2, sticky="nswe", padx=10, pady=10)

        try:
            exec('diagframe.' + self.my_tab + '_diags=diagframe.' + self.my_tab + '_diags+1')
        except:
            print("failure :(")

    def selector(self, switch, GUI):
        """ this method activates or deactivates a given diagnostic based on the
        value of the switch"""
        if switch.get():
            """activate"""
            """ we need to add the diagnostic to the list of active diags 
                and set its default options"""
            self.diag.__reload__(GUI.sim.data.avail_vars, GUI.sim.run[0].specnames)
            GUI.diag_active.append(self.diag)
        else:
            """ deactivate """
            GUI.diag_active.remove(self.diag)
