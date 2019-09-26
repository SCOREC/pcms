#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import filedialog, LabelFrame, Label, Entry, Button, END
import os
from pathlib import Path
import numpy as np

# TODO remove all the callbacks and make it more coincise
# Does the lambda require the inputs to be set before or not?

""" CLass for defining input, output folders and run extensions to process"""


class SimulationFrame:
    def __init__(self, parent, gui, sim):
        """ Container"""
        self.simframe = LabelFrame(parent, text="Simulation", labelanchor="ne")
        self.simframe.grid(row=0, column=1, sticky="nswe")

        """ Folder contaning the run """
        self.folderlabel = Label(self.simframe, text="Folder Path:")
        self.folderlabel.grid(row=1, column=0, sticky="nswe", padx=10, pady=10)

        self.folderentry = Entry(self.simframe)
        if Path('.oldrun.npz').exists():
            self.folderentry.insert(0, np.load('.oldrun.npz')['in_folder'])
            sim.in_folder = self.folderentry.get()
        self.folderentry.grid(row=1, column=1, columnspan=1, sticky="nswe", padx=10, pady=10)
        self.folderentry.bind('<Return>',
                              lambda event, sim=sim, self=self, gui=gui: self.folder_update(event,
                                                                                            self,
                                                                                            gui,
                                                                                            sim))

        self.folderbutton = Button(self.simframe, text=">",
                                   command=lambda gui=gui, self=self: self.set_folderentry(gui,
                                                                                           sim))
        self.folderbutton.grid(row=1, column=2, sticky="eswn", padx=10, pady=10)

        """ folder containing the output of the diagnostics """
        self.outputlabel = Label(self.simframe, text="Output Path:")
        self.outputlabel.grid(row=2, column=0, sticky="nswe", padx=10, pady=10)

        self.outputentry = Entry(self.simframe)
        if Path('.oldrun.npz').exists():
            self.outputentry.insert(0, np.load('.oldrun.npz')['out_folder'])
            sim.out_folder = self.outputentry.get()
        self.outputentry.grid(row=2, column=1, columnspan=1, sticky="nswe", padx=10, pady=10)
        self.outputentry.bind('<Return>',
                              lambda event, sim=sim, self=self, gui=gui: self.output_update(event,
                                                                                            self,
                                                                                            gui,
                                                                                            sim))

        self.outputbutton = Button(self.simframe, text=">",
                                   command=lambda: self.set_outputentry(gui, sim))
        self.outputbutton.grid(row=2, column=2, sticky="nswe", padx=10, pady=10)
        self.equalbutton = Button(self.simframe, text="=",
                                  command=lambda: self.set_equalentry(gui, sim))
        self.equalbutton.grid(row=2, column=3, sticky="nswe", padx=10, pady=10)

        """ runs to process """
        self.runslabel = Label(self.simframe, text="Runs:")
        self.runslabel.grid(row=3, column=0, sticky="nswe", padx=10, pady=10)

        self.runsentry = Entry(self.simframe)
        if Path('.oldrun.npz').exists():
            sim.extensions = [x for x in np.load('.oldrun.npz')['exts']]
            self.runsentry.insert(0, ''.join(str(elem[1:] + ',') for elem in sim.extensions))
            sim.Prepare()
        self.runsentry.grid(row=3, column=1, columnspan=1, sticky="nswe", padx=10, pady=10)
        self.runsentry.bind('<Return>',
                            lambda event, sim=sim, self=self, gui=gui: self.run_update(event, self,
                                                                                       gui, sim))

        self.runs_select_btn = Button(self.simframe, text=">",
                                      command=lambda gui=gui, sim=sim: self.set_runentry(gui, sim))
        self.runs_select_btn.grid(row=3, column=3, sticky="nswe", padx=10, pady=10)

        self.runs_autofeed_btn = Button(self.simframe, text="<", command=lambda self=self, gui=gui,
                                                                                sim=sim:
        self.autoselect_runs(
            self, gui, sim))
        self.runs_autofeed_btn.grid(row=3, column=2, sticky="nswe", padx=10, pady=10)

        self.simframe.columnconfigure(1, weight=1)

    def folder_update(self, event, parent, gui, sim):
        try:
            sim.in_folder = str(parent.folderentry.get())
            gui.status.info_txt.insert(END, "Simulation path is: " + sim.in_folder + "\n")
            gui.status.info_txt.see(END)
        except  ValueError:
            pass

    def output_update(self, event, parent, gui, sim):
        try:
            sim.out_folder = str(parent.outputentry.get())
            gui.status.info_txt.insert(END, "Output path is: " + sim.out_folder + "\n")
            gui.status.info_txt.see(END)
        except  ValueError:
            pass

    def set_folderentry(self, parent, sim):

        directory = filedialog.askdirectory()
        file_path = os.path.split(directory)
        file_path = file_path[0] + '/' + file_path[1]
        self.folderentry.delete(0, END)
        self.folderentry.insert(0, file_path)
        sim.in_folder = file_path

        parent.status.info_txt.insert(END, "Simulation path is: " + sim.in_folder + "\n")
        parent.status.info_txt.see(END)

    def set_outputentry(self, parent, sim):

        directory = filedialog.askdirectory()
        file_path = os.path.split(directory)
        file_path = file_path[0] + '/' + file_path[1]
        self.outputentry.delete(0, END)
        self.outputentry.insert(0, file_path)
        sim.out_folder = file_path
        parent.status.info_txt.insert(END, "Ouptput path is: " + sim.out_folder + "\n")
        parent.status.info_txt.see(END)

    def set_equalentry(self, parent, sim):

        self.outputentry.delete(0, END)
        self.outputentry.insert(0, sim.in_folder)
        sim.out_folder = sim.in_folder
        parent.status.info_txt.insert(END, "Output path is: " + sim.out_folder + "\n")
        parent.status.info_txt.see(END)

    def set_runentry(self, parent, sim):

        if not sim.in_folder:
            parent.status.info_txt.insert(END, "Select a simulation folder first\n")
            parent.status.info_txt.see(END)
        else:
            if not Path(sim.in_folder).is_dir():
                parent.status.info_txt.insert(END, "Simulation folder doesn't exist\n")
                parent.status.info_txt.see(END)
            else:
                #       try:
                FILEOPENOPTIONS = dict(initialdir=sim.in_folder,
                                       filetypes=[('parameters', 'parameters*')],
                                       title="Select files")
                file = filedialog.askopenfilenames(**FILEOPENOPTIONS)

                self.runsentry.delete(0, END)
                file_suffixes = []

                for item in file:
                    file_suffixes.append(os.path.split(item)[1][10:])
                    self.runsentry.insert(END, file_suffixes[file.index(item)][1:])
                    self.runsentry.insert(END, ',')

                sim.extensions = file_suffixes
                sim.Prepare()
                self.__write_history__(sim)
                parent.status.info_txt.insert(END, "Selected run: " + ''.join(
                        str(elem[1:]) for elem in sim.extensions) + "\n")
                parent.status.info_txt.see(END)

    # This exception for the moment is bad since do not allow us to debug
    #        except:
    #           parent.status.info_txt.insert(END, "You should have not reached this point\n")
    #           parent.status.info_txt.see(END)
    #           pass

    def run_update(self, event, parent, gui, sim):
        try:
            # need to be change to put the right extensions
            my_ext = str(parent.runsentry.get())
            my_ext = my_ext.split(',')
            sim.extensions = []
            for ext in my_ext:
                if Path(sim.in_folder + '/parameters_' + ext).is_file():
                    sim.extensions.append('_' + ext)
                else:
                    if Path(sim.in_folder + '/parameters.' + ext).is_file():
                        sim.extensions.append('.' + ext)

            sim.Prepare()
            self.__write_history__(sim)
            parent.runsentry.delete(0, END)
            parent.runsentry.insert(0,
                                    ''.join(str(elem[1:] + ',') for elem in sim.extensions)[0:-1])
            gui.status.info_txt.insert(END, "Selected run: " + ''.join(
                    str(elem[1:] + ' ') for elem in sim.extensions) + "\n")
            gui.status.info_txt.see(END)
        except  ValueError:
            pass

    def autoselect_runs(self, parent, gui, sim):

        if not sim.in_folder:
            gui.status.info_txt.insert(END, "Select a simulation folder first\n")
            gui.status.info_txt.see(END)
        else:
            if not Path(sim.in_folder).is_dir():
                gui.status.info_txt.insert(END, "Simulation folder doesn't exist\n")
                parent.status.info_txt.see(END)
            else:
                import glob
                sim.extensions = []
                for file in glob.glob(sim.in_folder + "/parameters*"):
                    if file[-3:] != '.h5':
                        sim.extensions.append(os.path.split(file)[1][10:])
                #                try:
                sim.Prepare()
                self.__write_history__(sim)
                parent.runsentry.delete(0, END)
                parent.runsentry.insert(0, ''.join(str(elem[1:] + ',') for elem in sim.extensions)[
                                           0:-1])
                gui.status.info_txt.insert(END, "Selected run: " + ''.join(
                        str(elem[1:] + ' ') for elem in sim.extensions) + "\n")
                gui.status.info_txt.see(END)

    #                except:
    #                    gui.status.info_txt.insert(END, "No valid files found\n")
    #                    gui.status.info_txt.see(END)
    #                    pass

    def __write_history__(self, sim):
        import numpy as np
        if sim.extensions and sim.in_folder and sim.out_folder:
            np.savez('.oldrun',
                     exts=[x[0:-3] for x in sim.extensions] if sim.run[0].is_h5 else sim.extensions,
                     in_folder=sim.in_folder, out_folder=sim.out_folder, allow_pickle=False)
