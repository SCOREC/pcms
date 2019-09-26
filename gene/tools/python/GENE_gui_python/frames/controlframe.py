#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import LabelFrame, Button, END
from utils.loader import Loader


class ControlFrame:

    def __init__(self, parent, gui, sim):

        self.controlframe = LabelFrame(parent, bd=0)
        self.controlframe.grid(row=2, column=1, sticky="nswe")

        self.start_btn = Button(self.controlframe, text="Start",
                                command=lambda: self.fire(gui, sim))
        self.start_btn.grid(row=0, column=3, sticky="nswe", padx=10, pady=10)

        self.stop_btn = Button(self.controlframe, text="Stop")
        self.stop_btn.grid(row=0, column=4, sticky="nswe", padx=10, pady=10)
        self.controlframe.grid_columnconfigure(0, weight=1, uniform="group2")
        self.controlframe.grid_columnconfigure(1, weight=1, uniform="group2")
        self.controlframe.grid_columnconfigure(2, weight=1, uniform="group2")
        self.controlframe.grid_columnconfigure(3, weight=1, uniform="group2")

    def fire(self, gui, sim):
        gui.status.info_txt.insert(END, "*******************\n")
        gui.status.info_txt.insert(END, "Starting postprocessing, buckle up!\n")
        if len(gui.diag_active) == 0:
            gui.status.info_txt.insert(END, " No diagnostic selected  \n")
            gui.status.info_txt.see(END)
        else:
            loader = Loader(gui.diag_active, sim.run[0], sim.data)
            loader.set_interval(sim.data, sim.starttime, sim.endtime)
            gui.status.info_txt.see(END)

            for it, time in enumerate(loader.processable_times):
                gui.status.info_txt.insert(END, " time {}\n".format(time))
                gui.status.info_txt.see(END)
                for i_d, diag in enumerate(gui.diag_active):
                    diag.execute(sim.data, loader.processable[it])

        gui.status.info_txt.insert(END, "\nVisualizing results\n")
        gui.status.info_txt.see(END)
        for i_d, diag in enumerate(gui.diag_active):
            diag.plot(loader.processable_times, gui.status)
