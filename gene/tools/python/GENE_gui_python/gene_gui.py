#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 17:22:27 2018
@author: GENE DEVELOPMENT TEAM
"""

import tkinter as tk
from tkinter import Frame

from frames.menubar import Menubar

from frames.simulationframe import SimulationFrame
from frames.diagframe import DiagFrame
from frames.toolsframe import ToolsFrame
from frames.statusframe import StatusFrame
from frames.timeframe import TimeFrame
from frames.controlframe import ControlFrame

from utils.run import Simulation

LARGE_FONT = ("Verdana", 12)
NORM_FONT = ("Verdana", 10)
SMALL_FONT = ("Verdana", 8)


class GENE_gui(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "GENE Diagnostic Tool 1.0")

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        self.minsize(width=640, height=360)
        self.maxsize(width=self.winfo_screenheight(), height=self.winfo_screenwidth())

        self.menubar = Menubar(self)

        self.sim = Simulation()
        self.diag_list = []
        self.diag_active = []

        self.frames = self.Create_Frames()

    def Create_Frames(self):
        self.top_frame = Frame(master=self)
        self.top_frame.grid(row=0, column=0, sticky="nswe")

        self.bottom_frame = Frame(master=self)
        self.bottom_frame.grid(row=1, column=0, sticky="nswe")

        self.grid_rowconfigure(0, weight=1, uniform="group1")
        self.grid_rowconfigure(1, weight=3, uniform="group1")
        self.grid_columnconfigure(0, weight=1)

        """ top half """
        self.toolsframe = ToolsFrame(self.top_frame, self, self.sim)
        self.simframe = SimulationFrame(self.top_frame, self, self.sim)

        self.control = ControlFrame(self.bottom_frame, self, self.sim)
        self.top_frame.grid_columnconfigure(1, weight=1, uniform="group1")

        """ bottom half"""
        self.status = StatusFrame(self.bottom_frame)
        self.diags = DiagFrame(self.bottom_frame, self)
        self.times = TimeFrame(self.bottom_frame, self, self.sim)

        self.bottom_frame.grid_columnconfigure(0, weight=3, uniform="group1")
        self.bottom_frame.grid_columnconfigure(1, weight=2, uniform="group1")
        self.bottom_frame.grid_rowconfigure(0, weight=2, uniform="group2")
        self.bottom_frame.grid_rowconfigure(1, weight=3, uniform="group2")
        self.bottom_frame.grid_rowconfigure(2, weight=1, uniform="group2")

        self.bottom_frame.grid_rowconfigure(0, weight=3, uniform="group2")
        self.bottom_frame.grid_rowconfigure(1, weight=5, uniform="group2")
        self.bottom_frame.grid_rowconfigure(2, weight=1, uniform="group2")

    def _close(self):
        self.quit()


if __name__ == '__main__':
    app = GENE_gui()
    app.resizable(width=True, height=True)
    app.geometry("1280x720")
    app.mainloop()
