#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import LabelFrame, Label, Entry, Button, END, StringVar
from bisect import bisect_right


class TimeFrame:

    def __init__(self, parent, gui, sim):
        self.timeframe = LabelFrame(parent, text="Time", labelanchor="ne")
        self.timeframe.grid(row=0, column=1, sticky="nswe")

        """ start time of the GUI """
        self.starttimelabel = Label(self.timeframe, text="Start Time:")
        self.starttimelabel.grid(row=0, column=0, sticky="nswe", padx=10, pady=10)

        self.start_time_entry = Entry(self.timeframe)
        self.start_time_entry.grid(row=0, column=1, sticky="nswe", padx=10, pady=10)
        self.start_time_entry.bind('<Return>', lambda event, sim=sim, gui=gui,
                                                      self=self: self.start_time_update(event, self,
                                                                                        gui, sim))

        self.starttimebutton = Button(self.timeframe, text="First",
                                      command=lambda: self.set_start_time(parent, gui, sim))
        self.starttimebutton.grid(row=0, column=2, sticky="nswe", padx=10, pady=10)

        self.var_start = StringVar()
        self.var_start.set('#')
        self.start_step = Label(self.timeframe, textvariable=self.var_start)
        self.start_step.grid(row=0, column=3, sticky="nswe", padx=10, pady=10)

        """ end time of the GUI """
        self.endtimelabel = Label(self.timeframe, text="End Time:")
        self.endtimelabel.grid(row=1, column=0, sticky="nswe", padx=10, pady=10)

        self.end_time_entry = Entry(self.timeframe)
        self.end_time_entry.grid(row=1, column=1, sticky="nswe", padx=10, pady=10)
        self.end_time_entry.bind('<Return>',
                                 lambda event, sim=sim, gui=gui, self=self: self.end_time_update(
                                     event, self, gui, sim))

        self.endtimebutton = Button(self.timeframe, text="Last",
                                    command=lambda: self.set_end_time(parent, gui, sim))
        self.endtimebutton.grid(row=1, column=2, sticky="nswe", padx=10, pady=10)

        self.var_stop = StringVar()
        self.var_stop.set('#')
        self.end_step = Label(self.timeframe, textvariable=self.var_stop)
        self.end_step.grid(row=1, column=3, sticky="nswe", padx=10, pady=10)

        """ stepping of the diagnostic """
        self.var_stepping = StringVar()
        self.var_stepping.set(str(sim.stepping))
        self.stepping_label = Label(self.timeframe, text="Stepping:")
        self.stepping_label.grid(row=2, column=0, sticky="nswe", padx=10, pady=10)
        self.stepping_entry = Entry(self.timeframe, width=10, textvariable=self.var_stepping)
        self.stepping_entry.grid(row=2, column=1, sticky="nswe", padx=10, pady=10)
        self.stepping_entry.bind('<Return>',
                                 lambda event, sim=sim, self=self: self.stepping_update(event, self,
                                                                                        gui, sim))

        """ max steps of the diagnostic """
        self.var_steps = StringVar()
        self.var_steps.set(str(sim.max_steps))
        self.steps_label = Label(self.timeframe, text="Max_steps:")
        self.steps_label.grid(row=2, column=2, sticky="nswe", padx=10, pady=10)
        self.steps_entry = Entry(self.timeframe, width=10, textvariable=self.var_steps)
        self.steps_entry.grid(row=2, column=3, sticky="nswe", padx=10, pady=10)
        self.steps_entry.bind('<Return>',
                              lambda event, sim=sim, self=self: self.steps_update(event, self, gui,
                                                                                  sim))

    def steps_update(self, event, parent, gui, sim):
        # try:
        steps = parent.steps_entry.get()
        if steps == 'all':
            self.var_steps.set(steps)
            sim.max_steps = steps
            gui.status.info_txt.insert(END, "Max steps set to " + sim.max_steps + "\n")
            gui.status.info_txt.see(END)
        else:
            if int(steps):
                sim.max_steps = steps
                gui.status.info_txt.insert(END, "Max steps set to " + sim.max_steps + "\n")
                gui.status.info_txt.see(END)

    # except  ValueError:
    #     pass

    def stepping_update(self, event, parent, gui, sim):
        # try:
        stepping = int(parent.stepping_entry.get())
        self.var_stepping.set(str(stepping))
        sim.stepping = stepping
        gui.status.info_txt.insert(END, "Stepping set to " + str(sim.stepping) + " steps \n")
        gui.status.info_txt.see(END)

    #  except  ValueError:
    #      pass

    def set_start_time(self, parent, gui, sim):
        whole_t = self.find_tab_times(gui, sim)
        if len(whole_t) > 0:
            act_start = min(whole_t)
            self.start_time_entry.delete(0, END)
            self.start_time_entry.insert(END, act_start)
            self.var_start.set(str(0))
            sim.starttime = act_start

    def set_end_time(self, parent, gui, sim):
        whole_t = self.find_tab_times(gui, sim)
        if len(whole_t) > 0:
            act_stop = max(whole_t)
            self.end_time_entry.delete(0, END)
            self.end_time_entry.insert(END, act_stop)
            self.var_stop.set(str(len(whole_t)))
            sim.endtime = act_stop

    def end_time_update(self, event, parent, gui, sim):
        whole_t = self.find_tab_times(gui, sim)
        # try:
        act_endtime = float(parent.end_time_entry.get())
        pos = self.__calc_nearest__(whole_t, act_endtime)
        self.end_time_entry.delete(0, END)
        self.end_time_entry.insert(END, whole_t[pos])
        self.var_stop.set(str(pos))
        sim.endtime = whole_t[pos]

    #       except  ValueError:
    #           pass

    def start_time_update(self, event, parent, gui, sim):
        whole_t = self.find_tab_times(gui, sim)
        # try:
        act_starttime = float(parent.start_time_entry.get())
        pos = self.__calc_nearest__(whole_t, act_starttime)
        self.start_time_entry.delete(0, END)
        self.start_time_entry.insert(END, whole_t[pos])
        self.var_start.set(str(pos))
        sim.starttime = whole_t[pos]

    #   except  ValueError:
    #       pass

    def find_tab_times(self, gui, sim):
        """ this function gets all the time for the specific tab being selected
            based on what set in the diagnostic frame"""

        act_tab = gui.diags.notebook.tab(gui.diags.notebook.select(), "text")

        if sim.extensions:
            return getattr(getattr(sim.data.avail_times, getattr(gui.diags, act_tab + "_fld")),
                           'times')
        else:
            gui.status.info_txt.insert(END, "Select a simulation first \n")
            gui.status.info_txt.see(END)
            return []

    def __calc_nearest__(self, timearray, target):
        """ For a single time, find element closest to given input time"""
        if target >= timearray[-1]:
            return len(timearray) - 1
        elif target <= timearray[0]:
            return 0
        else:
            return bisect_right(timearray, target)
