#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import LabelFrame, Text, Scrollbar


class StatusFrame:

    def __init__(self, parent):
        self.infoframe = LabelFrame(parent, text="Status", labelanchor="ne")
        self.infoframe.grid(row=1, column=1, sticky="nsew")

        self.info_txt = Text(self.infoframe)
        self.info_txt.grid(row=0, column=0, sticky="nswe", padx=10, pady=10)
        self.scrollb = Scrollbar(self.infoframe, command=self.info_txt.yview)
        self.scrollb.grid(row=0, column=1, sticky='nsew')
        self.info_txt['yscrollcommand'] = self.scrollb.set

        self.infoframe.grid_columnconfigure(0, weight=1, uniform="group1")
        self.infoframe.grid_rowconfigure(0, weight=1, uniform="group1")
