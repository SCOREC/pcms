#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import *
import tkinter as tk


class Menubar:

    def __init__(self, parent):
        self.menubar = tk.Menu(parent)

        file = tk.Menu(self.menubar, tearoff=0)

        file.add_command(label="Save setting")
        file.add_separator()
        file.add_command(label="Exit")

        # This all goes to File cascade
        self.menubar.add_cascade(label="File", menu=file)
        tk.Tk.config(parent, menu=self.menubar)
