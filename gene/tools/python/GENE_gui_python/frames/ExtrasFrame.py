#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tkinter as tk


class ExtrasFrame:

    def __init__(self, parent):
        self.ExtrasFrame = tk.Frame(parent)
        self.ExtrasFrame.grid(row=3, column=1, sticky="nswe")
