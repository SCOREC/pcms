#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from abc import abstractmethod
from tkinter import Label, END, Frame, Text, Entry, OptionMenu, StringVar, BooleanVar
from tkinter import Tk
import itertools


class Diagnostic(object):
    def __init__(self):
        """ When we initialize each diagnostic class this call the set buttons  """
        self.name = 'my_name'
        self.tab = 'my_tab'
        pass

    #    """" This method defines the default options for a specific diagnostic"""
    ##    TODO move this in init and get rid of this method
    #    @abstractmethod
    #    def default_options(self):
    #
    #        pass

    """ This method must be called before starting the diagnstotic loop. 
        Will set the optiopns for the specific run and return info to the loader."""

    @abstractmethod
    def set_options(self):
        pass

    """ The diagnostic itself"""

    @abstractmethod
    def execute(self):
        pass

    """ for diagnostics requring time averages, execute saves time dependent data
        in a relative (h5) file. For non time-average we just keep data in memory.
        Plot method plots the results"""

    @abstractmethod
    def plot(self):
        pass

    """ we set the default options, in case the option gui is not called"""

    def set_defaults(self):
        for opt_id, opt_info in self.opts.items():
            if not opt_info['values']:
                """option has no default, e.g. might be different for local/global"""
                self.opts[opt_id]['values'] = None
                self.opts[opt_id]['value'] = None
                self.opts[opt_id]['bol'] = False
            elif type(opt_info['values'][0]) == type(True):
                """option is a boolean"""
                self.opts[opt_id]['value'] = self.opts[opt_id]['values'][0]
                self.opts[opt_id]['bol'] = True
            else:
                """option is not a boolean"""
                self.opts[opt_id]['value'] = self.opts[opt_id]['values'][0]
                self.opts[opt_id]['bol'] = False

    def ListAvailVars(self, filtered):
        """ filters the available data accoridng to request"""
        if self.avail_vars:
            return list(itertools.chain.from_iterable(
                    [self.avail_vars[x].values() for x in self.avail_vars.keys() if x in filtered]))
        else:
            return None

    def ListSpecnames(self):
        """ returns a list of all species plus all if n_spec>1"""
        if self.specnames:
            return ['all'] + self.specnames if len(self.specnames) > 1 else self.specnames
        else:
            return None

    def Get_NeededVars(self, filtered=None):
        self.needed_vars = {}
        """ if we give a list of file names, set corresponding entries to true"""
        if filtered:
            for x in filtered:
                self.needed_vars[x] = {y: True for y in self.avail_vars[x].values()}
        else:
            """ use this to loop over quants from options"""
            for file_key in self.avail_vars.keys():
                self.needed_vars[file_key] = {}
                for data_value in self.avail_vars[file_key].values():
                    if (data_value in self.opts['quant']['value']):
                        self.needed_vars[file_key][data_value] = True
                    else:
                        self.needed_vars[file_key][data_value] = False

    def Get_Spec_From_Opts(self):
        spec_tmp = self.opts['spec']['value']
        self.specnames = [spec_tmp] if spec_tmp != 'all' else self.run_data.specnames

    def __reload__(self, avail_vars, specnames):
        """ when we select a fdiagnostic or we change the file keeping it selcted,
            we need to reload variables and species"""

        self.avail_vars = avail_vars
        self.specnames = specnames
        self.set_defaults()

    class Options_GUI(Tk):
        """ This manages the additional gui to set options or a give diagnostic"""

        def __init__(self):
            """ Nothing to do"""
            pass

        def startup(self, diagnostic):
            """ When we press the option button this creates a new GUI.
                Not sure this is ok or we need the toplevel instead
                The option gui ais a window broken into two pieces, left has help,
                right has options"""
            #                TODO check toplevel
            #                TODO block gui allowing only one options at the time
            Tk.__init__(self)
            Tk.wm_title(self, diagnostic.name + " options")
            print("Updating options for " + diagnostic.name)
            self.columnconfigure(0, weight=1)
            self.rowconfigure(0, weight=1)

            self.minsize(width=50, height=35)
            self.maxsize(width=self.winfo_screenheight(), height=self.winfo_screenwidth())

            self.helpframe = Frame(self)
            self.helpframe.grid(row=0, column=0, sticky="nswe")
            self.help_txt = Text(self.helpframe)
            self.help_txt.grid(row=0, column=0, sticky="nswe", padx=10, pady=10)
            self.help_txt.grid_columnconfigure(0, weight=1)
            self.help_txt.grid_rowconfigure(0, weight=1)
            self.help_txt.insert(END, diagnostic.help_txt)

            self.optframe = Frame(self)
            self.optframe.grid(row=0, column=1, sticky="nswe")

            self.helpframe.grid_columnconfigure(0, weight=1, uniform="group2")
            self.optframe.grid_columnconfigure(1, weight=1, uniform="group2")

            """ each time we creaate the gui we keep track of all options that 
                the specific diagnostic has"""
            self.opts_input = []
            for i_opt, opt_name in enumerate(diagnostic.opts.keys()):
                self.opts_input.append(
                    Diagnostic.Option(self.optframe, opt_name, i_opt, diagnostic))

    class Option(object):
        """ Here we crat the specific single option.
            Need to distinguish between paramters that can be set by the user or predefined ones"""

        def __init__(self, parent, opt_name, idx, diag):
            self.variable = StringVar(parent)  # BooleanVar(parent) if diag.opts[opt_name]['bol']
            # else StringVar(parent)
            self.label = Label(parent, text=diag.opts[opt_name]['tag'])
            self.label.grid(row=idx, column=0, sticky="nswe", padx=10, pady=10)

            if not diag.opts[opt_name]['values']:
                self.entry = Entry(parent, width=10, textvariable=self.variable)
                self.entry.grid(row=idx, column=1, sticky="nswe", padx=10, pady=10)
                self.entry.bind('<Return>',
                                lambda event, diag=diag, opt_name=opt_name: self.opt_updt(event,
                                                                                          diag,
                                                                                          opt_name))
            else:
                self.variable.set(
                    str(diag.opts[opt_name]['value']) if diag.opts[opt_name]['bol'] else
                    diag.opts[opt_name]['value'])
                self.menu = OptionMenu(parent, self.variable,
                                       *[str(x) for x in diag.opts[opt_name]['values']] if
                                       diag.opts[opt_name]['bol'] else diag.opts[opt_name][
                                           'values'], command=lambda event, diag=diag,
                                                                     opt_name=opt_name:
                    self.opt_updt(
                        event, diag, opt_name))
                self.menu.grid(row=idx, column=1, sticky="nswe", padx=10, pady=10)

        def opt_updt(self, event=None, diag=None, opt_name=None):
            """ Updates the option value"""
            value = self.variable.get()
            if diag.opts[opt_name]['bol']:
                diag.opts[opt_name]['value'] = value == 'True'
            else:
                diag.opts[opt_name]['value'] = value
            print("Updated " + opt_name + " to ", diag.opts[opt_name]['value'])
