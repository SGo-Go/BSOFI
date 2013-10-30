#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Parsing benchmarking results.
# ------------------------------------------------------------------

import os, sys
from numpy import array, zeros, cumsum, row_stack, multiply

class ReportTable(object):
    """ Class for parsing reports 
    """
    def __init__(self, file_name = None):
        self._params = {}
        if file_name is not None: self.load(file_name)

    def load(self, file_name, remove_comments = True):
        self.file_name = file_name

        fh = open(file_name)
        self._params = {}

        names = fh.readline().split()
        for name in names: self._params[name] = []

        for line in fh.readlines():
            try:
                try:
                    if remove_comments: line = line[0:line.index('#')] # remove comments
                    else: line = line.replace('#','')
                except ValueError: pass
                values = line.split()
                if len(values) == len(self._params):
                    for name, value in zip(names, values): 
                        self._params[name].append(float(value))
                else: pass
                del values
            except: pass
        for name in names: self._params[name] = array(self._params[name])

    def __getitem__(self, name):
        if name == 'N':
            return array([n*L for n, L in zip(self._params['n'], self._params['L'])])
        if self._params.has_key(name):
            return self._params[name]
        else:
            raise KeyError('Unknown attribute "%s"' % name)

    def sum(self, *args):
        result = array(list(self._params[args[0]]))
        if len(args) > 1:
            for name in args[1:]:
                for i in xrange(len(result)):
                    result[i] += self._params[name][i]
        return result

class BenchReportTable(ReportTable):
    """ Class for parsing benchmarking reports 
    """

    measure_name_format = r"{name}_{type}"
    
    def __init__(self, file_name = None):
        super(BenchReportTable, self).__init__(file_name)

    def value(self, names, val_type):
        if isinstance(names, (list, tuple)): 
            if len(names) == 1 and isinstance(names[0], (list, tuple)): 
                lstNames = names[0]
            else: lstNames = names
        else: lstNames = (names,)
        return array(self.sum(*(self.measure_name_format.format(name=name, type=val_type) for name in lstNames)))

    def time (self, *names): return self.value(names, val_type='T')
    def instr(self, *names): return self.value(names, val_type='I')

    def perf_scaled (self, names, instr_names):
        return [flops/T for T, flops in zip(self.time(name), self.instr(name))]

    def perf_scaled (self, names):
        return [flops/T for T, flops in zip(self.time(name), self.instr(name))]

    def percentage (self, names, val_type = 'T'):
        total = self.value(names, val_type)
        return [100*array(self.value(name, val_type))/array(total) for name in names]

#from stackplot import *
