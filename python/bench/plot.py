#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Protting benchmarking results.
# ------------------------------------------------------------------

from report import *
from stackplot import *

#import os, sys
try:    from matplotlib.pyplot import figure, show, savefig
except: raise

file_name_format_default    = r"./{target_platform}/mkl_{n:04d}_{threads:02d}.blk"
measure_name_format_default = r"{name}_{type}"

def main1(plot, x_val, y_val, n_list, target_platform, 
          threads = 12, y_type = 'G', xlim = None, ylim = None,
          plotting = "", legend_pos = 'upper left',
          measure_name_format = measure_name_format_default,
          file_name_format    = file_name_format_default):
    
    #reports = {'mkl_1024_06': ReportTable(file_name_format.(n = 256, threads = 12))}

    plot.grid(True)
    if xlim is not None : plot.set_xlim(*xlim)
    if ylim is not None : plot.set_ylim(*ylim)

    plot.set_ylabel({'G':'Performance, GFlop/s', 'I': '$I, 10^6$ Flops', 'T' : 'Time, ms'}[y_type])
    plot.set_xlabel({'N':r'$N=n \times P$', 'L': 'P (number of blocks)'}[x_val])
    legend_names = []
    plot_routine = {'logy' : plot.semilogy, 'logx' : plot.semilogx,}.get(plotting, plot.plot)
    for n in n_list:
        report = ReportTable(file_name_format.format(target_platform=target_platform, n = n, threads = threads))
        x, y = report[x_val], report[measure_name_format.format(name=y_val, type=y_type)]
        plot_routine(x, y, {True:'-s', False: '-'}[len(y) < 20])
        legend_names.append("$n = %s$" % n)
        del report
        pass
    plot.legend(legend_names, legend_pos, shadow = True)#, mode="expand")

def main2(plot, x_val, y_val_list, n, target_platform, 
          threads = 12, y_type = 'T', xlim = None, ylim = None,
          measure_name_format = measure_name_format_default,
          file_name_format    = file_name_format_default):
    
    plot.grid(True)

    plot.set_ylabel({'I': 'Floating-point operations', 'T' : 'Time'}[y_type])#(r'T\\Phase 1 Phase 2 Phase 3')#
    plot.set_xlabel({'N':r'$N=n \times P$', 'L': r'$P$ (number of blocks)'}[x_val])

    report = BenchReportTable(file_name_format.format(target_platform=target_platform, n = n, threads = threads))

    x, y_list = report[x_val], report.percentage(y_val_list, y_type)
    legend_names = y_val_list #(name for name in y_val_list)

    plots = []
    if len(x) > 10:
        if xlim is not None : plot.set_xlim(*xlim)
        plot.set_ylim(0,100)
        raxis = plot.twinx(); raxis.set_ylim(100, 0)
        plot.figure.canvas.draw()
        plot.set_yticklabels([r'{0}\%'.format(item.get_text()) 
                              for item in plot.get_yticklabels()])
        raxis.set_yticklabels([r'{0}\%'.format(item.get_text()) 
                               for item in raxis.get_yticklabels()])

        plots.append(stackplot(plot, x, *y_list, 
                               colors=("#999999", '#eeeeee', '#777777')))
        # y_stack = cumsum(row_stack(y_list), axis=0)
        # a = plot.fill_between(x, 0, y_stack[0,:], color="#CC6666", alpha=.7, hatch='//')
        # for i, name in zip(xrange(len(y_val_list)-1), y_val_list[1:]):
        #     plot.fill_between(x, y_stack[i], y_stack[i+1], color="#CC6666", alpha=.7)
        #     pass
    else:
        if xlim is None : x_max = max(x)
        else: x_max = xlim[1]
        x = array(x)
        width  = x_max/(len(x[x <= x_max])*2.5)
        y_stack = cumsum(row_stack(y_list), axis=0)
        for i, name in zip(range(len(y_val_list))[::-1], y_val_list[::-1]):
            plots.append(plot.bar(x, y_stack[i,:], width, 
                                  color='y', hatch='//'))#, bottom=y_list[i]
            pass
        plot.legend(plots, legend_names, 'upper left', shadow = True)

def main3(plot, x_val, y_val_frac, n_list, target_platform, 
          threads, y_type = 'T', xlim = None, ylim = None,
          plotting = "", legend_pos = 'upper left',
          measure_name_format = measure_name_format_default,
          file_name_format    = file_name_format_default):
    plot.grid(True)
    if xlim is not None : plot.set_xlim(*xlim)
    if ylim is not None : plot.set_ylim(*ylim)

    y_val_frac0, y_val_frac1 = y_val_frac
    if not isinstance(y_val_frac[0], (list,tuple)):
        y_val_frac0 = (y_val_frac0,)
        y_val_frac1 = (y_val_frac1,)
        pass
    names1 = '+'.join(y_val_frac0)
    names2 = '+'.join(y_val_frac1)
    if len(y_val_frac0) == 3: names1, names2 = 'BSOF/I', 'GEQRF/I'
    plot.set_ylabel(r'$T_{{%s}}/T_{{%s}}$' % (names1, names2))
    plot.set_xlabel({'N':r'$N=n \times P$', 'L': 'P (number of blocks)'}[x_val])
    legend_names = []
    for n in n_list:
        file_name = file_name_format.format(target_platform=target_platform, n = n, threads = threads)
        report = BenchReportTable(file_name)
        x, y = report[x_val], report.value(y_val_frac1, y_type)/report.value(y_val_frac0, y_type)
        plot.semilogy(x, y, {True:'-s', False: '-'}[len(y) < 20])
        legend_names.append("$n = %s$" % n)
        del report
        pass
    plot.legend(legend_names, 'upper left', shadow = True)
    pass


import math 
from numpy import arange
def main4(plot, N, routines, n_list = (32, 64, 100, 256, 512, 1024), 
          threads = 12, 
          plotting = "", legend_pos = 'upper left', width = 0.35, 
          file_name_format    = file_name_format_default):
    perf = {'gpu' : [], 'cpu 12':[], 'cpu  6':[],}
    if ('BSOF' not in routines and 'BSTRI' not in routines): perf['cpu/gpu'] = []

    plot.set_xlabel('$N \times P$')#r'$N=n \times P$')
    plot.set_ylabel('Performance, GFlop/s')

    ind = arange(len(n_list))
    xticklabels = [] #plot.set_xticklabels( [r"$%s \times %s$" %(n, N/n) for n in n_list])
    for n in n_list:
        L = 0
        for target in perf.keys():
            target_platform, threads = {'gpu': ('gpu/bsofi', 12), 
                                        'cpu 12': ('cpu/trsm', 12) , 
                                        'cpu  6': ('cpu/trsm', 6) , 
                                        'cpu/gpu' : ('hybrid',12)
                                        }[target]
            report = BenchReportTable(file_name_format.format(target_platform=target_platform, n = n, threads = threads))
            x = report['N']; idx = 0
            for i, valN in zip(xrange(len(x)), x):
                if math.fabs(valN - N) < n: idx = i; break
                pass
            if 'BSOFI' == routines[0] and target != 'cpu/gpu': 
                routines_list = ['BSOF', 'BSTRI']
                routines_list.extend(routines[1:])
            else: routines_list = routines
            perf[target].append(report.value(routines_list, 'I')[idx]/report.value(routines_list, 'T')[idx])
            L = int(report['L'][idx])
            del x, report
        xticklabels.append(r"$%s \times %s$" %(n, L))
    plot.set_xticks(ind+2*width)
    plot.set_xticklabels(xticklabels)
    plots = []
    plots.append(plot.bar(ind        , perf['cpu  6'], width, color='r', ))
    plots.append(plot.bar(ind+  width, perf['cpu 12'], width, color='y', ))
    plots.append(plot.bar(ind+2*width, perf['gpu'], width, color='g', ))
    if ('BSOF' not in routines and 'BSTRI' not in routines): 
        plots.append(plot.bar(ind+3*width, perf['cpu/gpu'], width, color='b', ))
    plot.legend(plots, ('CPU (6 cores)', 'CPU  (12 cores)', 
                        'GPU', 'CPU/GPU'), legend_pos)

    #plot.legend(plots, [name.upper() for name in perf.keys()], legend_pos)
    # pos = idx
    # for target in perf.keys():
    #     plot.bar(pos, perf[target], width)
    #     ind +=width

class Decorations(object):
    def __init__(self, decors = ['-']):
        self._decors = decors
        self.idx = -1
    def next(self): 
        self.idx = (self.idx + 1) % len(self._decors); return self._decors[self.idx]
