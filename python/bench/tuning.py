#!/usr/bin/python 

# ------------------------------------------------------------------
#  BSOF/I - Block structured factorization and inversion codes 
#  for CPU+GPU platforms
#  Copyright (c) 2013, Sergiy Gogolenko
#  e-mail: sgogolenko@ucdavis.edu
# ------------------------------------------------------------------
#  Description:
#   Classes for analysis results of performance tuning benchmarks.
# ------------------------------------------------------------------

class BSOFIBenchValues(object):
    def __init__(self, times, measures):
        super(BSOFIBenchValues, self).__init__()
        self.N, self.T = times, measures

    @property
    def NT(self): return self.N, self.T

    # def __getitem__(self, t):
    #     if k in self.keys():
    #         return self[k]
    #     raise AttributeError

class BSOFIBenchParameter(object):
    def __init__(self, times, measures):
        super(BSOFIBenchParameter, self).__init__()
        self._data = BSOFIBenchValues(times, measures)
        self._model   = None
        self._computedModel = False 

    def computeModel(self): 
        self._computedModel = True

    @property
    def data(self):
        return self._data

    @property
    def model(self):
        if not self._computedModel: self.computeModel(); 
        return self._model
        
class BSOFIDefaultBenchParameter(BSOFIBenchParameter):
    def __init__(self, times, measures):
        super(BSOFIDefaultBenchParameter, self).__init__(times, measures)

    def computeModel(self): 
        super(BSOFIDefaultBenchParameter, self).computeModel()
        self._model= self._data

class BSOFIFittedBenchParameter(BSOFIBenchParameter):
    def __init__(self, times, measures, func = None, p0 = (0.1, 1., 0.1, 200)):
        super(BSOFIFittedBenchParameter, self).__init__(times, measures)
        self.params0  = p0
        self._params  = None 
        self.fitFunction = func
        
    @property
    def params(self):
        if not self._computedModel: self.computeModel(); 
        return self._params

    def __str__(self): 
        a, b, c, d = self.params
        d   = int(d)
        f_d = self.fitFunction(d + 1, a, b, c, d)
        func_format = r"((__n) < {d} ? ({a:.4e} + {b:.4e}/((__n) - ({c:.4e}))) : {f_d:.4e})"
        return func_format.format(a = a, b = b, c = c, d = int(d), f_d = f_d)

    def computeModel(self): 
        super(BSOFIFittedBenchParameter, self).computeModel()
        from scipy.optimize import curve_fit
        self._params, self._pcov = curve_fit(self.fitFunction, self._data.N, self._data.T, p0 = self.params0)
        self._model = BSOFIBenchValues(self._data.N, self.fitFunction(self.data.N, *self._params))

import numpy 
class BSOFISteppedBenchParameter(BSOFIBenchParameter):
    def __init__(self, times, measures, 
                 approx = numpy.round, window = numpy.hanning,):
        super(BSOFISteppedBenchParameter, self).__init__(times, measures)
        self.window  = window
        self.approx  = approx
        self._smooth = None

    @property
    def smooth(self):
        if not self._computedModel: self.computeModel(); 
        return self._smooth

    def __str__(self): 
        N, T = self.model.NT
        # ifPart = lambda n : {True  : r"(n > {n} ? {t} : {part})".\
        #                          format(part = ifPart(n-2), n = int(N[n-2]), t=int(T[n-2])), 
        #                      False : ""}[n>2]
        if_format = r"((__n) > {n} ? {t} : {part})"
        def ifPart(n):
            if n>2:
                return if_format.format(part = ifPart(n-2), n = int(N[n-2]), t=int(T[n-2]))
            else: return "{t}".format(t=int(T[n-2]))
        return ifPart(len(N))

    def computeModel(self): 
        super(BSOFISteppedBenchParameter, self).computeModel()
        self._model  = None
        self._smooth = None

        def smoothener(x, win_length, window = self.window):
            # Windows: 'hanning', 'hamming', 'bartlett', 'blackman'
            s=numpy.r_[x[win_length-1:0:-1],x,x[-1:-win_length:-1]]
            w = window(win_length)
            return numpy.convolve(w/w.sum(), s, mode='valid')
    
        for wl in (2**numpy.array(range(5,10))):
            xs, ys = \
                smoothener(self._data.N, win_length = wl), \
                smoothener(self._data.T, win_length = wl)
            imax, imin = numpy.argmax(xs), numpy.argmin(xs)
            xs, ys = xs[imax:imin], ys[imax:imin]
            ym     = self.approx(ys)
            y_diff = numpy.diff(ym)
            if min(y_diff) == 0 or max(y_diff) == 0:
                xc, yc = [self._data.N[0]], [ym[0]]
                for i in xrange(len(ym)): 
                    if(ym[i] > yc[-1]):
                        xc.append(xs[i]); xc.append(xs[i]); 
                        yc.append(yc[-1]); yc.append(ym[i])
                xc.append(self._data.N[-1]); yc.append(yc[-1])
                self._win_length = wl
                self._smooth = BSOFIBenchValues(xs[::-1], ys[::-1])
                self._model  = BSOFIBenchValues(numpy.array(xc[::-1]), numpy.array(yc[::-1]))
                return
        return

class BSOFIPerformanceModel(dict):
    def __init__(self, report, measure = 'time', n0 = None, n_end = None, n_step = 1):#,*arg,**kw):
        if n0 is None: n0 = 0
        else: n0 = n0
        if n_end is None: n_end = len(report['n'])
        else: n_end = n_end
        N = report['n'][n0:n_end:n_step]

        #eta_F = 0.5; eta_I = 1. - eta_F
        eta_I, eta_F = 0.5, 0.5

        funcFitKapa = lambda x, a,b,c,d: a + b/(numpy.minimum(x,d)-c)

        time = { 'time' : lambda name: report.time(name)[n0:n_end:n_step],
                 'perf' : lambda name: report.instr(name)[n0:n_end:n_step]/report.time(name)[n0:n_end:n_step]/1e3,
                 'flops': lambda name: report.instr(name)[n0:n_end:n_step]}[measure]
        kappa_C = time('GPU_GEMM')/time('CPU_GEMM')
        kappa_R = time('GPU_GEMM')/time('CPU_GEMM')
        kappa_Q = time('GPU_GEMM2')/time('CPU_GEMM2')
        alpha_C = (time('CPU_TRTRI')+2*time('CPU_TRMML'))/(2*time('CPU_GEMM'))
        alpha_R = (time('CPU_TRTRI')+3*time('CPU_TRMML'))/(2*time('CPU_GEMM'))
        beta_1  = (time('CPU_ORGQR'))/(5*time('CPU_GEMM'))
        beta_2  = abs(2*time('GPU_GEMM') - time('GPU_GEMM2'))
        alpha_F = (time('CPU_TRTRI2'))/(2*time('CPU_GEMM'))
        beta_F  = (time('CPU_GEQRF') + time('CPU_ORMQR'))/(5*time('CPU_GEMM'))
        l_min   = (1 - 2*alpha_C) + 5*beta_F*eta_F/eta_I

        super(BSOFIPerformanceModel, self).__init__({
                'kappa_C' : BSOFIFittedBenchParameter(N, kappa_C, func = funcFitKapa),
                'kappa_R' : BSOFIFittedBenchParameter(N, kappa_R, func = funcFitKapa),
                'kappa_Q' : BSOFIFittedBenchParameter(N, kappa_Q, func = funcFitKapa),
                'alpha_C' : BSOFIDefaultBenchParameter(N, alpha_C),
                'alpha_R' : BSOFIDefaultBenchParameter(N, alpha_R),
                'beta_1'  : BSOFIDefaultBenchParameter(N, beta_1),
                'beta_2'  : BSOFIDefaultBenchParameter(N, beta_2),
                'alpha_F' : BSOFIDefaultBenchParameter(N, alpha_F),
                'beta_F'  : BSOFIDefaultBenchParameter(N, beta_F),
                'l_min'   : BSOFIDefaultBenchParameter(N, l_min),
                })

        self['c_i']  = BSOFISteppedBenchParameter(N, 2*(alpha_R - 1.))
        self['c_j']  = BSOFISteppedBenchParameter(N, 2*(alpha_C - 1.))
        self['c_k1'] = BSOFISteppedBenchParameter(N, (5*beta_1/2 - 2 - kappa_Q*beta_2))
        self['c_k2'] = BSOFISteppedBenchParameter(N, ((5*beta_1/2 - 2)/2 - kappa_Q*beta_2))
        
        self['l_F']  = BSOFISteppedBenchParameter(N, l_min - eta_I/(5.*beta_F*eta_F)* \
                                                      ((l_min*(l_min+1))/2 - 1 - \
                                                           2*(l_min-2)*(1-alpha_C) - 2*(1 - alpha_F)), \
                                                      approx = numpy.ceil)

        pass

    def __getattr__(self, k):
        if k in self.keys():
            return self[k]
        raise AttributeError

    def __str__(self):
        strModel = "\n"
        for param in ('kappa_R', 'kappa_Q', 'l_F', 
                      r'c_i', r'c_j', r'c_k1', r'c_k2',):
            strModel += "#define {param_name}(__n) {param_expr}\n".format(
                param_name = param.upper(), 
                param_expr = str(self[param]))
        return strModel
