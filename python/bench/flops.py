#!/usr/bin/python

out = 'console'#'TeX'#

from sympy import init_printing, Symbol, expand, latex, symbols, simplify, oo, Function
from sympy import Rational as R
import sympy
from sympy.parsing.sympy_parser import parse_expr#, standard_transformations, implicit_application

if out=='TeX':
    headFlopsStr = r"""  \begin{tabular}{l c c c}
    \hline
    Routine & Fadds & Fmuls & Flops \\
    \hline"""
    tailFlopsStr = r"""    \hline  
  \end{tabular}"""
    formatFlopsStr = '{routine} & \n  {fadds} &\n  {fmuls} &\n  {total}\\\\'
    def expr_format(s): return "$%s$" % latex(s)
    #expr_format = latex
else:
    headFlopsStr = r""
    tailFlopsStr = r""
    formatFlopsStr = 'Flops {routine}:\n  fadds = {fadds} \n  fmuls = {fmuls} \n  flops = {total}\n'
    def expr_format(s): return s


init_printing()

# FMULS_GEQRF(__m, __n) (((__m) > (__n)) ? ((__n) * ((__n) * (  0.5-(1./3.) * (__n) + (__m)) +    (__m) + 23. / 6.)) \
#                                :                 ((__m) * ((__m) * ( -0.5-(1./3.) * (__m) + (__n)) + 2.*(__n) + 23. / 6.)) )
# FADDS_GEQRF(__m, __n) (((__m) > (__n)) ? ((__n) * ((__n) * (  0.5-(1./3.) * (__n) + (__m))            +  5. / 6.)) \
#                                :                 ((__m) * ((__m) * ( -0.5-(1./3.) * (__m) + (__n)) +    (__n) +  5. / 6.)) )
class FMULS_GEQRF(Function):
    nargs = 2
    @classmethod
    def eval(cls, __m, __n):
        return ((__n) * ((__n) * (  R(1,2)-(R(1,3)) * (__n) + __m) +    (__m) + R(23,6))) 
class FADDS_GEQRF(Function):
    nargs = 2
    @classmethod
    def eval(cls, __m, __n):
        return ((__n) * ((__n) * (  R(1,2)-(R(1,3)) * (__n) + (__m))            +  R(5,6)))
class FLOPS_DGEQRF(Function):
    nargs = 2
    @classmethod
    def eval(cls, __m, __n):
        return (FMULS_GEQRF(__m,__n) + FADDS_GEQRF(__m,__n)).simplify()

 # FMULS_UNGQR(__m, __n, __k) ((__k) * (2.* (__m) * (__n) +    2. * (__n) - 5./3. + (__k) * ( 2./3. * (__k) - ((__m) + (__n)) - 1.)))
 # FADDS_UNGQR(__m, __n, __k) ((__k) * (2.* (__m) * (__n) + (__n) - (__m) + 1./3. + (__k) * ( 2./3. * (__k) - ((__m) + (__n))     )))
class FMULS_UNGQR(Function):
    nargs = 3
    @classmethod
    def eval(cls, __m, __n,__k):
        return ((__k) * (2* (__m) * (__n) +    2 * (__n) - R(5,3) + (__k) * ( R(2,3) * (__k) - ((__m) + (__n)) - 1)))
class FADDS_UNGQR(Function):
    nargs = 3
    @classmethod
    def eval(cls, __m, __n,__k):
        return ((__k) * (2* (__m) * (__n) + (__n) - (__m) + R(1,3) + (__k) * ( R(2,3) * (__k) - ((__m) + (__n))     )))
class FLOPS_DORGQR(Function):
    nargs = 3
    @classmethod
    def eval(cls, __m, __n,__k):
        return (FMULS_UNGQR(__m,__n,__k) + FADDS_UNGQR(__m,__n,__k)).simplify()

 # FMULS_UNMQR(__m, __n, __k, __side) (( (__side) == MagmaLeft ) \
 #    ?  (2.*(__n)*(__m)*(__k) - (__n)*(__k)*(__k) + 2.*(__n)*(__k)) \
 #    :  (2.*(__n)*(__m)*(__k) - (__m)*(__k)*(__k) + (__m)*(__k) + (__n)*(__k) - 0.5*(__k)*(__k) + 0.5*(__k)))
 # FADDS_UNMQR(__m, __n, __k, __side) (( ((__side)) == MagmaLeft ) \
 #    ?  (2.*(__n)*(__m)*(__k) - (__n)*(__k)*(__k) + (__n)*(__k)) \
 #    :  (2.*(__n)*(__m)*(__k) - (__m)*(__k)*(__k) + (__m)*(__k)))

class FMULS_UNMQR(Function):
    nargs = 3
    side='L'
    @classmethod
    def eval(cls, __m, __n,__k):
        if cls.side == 'L' :  return (2*(__n)*(__m)*(__k) - (__n)*(__k)*(__k) + 2*(__n)*(__k)) 
        else:                 return (2*(__n)*(__m)*(__k) - (__m)*(__k)*(__k) + (__m)*(__k) + (__n)*(__k) - R(1,2)*(__k)*(__k) + R(1,2)*(__k))
class FADDS_UNMQR(Function):
    nargs = 3
    side='L'
    @classmethod
    def eval(cls, __m, __n,__k):
        if cls.side == 'L' :  return (2*(__n)*(__m)*(__k) - (__n)*(__k)*(__k) + (__n)*(__k))
        else:                 return (2*(__n)*(__m)*(__k) - (__m)*(__k)*(__k) + (__m)*(__k))
class FLOPS_DORMQR(Function):
    nargs = 3
    side='L'
    @classmethod
    def eval(cls, __m, __n,__k):
        return (FMULS_UNMQR(__m,__n,__k) + FADDS_UNMQR(__m,__n,__k)).simplify()

# FMULS_TRTRI(__n) ((__n) * ((__n) * ( R(1,6) * (__n) + 0.5 ) + R(1,3)))
# FADDS_TRTRI(__n) ((__n) * ((__n) * ( R(1,6) * (__n) - 0.5 ) + R(1,3)))
class FMULS_TRTRI(Function):
    nargs = 1
    @classmethod
    def eval(cls, __n):
        return ((__n) * ((__n) * ( R(1,6) * (__n) + R(1,2) ) + R(1,3)))
class FADDS_TRTRI(Function):
    nargs = 1
    @classmethod
    def eval(cls,__n):
        return ((__n) * ((__n) * ( R(1,6) * (__n) - R(1,2) ) + R(1,3)))
class FLOPS_DTRTRI(Function):
    nargs = 1
    @classmethod
    def eval(cls, __n):
        return (FMULS_TRTRI(__n) + FADDS_TRTRI(__n)).simplify()

# FMULS_GEMM(__m, __n, __k) ((__m) * (__n) * (__k))
# FADDS_GEMM(__m, __n, __k) ((__m) * (__n) * (__k))
class FMULS_GEMM(Function):
    nargs = 3
    @classmethod
    def eval(cls, __m, __n, __k):
        return ((__m) * (__n) * (__k))
class FADDS_GEMM(Function):
    nargs = 3
    @classmethod
    def eval(cls, __m, __n, __k):
        return ((__m) * (__n) * (__k))
class FLOPS_DGEMM(Function):
    nargs = 3
    @classmethod
    def eval(cls, __m, __n, __k):
        return (FMULS_GEMM(__m,__n,__k) + FADDS_GEMM(__m,__n,__k)).simplify()

# FMULS_TRMM_2(__m, __n) (0.5 * (__n) * (__m) * ((__m)+1))
# FADDS_TRMM_2(__m, __n) (0.5 * (__n) * (__m) * ((__m)-1))
# FMULS_TRMM(__side, __m, __n) ( ( (__side) == MagmaLeft ) ? FMULS_TRMM_2((__m), (__n)) : FMULS_TRMM_2((__n), (__m)) )
# FADDS_TRMM(__side, __m, __n) ( ( (__side) == MagmaLeft ) ? FADDS_TRMM_2((__m), (__n)) : FADDS_TRMM_2((__n), (__m)) )
class FMULS_TRMM(Function):
    nargs = 2
    @classmethod
    def eval(cls, __m, __n):
        return (R(1,2) * (__n) * (__m) * ((__m)+1))
class FADDS_TRMM(Function):
    nargs = 2
    @classmethod
    def eval(cls, __m, __n):
        return (R(1,2) * (__n) * (__m) * ((__m)-1))
class FLOPS_DTRMM(Function):
    nargs = 2
    @classmethod
    def eval(cls, __m, __n):
        return (FMULS_TRMM(__m,__n) + FADDS_TRMM(__m,__n)).simplify()

n, m, k = symbols(('n', 'm', 'k'))
print ((FLOPS_DORGQR(2*n,2*n,n) - 3*n**2 +4*n/3)/(5*FLOPS_DGEMM(n,n,n))).simplify()

"""
n, m, k = symbols(('n', 'm', 'k'))

print headFlopsStr
print formatFlopsStr.format(routine='DGEQRF', 
                            fadds=expr_format(FADDS_GEQRF(m,n)), 
                            fmuls=expr_format(FMULS_GEQRF(m,n)), 
                            total=expr_format(FLOPS_DGEQRF(m,n)))
print formatFlopsStr.format(routine='DORMQR', 
                            fadds=expr_format(FADDS_UNMQR(m,n,k)), 
                            fmuls=expr_format(FMULS_UNMQR(m,n,k)), 
                            total=expr_format(FLOPS_DORMQR(m,n,k)))
print formatFlopsStr.format(routine='DORGQR', 
                            fadds=expr_format(FADDS_UNGQR(m,n,k)), 
                            fmuls=expr_format(FMULS_UNGQR(m,n,k)), 
                            total=expr_format(FLOPS_DORGQR(m,n,k)))
print formatFlopsStr.format(routine='DTRTRI', 
                            fadds=expr_format(FADDS_TRTRI(n)), 
                            fmuls=expr_format(FMULS_TRTRI(n)), 
                            total=expr_format(FLOPS_DTRTRI(n)))
print formatFlopsStr.format(routine='DTRMM', 
                            fadds=expr_format(FADDS_TRMM(m,n)), 
                            fmuls=expr_format(FMULS_TRMM(m,n)), 
                            total=expr_format(FLOPS_DTRMM(m,n)))
print formatFlopsStr.format(routine='DGEMM', 
                            fadds=expr_format(FADDS_GEMM(m,n,k)), 
                            fmuls=expr_format(FMULS_GEMM(m,n,k)), 
                            total=expr_format(FLOPS_DGEMM(m,n,k)))
print tailFlopsStr

__n, __L = symbols(('n', 'p'))
# formatFlopsBSOF = r"((((__L)-2)*((2*{prefix}ORMQR(2*(__n), (__n), (__n))) \
#                               + ({prefix}GEQRF(2*(__n), (__n)))) \
#                    + ({prefix}GEQRF(2*(__n), 2*(__n)))))"
# print expr_format(simplify(parse_expr(formatFlopsBSOF.format(prefix='FLOPS_D'), 
#                                 transformations = standard_transformations + (implicit_application,))))
def fBSOF(fGEQRF,fORMQR):
    fORMQR.side='L'
    return ((((__L)-2)*((2*fORMQR(2*(__n), (__n), (__n))) + (fGEQRF(2*(__n), (__n)))) + (fGEQRF(2*(__n), 2*(__n)))))
FADDS_BSOF=simplify(fBSOF(FADDS_GEQRF,FADDS_UNMQR))
FMULS_BSOF=simplify(fBSOF(FMULS_GEQRF,FMULS_UNMQR))
FLOPS_BSOF=simplify(fBSOF(FLOPS_DGEQRF,FLOPS_DORMQR))
FLOPS_DBSOF=simplify(FADDS_BSOF + FMULS_BSOF)

def fBSTRI(fGEMM,fTRTRI,fTRMM):
    return (fTRTRI(3*(__n)) + fTRMM(__n,(__n)*(__L-3)) \
                + ((__L)-3)*(fTRTRI(__n) + 2*fTRMM(__n, __n)) \
                +(((__L)-3)*((__L)+2)/2)*(fGEMM(__n, __n, __n)))
FADDS_BSTRI=simplify(fBSTRI(FADDS_GEMM,FADDS_TRTRI,FADDS_TRMM))
FMULS_BSTRI=simplify(fBSTRI(FMULS_GEMM,FMULS_TRTRI,FMULS_TRMM))
FLOPS_DBSTRI=simplify(FADDS_BSTRI + FMULS_BSTRI)

def fBSOI(fGEQRF,fORMQR):
    fORMQR.side='R'
    return (((__L)-2)*(fORMQR((__L)*(__n), 2*(__n), (__n))) \
			      + (fORMQR((__L)*(__n), 2*(__n), 2*(__n))))
FADDS_BSOI=simplify(fBSOI(FADDS_GEQRF,FADDS_UNMQR))
FMULS_BSOI=simplify(fBSOI(FMULS_GEQRF,FMULS_UNMQR))
FLOPS_DBSOI=simplify(FADDS_BSOI + FMULS_BSOI)

print headFlopsStr
print formatFlopsStr.format(routine='BSOF', 
                            fadds=expr_format(FADDS_BSOF), 
                            fmuls=expr_format(FMULS_BSOF), 
                            total=expr_format(FLOPS_DBSOF))

print formatFlopsStr.format(routine='BSTRI', 
                            fadds=expr_format(FADDS_BSTRI), 
                            fmuls=expr_format(FMULS_BSTRI), 
                            total=expr_format(FLOPS_DBSTRI))

print formatFlopsStr.format(routine='BSOI', 
                            fadds=expr_format(FADDS_BSOI), 
                            fmuls=expr_format(FMULS_BSOI), 
                            total=expr_format(FLOPS_DBSOI))
print tailFlopsStr

FLOPS_BSOFI = simplify(FLOPS_DBSOF + FLOPS_DBSTRI + FLOPS_DBSOI)
print 'Total:', expr_format(FLOPS_BSOFI)
"""
