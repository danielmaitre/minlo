# distutils: language=c++

cimport clhapdf as c
cimport lhapdf 

from libcpp cimport bool

cdef extern from "../sudakovs.h":
     double SherpaIntegrand(double q2,double Q2,char fl,c.PDF* pdf,bool m_fo,int m_mode,double m_mur2)
     double KeithIntegrand(double q2,double Q2,char fl)
     double K(const double &nf,bool m_fo,int m_mode)

cdef extern from "../sudakovs.h":
     double KeithExponentQuark(double q02,double q2)
     double KeithExponentGluon(double q02,double q2)
     double SherpaExponentQuark(double q02, double q2,c.PDF* pdf,int mode);
     double SherpaExponentGluon(double q02, double q2,c.PDF* pdf,int mode);
     double SherpaSudakov(double q20,double q2h,double q2l,int flav,c.PDF* pdf ,int mode);
     double KeithSudakov(double q20,double q2h,double q2l,int flav,c.PDF* pdf);

cdef extern from "../nll.h":
     nll_sudakov_withInfo(const  double &q20,const double &q2h,const double &q2l,const int &flav);
cdef extern from "../nll.h":
     nll_exponent_withInfo(const  double &q20,const double &q2h,const double &q2l,const int &flav);



def SherpaIntegrandQuark(q2, Q2,lhapdf.PDF pdf,mode=3 ):
    return SherpaIntegrand(q2,Q2,'q',pdf._ptr,False,mode,1.0)

def SherpaIntegrandGluon(q2, Q2,lhapdf.PDF pdf ,mode=3):
    return SherpaIntegrand(q2,Q2,'g',pdf._ptr,False,mode,1.0)

def KeithIntegrandQuark(q2, Q2 ):
    return KeithIntegrand(q2,Q2,'q')

def KeithIntegrandGluon(q2, Q2 ):
    return KeithIntegrand(q2,Q2,'g')

def KeithExpQuark(q02,q2):
    return KeithExponentQuark(q02, q2)

def KeithExpGluon(q02,q2):
    return KeithExponentGluon(q02, q2)

def SherpaExpQuark(q02,q2,lhapdf.PDF pdf,int mode):
    return SherpaExponentQuark(q02, q2,pdf._ptr,mode)

def SherpaExpGluon(q02,q2,lhapdf.PDF pdf,int mode):
    return SherpaExponentGluon(q02, q2,pdf._ptr,mode)

def SherpaSudakovFactor(q20,q2h,q2l,int fl,lhapdf.PDF pdf,int mode):
    return SherpaSudakov(q20,q2h,q2l,fl,pdf._ptr,mode)

def KeithSudakovFactor(q20,q2h,q2l,int fl,lhapdf.PDF pdf):
    return KeithSudakov(q20,q2h,q2l,fl,pdf._ptr)

cdef extern from "../sudakovs.h":
     double alphasKeith(double Q2,int nf)

def alphasQ2Keith(Q2,nf=5):
    return  alphasKeith(Q2,nf)

def alphasQKeith(Q,nf=5):
    return  alphasKeith(Q*Q,nf)

def k(nf,mode):
    return  K(nf,False,mode)
