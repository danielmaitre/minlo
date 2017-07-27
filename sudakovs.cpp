/*
 * sudakov.cpp
 *
 *  Created on: 30 Mar 2017
 *      Author: daniel
 */
#include "LHAPDF/LHAPDF.h"
#include "sudakovs.h"
#include "gsl/gsl_integration.h"
#include "MinloKeithDeclarations.h"
#include "pdf.h"

const double CA=3.0;
const double CF=4.0/3.0;

/* taken from Sherpa MINLO_KFactor_Setter.C */

inline double sqr(double x){ return (x*x);}

double K(const double &nf,bool m_fo,int m_mode)
{
  if (m_fo || !(m_mode&2)) return 0.0;
  return 3.*(67./18.-sqr(M_PI)/6.)-10./9.*nf/2.;
}

double Ggq(const double &e,const double &q2,double m=0.0)
{
  if (m*m>q2) return 0.0;
  return 0.5/(q2+m*m)*sqr(1.0-e)*(1.0-(1.0-e)*(1.0+3.0*e)/3.0*q2/(q2+m*m));
}

double SherpaIntegrand(double q2,double Q2,char fl,LHAPDF::PDF* pdf,bool m_fo,int m_mode,double m_mur2,int nfgs){
	double eps(sqrt(q2/Q2));
	double e((m_mode&1)?eps:0.0);
	double as2pi=pdf->alphasQ2((m_fo?m_mur2:q2))/(2.0*M_PI);
	int nf_pdf=pdf->alphaS().numFlavorsQ2(m_fo?m_mur2:q2);
	double nf(nfgs<=nf_pdf?nfgs:nf_pdf);
	if (fl=='q') {
	  double gam=as2pi/q2*4.0/3.0*
		(2.0*log(1.0/eps)*(1.0+as2pi*K(nf,m_fo,m_mode))
		 -3.0/2.0*sqr(1.0-e));
	  /*

	  if (!m_fl.IsMassive()) return gam;
	  double k=sqrt(q2)/m_fl.Mass();
	  gam+=as2pi/q2*4.0/3.0*as2pi*K(nf)*
		(-(1.0-e)*k*k/(e*e+k*k)
		 +k*atan((1-e)*k/(e+k*k))+log((e*e+k*k)/(1.0+k*k)));
	  gam+=as2pi/q2*4.0/3.0*
		((1.0-e*e)/2.0+e*sqr(1.0-e)*(e/(k*k+e*e)-(1.0-e)/(k*k+sqr(1.0-e)))
		 -k*(atan(1.0/k)+(1.0-k*k)*atan(e*k/(1.0-e+k*k)))
		 -(1.0-k*k/2.0)*(log((k*k+1.0)/(k*k+e*e))-2.0*k*atan(e/k)));

	  return gam;
	  */
	  return gam;
	}
	if (fl=='g') {
	  double gam0=Ggq(e,q2), gam=3.0*gam0;
	  for (long int i(4);i<=nfgs;++i){
		if (false /*Flavour(i).Mass()*/) {/*gam+=Ggq(e,q2,Flavour(i).Mass());*/}
		else if (nf>=i) gam+=gam0;}
	  return as2pi/q2*3.0*
		(2.0*log(1.0/eps)*(1.0+as2pi*K(nf,m_fo,m_mode))
		 -sqr(1.0-e)/6.0*(11.0-e*(2.0-3.0*e)))
		+as2pi*gam;
	}
	return 0.0;
}

double alphasKeith(double Q2,int nf,double lambda2){
	const double b0=(33.0-2*nf)/(12*M_PI);
	const double b1=(153.0-19*nf)/(24*M_PI*M_PI);
	double l=log(Q2/lambda2);
	double alphas=1.0/(b0*l)*(1.0 - b1*log(l)/(b0*b0*l));
	return alphas;
};

double alphasLOKeith(double Q2,int nf,double lambda2){
	const double b0=(33.0-2*nf)/(12*M_PI);
	double l=log(Q2/lambda2);
	double alphas=1.0/(b0*l);
	return alphas;
};


double KeithIntegrand(double q2,double Q2,char fl){
	double A1,A2,B1;
	const double nf=5.0;
	const double b0=(33.0-2*nf)/(12*M_PI);
	double K=(67.0/18.0-M_PI*M_PI/6.0)*CA-5.0/9.0*nf;
	if (fl=='g') {
		A1=CA;
		A2=CA*K;
		B1=-2*M_PI*b0;
	} else {
		A1=CF;
		A2=CF*K;
		B1=-3.0/2.0*CF;
	}
	double a=alphasKeith(q2, nf)/(2*M_PI);
	return 1.0/q2 * ( (a*A1 +a*a*A2)*log(Q2/q2) +a*B1);
}


double KeithQuarkIntegrand(double q2,void *pQ2){
	return  KeithIntegrand(q2,*(double*)pQ2,'q');
}
double KeithGluonIntegrand(double q2,void *pQ2){
	return  KeithIntegrand(q2,*(double*)pQ2,'g');
}

struct SherpaParams {
	double Q2;
	LHAPDF::PDF* pdf;
	int mode;
	bool fo;
	int nfgs;
	double q2ren;
};

double SherpaIntegrand(double q2,void *params,char fl){
	SherpaParams* sp=(SherpaParams*)(params);
	double Q2=(*sp).Q2;
	LHAPDF::PDF* pdf =(*sp).pdf;
	int mode=(*sp).mode;

	bool fo=(*sp).fo;
	double mur2=(*sp).q2ren;
	int nfgs=(*sp).nfgs;
	return  SherpaIntegrand(q2,Q2,fl,pdf,fo,mode,mur2,nfgs);
}

double SherpaIntegrandQuark(double q2,void *params){
	return  SherpaIntegrand(q2,params,'q');
}

double SherpaIntegrandGluon(double q2,void *params){
	return  SherpaIntegrand(q2,params,'g');
}

double myIntegrate(gsl_function *F,double ql,double qh){
	int limit=5000;
	static gsl_integration_workspace *wsp=gsl_integration_workspace_alloc(limit);

	double epsabs=1e-8;
	double epsrel=1e-8;
	double result;
	double abserror;
	gsl_integration_qag (F, ql, qh, epsabs,epsrel, limit , GSL_INTEG_GAUSS21, wsp, &result, &abserror);
	return result;

}


double KeithExponentQuark(double q02, double q2){

	gsl_function F;
	F.function=&KeithQuarkIntegrand;
	F.params=&q2;

	return -myIntegrate(&F, q02, q2);

}

double KeithExponentGluon(double q02,double q2){

	gsl_function F;
	F.function=&KeithGluonIntegrand;
	F.params=&q2;

	return -myIntegrate(&F, q02, q2);

}

double SherpaExponentQuark(double q02, double q2,LHAPDF::PDF* pdf,int mode,bool fo,int nfgs,double q2ren){
	SherpaParams sp;
	sp.Q2=q2;
	sp.pdf=pdf;
	sp.mode=mode;
	sp.fo=fo;
	sp.nfgs=nfgs;
	sp.q2ren=q2ren;
	gsl_function F;
	F.function=&SherpaIntegrandQuark;
	F.params=&sp;

	return -myIntegrate(&F, q02, q2);

}

double SherpaExponentGluon(double q02, double q2,LHAPDF::PDF* pdf,int mode,bool fo,int nfgs,double q2ren){
	SherpaParams sp;
	sp.Q2=q2;
	sp.pdf=pdf;
	sp.mode=mode;
	sp.fo=fo;
	sp.nfgs=nfgs;
	sp.q2ren=q2ren;

	gsl_function F;
	F.function=&SherpaIntegrandGluon;
	F.params=&sp;

	return -myIntegrate(&F, q02, q2);

}

double SherpaSudakov(double q20,double q2h,double q2l,int flav,LHAPDF::PDF* pdf ,int mode, int ngfs){
	double num,den;
	if (flav==0){
		num=SherpaExponentGluon(q20,q2h,pdf,mode,false,ngfs,q2h /*not relevant*/);
		den=SherpaExponentGluon(q20,q2l,pdf,mode,false,ngfs,q2h /*not relevant*/);
	} else {
		num=SherpaExponentQuark(q20,q2h,pdf,mode,false,ngfs,q2h /*not relevant*/);
		den=SherpaExponentQuark(q20,q2l,pdf,mode,false,ngfs,q2h /*not relevant*/);
	}
	return exp(num-den);
}

double KeithSudakov(double q20,double q2h,double q2l,int flav,LHAPDF::PDF* pdf){
	currentPDF::setCurrent(pdf);
	return NLL_SUDAKOV(q20,q2h,q2l,flav);
}
