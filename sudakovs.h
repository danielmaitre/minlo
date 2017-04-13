/*
 * sudakovs.h
 *
 *  Created on: 30 Mar 2017
 *      Author: daniel
 */

#ifndef SUDAKOVS_H_
#define SUDAKOVS_H_

#include "LHAPDF/LHAPDF.h"

double SherpaIntegrand(double q2,double Q2,char fl,LHAPDF::PDF* pdf,bool m_fo=false,int m_mode=3,double m_mur2=1.0);
double KeithIntegrand(double q2,double Q2,char fl);
double K(const double &nf,bool m_fo,int m_mode);

double KeithExponentQuark(double q02,double q2);
double KeithExponentGluon(double q02,double q2);

double SherpaExponentQuark(double q02, double q2,LHAPDF::PDF* pdf,int mode);
double SherpaExponentGluon(double q02, double q2,LHAPDF::PDF* pdf,int mode);

double SherpaSudakov(double q20,double q2h,double q2l,int flav,LHAPDF::PDF* pdf ,int mode);
double KeithSudakov(double q20,double q2h,double q2l,int flav,LHAPDF::PDF* pdf);

double alphasKeith(double Q2,int nf=5,double lambda2=0.226*0.226);
double alphasLOKeith(double Q2,int nf=5,double lambda2=0.226*0.226);



#endif /* SUDAKOVS_H_ */
