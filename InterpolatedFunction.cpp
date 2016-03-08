/*
 * InterpolatedFunction.cpp
 *
 *  Created on: 4 Sep 2010
 *      Author: daniel
 */

#include "InterpolatedFunction.h"
#include <fstream>
#include "debug.h"
#include <iostream>
#include "LHAPDF/LHAPDF.h"

#include "pdf.h"

using namespace std;

void InterpolatedFunction::setup(const std::vector<double>& x,const std::vector<double>& fx){
	d_acc = gsl_interp_accel_alloc ();
	d_spline = gsl_spline_alloc (gsl_interp_cspline, x.size());;
	gsl_spline_init (d_spline, &x[0], &fx[0], x.size());
}




QofAlphasInterpolated::QofAlphasInterpolated(){
	double q,alpha;
	std::vector<double> xs; std::vector<double> fxs;
	

	
#ifdef LHAPDF_NEW_VERSION
	
	d_min=currentPDF::s_PDF->alphasQ(14000);
	d_max=currentPDF::s_PDF->alphasQ(1);
#else
	d_min=LHAPDF::alphasPDF(14000);
	d_max=LHAPDF::alphasPDF(1);
#endif



	
	for (int i=14000;i>10000;i-=500){
		q=i;
#ifdef LHAPDF_NEW_VERSION
		alpha=currentPDF::s_PDF->alphasQ(q);
#else
		alpha=LHAPDF::alphasPDF(q);
#endif
		xs.push_back(alpha);
		fxs.push_back(q);
	}
	for (int i=10000;i>2000;i-=20){
		q=i;
#ifdef LHAPDF_NEW_VERSION
		alpha=currentPDF::s_PDF->alphasQ(q);
#else
		alpha=LHAPDF::alphasPDF(q);
#endif
		xs.push_back(alpha);
		fxs.push_back(q);
	}
	for (int i=2000;i>90;i-=5){
		q=i;
#ifdef LHAPDF_NEW_VERSION
		alpha=currentPDF::s_PDF->alphasQ(q);
#else
		alpha=LHAPDF::alphasPDF(q);
#endif
		xs.push_back(alpha);
		fxs.push_back(q);
	}
	for (int i=90;i>=1;i--){
		q=i;
#ifdef LHAPDF_NEW_VERSION
		alpha=currentPDF::s_PDF->alphasQ(q);
#else
		alpha=LHAPDF::alphasPDF(q);
#endif
		xs.push_back(alpha);
		fxs.push_back(q);
	}
	setup(xs,fxs);
}

double QofAlphasInterpolated::operator()(double x){
	if (x<d_min){
		return InterpolatedFunction::operator ()(d_min);
	}
	if (x>d_max){
		return InterpolatedFunction::operator ()(d_max);
	}
	return InterpolatedFunction::operator ()(x);

}

