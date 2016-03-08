/*
 * InterpolatedFunction.h
 *
 *  Created on: 4 Sep 2010
 *      Author: daniel
 */

#ifndef INTERPOLATEDFUNCTION_H_
#define INTERPOLATEDFUNCTION_H_

#include <gsl/gsl_spline.h>
#include <vector>
#include <string>
#include <exception>
#include "LHAPDF/LHAPDF.h"

#ifdef LHAPDF_MAJOR_VERSION
#if LHAPDF_MAJOR_VERSION == 6
#define LHAPDF_NEW_VERSION
#endif
#endif


class InterpolatedFunctionNoFile : public std::exception {
  LHAPDF::PDF* d_PDF;
public:
	virtual const char* what() const throw()
	  {
	    return "File not found in InterpolatedFunction";
	  }

};

class InterpolatedFunction {
	gsl_interp_accel*	d_acc;
	gsl_spline* d_spline;
public:
	InterpolatedFunction(){};
	virtual double operator()(double x){return gsl_spline_eval (d_spline, x, d_acc);}
	virtual ~InterpolatedFunction(){}
protected:
	void setup(const std::vector<double>& x,const std::vector<double>& fx);
};



class QofAlphasInterpolated: public InterpolatedFunction {
	double d_max;
	double d_min;

public:
	QofAlphasInterpolated();
	virtual double operator()(double x);
};



#endif /* INTERPOLATEDFUNCTION_H_ */
