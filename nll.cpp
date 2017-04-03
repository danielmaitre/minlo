/*
 * nll.cpp
 *
 *  Created on: 31 Mar 2017
 *      Author: daniel
 */

#include "MinloKeithDeclarations.h"
#include "debug.h"
#include "sudakovs.h"
#include "MinloInfo.h"
#include <cmath>
#include "pdf.h"

double nll_sudakov_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav){
	NAMED_DEBUG("SUDAKOV_ARGUMENTS",std::cout <<"sudakov args q0: " << sqrt(q20)  <<" qh: " << sqrt(q2h)  <<" ql: " << sqrt(q2l)  <<" flav: " << flav <<  std::endl;)
	if (!MI.d_useSherpa){
			return NLL_SUDAKOV(q20,q2h,q2l,flav);
	} else {
		double q20l=q20;
		double q2ll=q2l;
		double q2hl=q2h;
		int fl=flav;
		return SherpaSudakov(q20, q2h, q2l, flav, currentPDF::s_PDF, MI.d_sherpaMode);
	}
};

double nll_exponent_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav){
	NAMED_DEBUG("SUDAKOV_EXPONENT",std::cout <<"sudakov epxponent args q0: " << sqrt(q20)  <<" qh: " << sqrt(q2h)  <<" ql: " << sqrt(q2l)  <<" flav: " << flav <<  std::endl;)
	return SUDAKOV_EXPONENT(q20,q2h,q2l,flav);
};
