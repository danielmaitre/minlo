/*
 * nll.h
 *
 *  Created on: 31 Mar 2017
 *      Author: daniel
 */

#ifndef NLL_H_
#define NLL_H_

#include "MinloInfo.h"

double nll_sudakov_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav);
double nll_exponent_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav);
double fo_exponent_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const double &q2ren,const int &flav);



#endif /* NLL_H_ */
