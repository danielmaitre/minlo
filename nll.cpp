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
#include "ExpIntegral.h"

double nll_sudakov_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav){
	NAMED_DEBUG("SUDAKOV_ARGUMENTS",std::cout <<"sudakov args q0: " << sqrt(q20)  <<" qh: " << sqrt(q2h)  <<" ql: " << sqrt(q2l)  <<" flav: " << flav <<  std::endl;)
	if (!MI.d_useSherpa){
			return NLL_SUDAKOV(q20,q2h,q2l,flav);
	} else {
		if (MI.d_useAnalyticalSherpa){
			double lambda2=MI.d_lambda*MI.d_lambda;
			if (flav==0){
				if (MI.d_sherpaMode==2){
					return analyticSudakovGluonMode2(q20, q2h, q2l, lambda2);
				}
				if (MI.d_sherpaMode==3){
					return analyticSudakovGluonMode3(q20, q2h, q2l, lambda2);
				}
			} else {
				if (MI.d_sherpaMode==2){
					return analyticSudakovQuarkMode2(q20, q2h, q2l, lambda2);
				}
				if (MI.d_sherpaMode==3){
					return analyticSudakovQuarkMode3(q20, q2h, q2l, lambda2);
				}
			}
		} else {
			double q20l=q20;
			double q2ll=q2l;
			double q2hl=q2h;
			int fl=flav;
			if (q2h<=q20){
				return 1;
			}
			return SherpaSudakov(q20, q2h, q2l, flav, currentPDF::s_PDF, MI.d_sherpaMode);
		}
	}
};

double nll_exponent_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav){
	NAMED_DEBUG("SUDAKOV_EXPONENT",std::cout <<"sudakov epxponent args q0: " << sqrt(q20)  <<" qh: " << sqrt(q2h)  <<" ql: " << sqrt(q2l)  <<" flav: " << flav <<  std::endl;)
	if (!MI.d_useSherpa){
		return SUDAKOV_EXPONENT(q20,q2h,q2l,flav);
	} else{
		double q20l=q20;
		double q2ll=q2l;
		double q2hl=q2h;
		int fl=flav;
		if (flav==0){
			return SherpaExponentGluon(q20,q2h, currentPDF::s_PDF,MI.d_sherpaMode,false,MI.d_nfgs,q2h); /*might need to add q2ren dependence...*/
		} else {
			return SherpaExponentQuark(q20,q2h, currentPDF::s_PDF,MI.d_sherpaMode,false,MI.d_nfgs,q2h);
		}
	}
};

double fo_exponent_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const double &q2ren,const int &flav){
	NAMED_DEBUG("SUDAKOV_FO_EXPONENT",std::cout <<"sudakov fo exponent args q0: " << sqrt(q20)  <<" qh: " << sqrt(q2h)  <<" ql: " << sqrt(q2l)  <<" flav: " << flav << "q2ren: " << q2ren << std::endl;)
	if (!MI.d_useSherpa){
		return EXPSUDAKOV(q20,q2h,q2l,flav);
	} else{
		if (MI.d_useAnalyticalSherpa){
			double lambda2=MI.d_lambda*MI.d_lambda;
			if (flav==0){
				if (MI.d_sherpaMode==2){
					return analyticSudakovSubGluonMode2(q20, q2h, q2l,q2ren, lambda2);
				}
				if (MI.d_sherpaMode==3){
					return analyticSudakovSubGluonMode3(q20, q2h, q2l,q2ren, lambda2);
				}
			} else {
				if (MI.d_sherpaMode==2){
					return analyticSudakovSubQuarkMode2(q20, q2h, q2l,q2ren, lambda2);
				}
				if (MI.d_sherpaMode==3){
					return analyticSudakovSubQuarkMode3(q20, q2h, q2l,q2ren, lambda2);
				}
			}
		} else {
			if (flav==0){
				double d1=-SherpaExponentGluon(q20,q2h, currentPDF::s_PDF,MI.d_sherpaMode,true,MI.d_nfgs,q2ren);
				double d2=-SherpaExponentGluon(q20,q2l, currentPDF::s_PDF,MI.d_sherpaMode,true,MI.d_nfgs,q2ren);
				NAMED_DEBUG("SUDAKOV_FO_EXPONENT",std::cout << "d1: " << d1 << " d2: " << d2 << std::endl;)
				NAMED_DEBUG("SUDAKOV_FO_EXPONENT",std::cout << "stefan: " << (d1-d2)/currentPDF::s_PDF->alphasQ2(q2ren) << " keith/stefan: " <<  EXPSUDAKOV(q20,q2h,q2l,flav) / ((d1-d2)/currentPDF::s_PDF->alphasQ2(q2ren))<< " alphas(q2ren):" << currentPDF::s_PDF->alphasQ2(q2ren) << std::endl;)
				return d1-d2;
			} else {
				double d1=-SherpaExponentQuark(q20,q2h, currentPDF::s_PDF,MI.d_sherpaMode,true,MI.d_nfgs,q2ren);
				double d2=-SherpaExponentQuark(q20,q2l, currentPDF::s_PDF,MI.d_sherpaMode,true,MI.d_nfgs,q2ren);
				NAMED_DEBUG("SUDAKOV_FO_EXPONENT",std::cout << "d1: " << d1 << " d2: " << d2 << std::endl;)
				NAMED_DEBUG("SUDAKOV_FO_EXPONENT",std::cout << "stefan: " << (d1-d2)/currentPDF::s_PDF->alphasQ2(q2ren) << " keith/stefan: " << EXPSUDAKOV(q20,q2h,q2l,flav) / ((d1-d2)/currentPDF::s_PDF->alphasQ2(q2ren))<< " alphas(q2ren):" << currentPDF::s_PDF->alphasQ2(q2ren) << std::endl;)
				return d1-d2;
			}
		}
	}
};
