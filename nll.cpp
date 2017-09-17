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
#include "cachedFunction.h"

double nll_sudakov_withInfo_fn(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav){
	NAMED_DEBUG("SUDAKOV_ARGUMENTS",std::cout <<"sudakov args q0: " << sqrt(q20)  <<" qh: " << sqrt(q2h)  <<" ql: " << sqrt(q2l)  <<" flav: " << flav <<  std::endl;)
	if (!MI.d_useSherpa){
			return NLL_SUDAKOV(q20,q2h,q2l,flav);
	} else {
		if (MI.d_useAnalyticalSherpa){
			if (q2h<=q20){
				return 1;
			}

			double lambda2=MI.d_lambda*MI.d_lambda;
			double res;
			if (flav==0){
				if (MI.d_sherpaMode==2){
					res=analyticSudakovGluonMode2(q20, q2h, q2l, lambda2);
				}
				if (MI.d_sherpaMode==3){
					res=analyticSudakovGluonMode3(q20, q2h, q2l, lambda2);
				}
			} else {
				if (MI.d_sherpaMode==2){
					res=analyticSudakovQuarkMode2(q20, q2h, q2l, lambda2);
				}
				if (MI.d_sherpaMode==3){
					res=analyticSudakovQuarkMode3(q20, q2h, q2l, lambda2);
				}
			}
			NAMED_DEBUG("NLL_ANALYTIC",
			double resSherpa=SherpaSudakov(q20, q2h, q2l, flav, currentPDF::s_PDF, MI.d_sherpaMode);
			if (abs((res-resSherpa)/(res+resSherpa))>1e-3){
				std::cout << "Discrepancy for " <<
						"rel diff: " << abs((res-resSherpa)/(res+resSherpa)) << " " <<
						"res: " << res << " " <<
						"resSherpa: " << resSherpa << " " <<
					 	"q20: " << q20 << " " <<
					 	"q2h: " << q2h << " " <<
					 	"q2l: " << q2l << " " <<
						"flav: " << flav << " " <<
						std::endl;
			}
			NAMED_DEBUG("PRINT_ALL",std::cout << "analytic:" << res << " numeric: " << resSherpa << std::endl;
			std::cout << "Discrepancy for " <<
					"rel diff: " << abs((res-resSherpa)/(res+resSherpa)) << " " <<
					"res: " << res << " " <<
					"resSherpa: " << resSherpa << " " <<
				 	"q20: " << q20 << " " <<
				 	"q2h: " << q2h << " " <<
				 	"q2l: " << q2l << " " <<
					"flav: " << flav << " " <<
					std::endl;
)

			)
			return res;
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
	NAMED_DEBUG("SUDAKOV_EXPONENT",std::cout <<"sudakov exponent args q0: " << sqrt(q20)  <<" qh: " << sqrt(q2h)  <<" ql: " << sqrt(q2l)  <<" flav: " << flav <<  std::endl;)
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

double nll_sudakov_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav){
	static cache<args_t> valueCache(100);
	static double n_used=0;
	static double n_notused=0;
	args_t a(MI.d_lambda,MI.d_sherpaMode,q20,q2h,q2l,flav);
	double cachedRes;
	bool cached=valueCache.get(a,cachedRes);
	if (cached){
		n_used+=1;
		NAMED_DEBUG("CACHE_STATS", std::cout << "reusing! ratio:" << n_used/(n_used+n_notused) << std::endl;  )
		return cachedRes;
	}
	n_notused+=1;
	NAMED_DEBUG("CACHE_STATS", std::cout << "not reusing!" << n_used/(n_used+n_notused) << std::endl;  )
	cachedRes=nll_sudakov_withInfo_fn(MI,q20,q2h,q2l,flav);
	valueCache.set(a,cachedRes);
	return cachedRes;
}



double fo_exponent_withInfo_fn(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const double &q2ren,const int &flav){
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

double fo_exponent_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const double &q2ren,const int &flav){
	static cache<arge_t> valueCache(100);

	arge_t a(MI.d_lambda,MI.d_sherpaMode,q20,q2h,q2l,q2ren,flav);
	double cachedRes;
	bool cached=valueCache.get(a,cachedRes);
	if (cached){
		NAMED_DEBUG("CACHE_STATS", std::cout << "reusing!" << std::endl;  )
		return cachedRes;
	}
	NAMED_DEBUG("CACHE_STATS", std::cout << "not reusing!" << std::endl;  )
	cachedRes=fo_exponent_withInfo_fn(MI,q20,q2h,q2l,q2ren,flav);
	valueCache.set(a,cachedRes);
	return cachedRes;
}

