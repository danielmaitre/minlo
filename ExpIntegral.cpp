/*
 * ExpIntegral.cpp
 *
 *  Created on: 25 Jul 2017
 *      Author: daniel
 */

#include <math.h>
#include "ExpIntegral.h"

double expIntegral(double x){
	if (x<0){
		throw;
	}

	double tot=log(x)+EulerGamma;

	double term=x;
	int max=floor(15+x*4);
	for (int i=1;i<max;i++){
		tot+=term;
		term*=x;
		term*=i;
		term/=((i+1)*(i+1));
	}
	return tot;
}

double HF111222(double x){
	double tot=1;

	double term=x/8.0;
	int max=floor(15+x*4);
	for (int i=1;i<max;i++){
		tot+=term;
		term*=x;
		term*=(i+1)*(i+1);
		term/=((i+2)*(i+2)*(i+2));
	}
	return tot;
}

#define Power(x,n) pow(x,n)
#define Log(x) log(x)

const double Pi=M_PI;

double K(int nf){
		return 3.*(67./18.-M_PI*M_PI/6.)-10./9.*nf/2.;
};

double numExpSudakovMode2(double q2,double Q2,int nf,double lambda2){
	double b0=(33-2*nf)/(12.0*M_PI);
	double b1=(153 - 19*nf) / (24*M_PI*M_PI);
	double logQ2Lambda2=log(Q2/lambda2);
	double t=log(q2/lambda2);
	double logt=log(t);

	return (-108*Power(b0,3)*Pi*
		      Power(t,2)*
		      (3*b1 -
		        2*b1*logQ2Lambda2 -
		        b1*Power(logt,2)*t +
		        2*Power(b0,2)*
		         Power(t,2) -
		        (-3 + 2*logQ2Lambda2)*
		         logt*
		         (b1 + Power(b0,2)*t))
		       + (18*Power(b1,2)*
		         Power(logt,2)*
		         (-2*logQ2Lambda2 +
		           3*t) +
		        27*b1*t*
		         (b1 -
		           8*Power(b0,2)*t) -
		        2*logQ2Lambda2*
		         (4*Power(b1,2) -
		           27*Power(b0,2)*b1*
		           t +
		           54*Power(b0,4)*
		            Power(t,2)) +
		        6*logt*
		         (2*b1*logQ2Lambda2*
		            (-2*b1 +
		            9*Power(b0,2)*t)\
		            + 9*t*
		            (Power(b1,2) -
		            4*Power(b0,2)*b1*
		           t -
		            2*Power(b0,4)*
		            Power(t,2))))*
		      K(nf))/
		   (324.*Power(b0,6)*
		     Power(Pi,2)*Power(t,3));
}

double numExpSudakovGluonMode2(double q2,double Q2,int nf,double lambda2){
	double b0=(33-2*nf)/(12.0*M_PI);
	double b1=(153 - 19*nf) / (24*M_PI*M_PI);
	double logQ2Lambda2=log(Q2/lambda2);
	double t=log(q2/lambda2);
	double logt=log(t);
	return (-12*Power(b0,3)*Pi*
		      Power(t,2)*
		      (23*b1 -
		        18*b1*logQ2Lambda2 -
		        9*b1*Power(logt,2)*
		         t +
		        18*Power(b0,2)*
		         Power(t,2) -
		        (-23 +
		           18*logQ2Lambda2)*
		         logt*
		         (b1 + Power(b0,2)*t))
		       + (18*Power(b1,2)*
		         Power(logt,2)*
		         (-2*logQ2Lambda2 +
		           3*t) +
		        27*b1*t*
		         (b1 -
		           8*Power(b0,2)*t) -
		        2*logQ2Lambda2*
		         (4*Power(b1,2) -
		           27*Power(b0,2)*b1*
		           t +
		           54*Power(b0,4)*
		            Power(t,2)) +
		        6*logt*
		         (2*b1*logQ2Lambda2*
		            (-2*b1 +
		            9*Power(b0,2)*t)\
		            + 9*t*
		            (Power(b1,2) -
		            4*Power(b0,2)*b1*
		           t -
		            2*Power(b0,4)*
		            Power(t,2))))*
		      K(nf))/
		   (144.*Power(b0,6)*
		     Power(Pi,2)*Power(t,3));
}

double numExpSudakovMode3(double q2,double Q2,int nf,double lambda2){
	double b0=(33-2*nf)/(12.0*M_PI);
		double b1=(153 - 19*nf) / (24*M_PI*M_PI);
		double logQ2Lambda2=log(Q2/lambda2);
		double t=log(q2/lambda2);
		double logt=log(t);
		double EIt=expIntegral(t);
		double EIt2=expIntegral(t/2.0);
		double HFt=HF111222(t);
		double HFt2=HF111222(t/2);
		double lambda=sqrt(lambda2);
		double Q=sqrt(Q2);
		double expt=exp(t);
		const double E=2.7182818284590452354;
		return (-324*Power(b0,5)*EIt*lambda2*
			      Pi*Power(t,3) +
			     648*Power(b0,5)*EIt2*
			      lambda*Pi*Q*Power(t,3)\
			      + 108*Power(b0,3)*b1*
			      Power(logt,2)*Pi*
			      Power(Q,2)*Power(t,3) -
			     216*Power(b0,5)*Pi*
			      Power(Q,2)*Power(t,4) -
			     162*Power(b0,3)*b1*
			      lambda2*Pi*Power(t,2)*
			      (2*expt*(1 + logt) +
			        t*(-2 +
			           2*EulerGamma +
			           2*EulerGamma*
			           logt +
			           Power(logt,2) -
			           2*EIt*(1 + logt) +
			           2*HFt*t)) +
			     108*Power(b0,2)*b1*
			      (1 + logt)*Power(Q,2)*
			      Power(t,2)*
			      (b0*(-3 +
			           2*logQ2Lambda2)*Pi\
			         - 2*K(nf)) -
			     4*Power(b1,2)*
			      logQ2Lambda2*
			      (2 + 6*logt +
			        9*Power(logt,2))*
			      Power(Q,2)*K(nf) +
			     54*Power(b0,2)*b1*
			      logQ2Lambda2*
			      (1 + 2*logt)*Power(Q,2)*
			      t*K(nf) +
			     27*Power(b1,2)*
			      (1 + 2*logt*(1 + logt))*
			      Power(Q,2)*t*K(nf) -
			     108*Power(b0,4)*
			      logQ2Lambda2*Power(Q,2)*
			      Power(t,2)*K(nf) -
			     108*Power(b0,4)*logt*
			      Power(Q,2)*Power(t,3)*
			      (b0*(3 -
			           2*logQ2Lambda2)*Pi\
			         + K(nf)) +
			     162*Power(b0,3)*b1*
			      lambda*Pi*Q*Power(t,2)*
			      (4*Power(E,t/2.)*
			         (1 + logt) -
			        t*(2 - 2*EulerGamma +
			           2*EIt2*
			           (1 + logt) -
			           HFt2*t + Log(4) -
			           logt*
			            (2*EulerGamma +
			            Log(t/4.)))))/
			   (324.*Power(b0,6)*
			     Power(Pi,2)*Power(Q,2)*
			     Power(t,3));
}




double numExpSudakovGluonMode3(double q2,double Q2,int nf,double lambda2){
	double b0=(33-2*nf)/(12.0*M_PI);
		double b1=(153 - 19*nf) / (24*M_PI*M_PI);
		double logQ2Lambda2=log(Q2/lambda2);
		double t=log(q2/lambda2);
		double logt=log(t);
		double EIt=expIntegral(t);
		double EIt2=expIntegral(t/2.0);
		double EI3t2=expIntegral(3*t/2.0);
		double EI2t=expIntegral(2*t);
		double HFt=HF111222(t);
		double HFt2=HF111222(t/2);
		double HF3t2=HF111222(3*t/2);
		double HF2t=HF111222(2*t);
		double lambda=sqrt(lambda2);
		double Q=sqrt(Q2);
		double expt=exp(t);
		const double E=2.7182818284590452354;
return  (108*Power(b0,3)*b1*
	      Power(logt,2)*Pi +
	     (72*Power(b0,5)*EI2t*
	        Power(lambda2,2)*Pi)/
	      Power(Q,4) -
	     (192*Power(b0,5)*EI3t2*
	        Power(lambda2,1.5)*Pi)
	       /Power(Q,3) -
	     (108*Power(b0,5)*EIt*
	        lambda2*Pi)/Power(Q,2)
	       + (504*Power(b0,5)*
	        EIt2*lambda*Pi)/Q -
	     216*Power(b0,5)*Pi*t -
	     (54*Power(b0,3)*b1*
	        lambda2*Pi*
	        (2*expt*(1 + logt) +
	          t*(-2 +
	            2*EulerGamma +
	            2*EulerGamma*
	           logt +
	            Power(logt,2) -
	            2*EIt*
	           (1 + logt) +
	            2*HFt*t)))/
	      (Power(Q,2)*t) +
	     (12*Power(b0,2)*b1*
	        (1 + logt)*
	        (b0*(-23 +
	            18*logQ2Lambda2)*
	           Pi - 18*K(nf)))/t\
	      + 12*Power(b0,4)*logt*
	      (b0*(-23 +
	           18*logQ2Lambda2)*Pi
	          - 9*K(nf)) -
	     (4*Power(b1,2)*
	        logQ2Lambda2*
	        (2 + 6*logt +
	          9*Power(logt,2))*
	        K(nf))/Power(t,3) +
	     (54*Power(b0,2)*b1*
	        logQ2Lambda2*
	        (1 + 2*logt)*K(nf))/
	      Power(t,2) +
	     (27*Power(b1,2)*
	        (1 +
	          2*logt*(1 + logt))*
	        K(nf))/Power(t,2) -
	     (108*Power(b0,4)*
	        logQ2Lambda2*K(nf))/t\
	      - (48*Power(b0,3)*b1*
	        Power(lambda2,1.5)*Pi*
	        (4*Power(E,(3*t)/2.)*
	           (1 + logt) +
	          3*t*
	           (-2 +
	            2*EulerGamma -
	            2*EI3t2*
	           (1 + logt) +
	            3*HF3t2*t +
	            Log(2.25) +
	            logt*
	            (2*EulerGamma +
	            logt + Log(2.25)))
	          ))/(Power(Q,3)*t) +
	     (72*Power(b0,3)*b1*
	        Power(lambda2,2)*Pi*
	        (Power(E,2*t)*
	           (1 + logt) +
	          t*(-2 +
	            2*EulerGamma -
	            2*EI2t*
	           (1 + logt) +
	            4*HF2t*t +
	            Log(4) +
	            logt*
	            (2*EulerGamma +
	            logt + Log(4)))))/
	      (Power(Q,4)*t) -
	     (126*Power(b0,3)*b1*
	        lambda*Pi*
	        (-4*Power(E,t/2.)*
	           (1 + logt) +
	          t*(2 -
	            2*EulerGamma +
	            2*EIt2*
	           (1 + logt) -
	            HFt2*t + Log(4) -
	            logt*
	            (2*EulerGamma +
	            Log(t/4.)))))/
	      (Q*t))/
	   (144.*Power(b0,6)*
	     Power(Pi,2));

}

double numExpSudakovSubMode3(double q2,double Q2,double Q2ren,int nf,double lambda2){
	double b0=(33-2*nf)/(12.0*M_PI);
		double b1=(153 - 19*nf) / (24*M_PI*M_PI);
		double logQ2Lambda2=log(Q2/lambda2);
		double logQ2renLambda2=log(Q2ren/lambda2);
		double t=log(q2/lambda2);
		double logt=log(t);
		double lambda=sqrt(lambda2);
		double Q=sqrt(Q2);
		double expt=exp(t);
		const double E=2.7182818284590452354;
		return -((3*expt*lambda2 -
		        12*Power(E,t/2.)*
		         lambda*Q -
		        2*logQ2Lambda2*
		         Power(Q,2)*t +
		        Power(Q,2)*t*(3 + t))*
		      (Power(b0,2)*
		         logQ2renLambda2 -
		        b1*Log(logQ2renLambda2)))/
		   (3.*Power(b0,3)*
		     Power(logQ2renLambda2,2)*
		     Pi*Power(Q,2));
}



double numExpSudakovSubMode2(double q2,double Q2,double Q2ren,int nf,double lambda2){
	double b0=(33-2*nf)/(12.0*M_PI);
	double b1=(153 - 19*nf) / (24*M_PI*M_PI);
	double logQ2Lambda2=log(Q2/lambda2);
	double logQ2renLambda2=log(Q2ren/lambda2);
	double t=log(q2/lambda2);
	//double logt=log(t);
	//double lambda=sqrt(lambda2);
	//double Q=sqrt(Q2);

	return	-(t*(3 - 2*logQ2Lambda2 + t)*
		      (Power(b0,2)*
		         logQ2renLambda2 -
		        b1*Log(logQ2renLambda2)))/
		   (3.*Power(b0,3)*
		     Power(logQ2renLambda2,2)*
		     Pi);
}

double numExpSudakovGluonSubMode3(double q2,double Q2,double Q2ren,int nf,double lambda2){
	double b0=(33-2*nf)/(12.0*M_PI);
		double b1=(153 - 19*nf) / (24*M_PI*M_PI);
		double logQ2Lambda2=log(Q2/lambda2);
		double logQ2renLambda2=log(Q2ren/lambda2);
		double t=log(q2/lambda2);
		double logt=log(t);
		double lambda=sqrt(lambda2);
		double Q=sqrt(Q2);
		double expt=exp(t);
		const double E=2.7182818284590452354;
return -((-9*Power(E,2*t)*
        Power(lambda2,2) +
       32*Power(E,(3*t)/2.)*
        Power(lambda2,1.5)*Q\
        + 27*expt*lambda2*
        Power(Q,2) -
       252*Power(E,t/2.)*
        lambda*Power(Q,3) -
       54*logQ2Lambda2*
        Power(Q,4)*t +
       3*Power(Q,4)*t*
        (23 + 9*t))*
     (Power(b0,2)*
        logQ2renLambda2 -
       b1*Log(logQ2renLambda2)))/
  (36.*Power(b0,3)*
    Power(logQ2renLambda2,2)*
    Pi*Power(Q,4));}

double numExpSudakovGluonSubMode2(double q2,double Q2,double Q2ren,int nf,double lambda2){
	double b0=(33-2*nf)/(12.0*M_PI);
		double b1=(153 - 19*nf) / (24*M_PI*M_PI);
		double logQ2Lambda2=log(Q2/lambda2);
		double logQ2renLambda2=log(Q2ren/lambda2);
		double t=log(q2/lambda2);

return -(t*(23 - 18*logQ2Lambda2 +
        9*t)*
      (Power(b0,2)*
         logQ2renLambda2 -
        b1*Log(logQ2renLambda2)))/
   (12.*Power(b0,3)*
     Power(logQ2renLambda2,2)*
     Pi);
}


double analyticSudakovGluonMode2(double q02,double Q2,double q2,double lambda2){
	int nf=5;
	double g1=numExpSudakovGluonMode2(q02,q2,nf,lambda2);
	double g2=numExpSudakovGluonMode2(q2,q2,nf,lambda2);
	double g3=numExpSudakovGluonMode2(q02,Q2,nf,lambda2);
	double g4=numExpSudakovGluonMode2(Q2,Q2,nf,lambda2);
	return exp(g3-g4-g1+g2);
}
double analyticSudakovQuarkMode2(double q02,double Q2,double q2,double lambda2){
	int nf=5;
	double g1=numExpSudakovMode2(q02,q2,nf,lambda2);
	double g2=numExpSudakovMode2(q2,q2,nf,lambda2);
	double g3=numExpSudakovMode2(q02,Q2,nf,lambda2);
	double g4=numExpSudakovMode2(Q2,Q2,nf,lambda2);
	return exp(g3-g4-g1+g2);
}

double analyticSudakovGluonMode3(double q02,double Q2,double q2,double lambda2){
	int nf=5;
	double g1=numExpSudakovGluonMode3(q02,q2,nf,lambda2);
	double g2=numExpSudakovGluonMode3(q2,q2,nf,lambda2);
	double g3=numExpSudakovGluonMode3(q02,Q2,nf,lambda2);
	double g4=numExpSudakovGluonMode3(Q2,Q2,nf,lambda2);
	return exp(g3-g4-g1+g2);
}
double analyticSudakovQuarkMode3(double q02,double Q2,double q2,double lambda2){
	int nf=5;
	double g1=numExpSudakovMode3(q02,q2,nf,lambda2);
	double g2=numExpSudakovMode3(q2,q2,nf,lambda2);
	double g3=numExpSudakovMode3(q02,Q2,nf,lambda2);
	double g4=numExpSudakovMode3(Q2,Q2,nf,lambda2);
	return exp(g3-g4-g1+g2);
}

double analyticSudakovSubGluonMode2(double q02,double Q2,double q2,double Q2ren,double lambda2){
	int nf=5;
	double g1=numExpSudakovGluonSubMode2(q02,q2,Q2ren,nf,lambda2);
	double g2=numExpSudakovGluonSubMode2(q2,q2,Q2ren,nf,lambda2);
	double g3=numExpSudakovGluonSubMode2(q02,Q2,Q2ren,nf,lambda2);
	double g4=numExpSudakovGluonSubMode2(Q2,Q2,Q2ren,nf,lambda2);
	return g4-g3-g2+g1;
}
double analyticSudakovSubQuarkMode2(double q02,double Q2,double q2,double Q2ren,double lambda2){
	int nf=5;
	double g1=numExpSudakovSubMode2(q02,q2,Q2ren,nf,lambda2);
	double g2=numExpSudakovSubMode2(q2,q2,Q2ren,nf,lambda2);
	double g3=numExpSudakovSubMode2(q02,Q2,Q2ren,nf,lambda2);
	double g4=numExpSudakovSubMode2(Q2,Q2,Q2ren,nf,lambda2);
	return g4-g3-g2+g1;
}

double analyticSudakovSubGluonMode3(double q02,double Q2,double q2,double Q2ren,double lambda2){
	int nf=5;
	double g1=numExpSudakovGluonSubMode3(q02,q2,Q2ren,nf,lambda2);
	double g2=numExpSudakovGluonSubMode3(q2,q2,Q2ren,nf,lambda2);
	double g3=numExpSudakovGluonSubMode3(q02,Q2,Q2ren,nf,lambda2);
	double g4=numExpSudakovGluonSubMode3(Q2,Q2,Q2ren,nf,lambda2);
	return g4-g3-g2+g1;
}
double analyticSudakovSubQuarkMode3(double q02,double Q2,double q2,double Q2ren,double lambda2){
	int nf=5;
	double g1=numExpSudakovSubMode3(q02,q2,Q2ren,nf,lambda2);
	double g2=numExpSudakovSubMode3(q2,q2,Q2ren,nf,lambda2);
	double g3=numExpSudakovSubMode3(q02,Q2,Q2ren,nf,lambda2);
	double g4=numExpSudakovSubMode3(Q2,Q2,Q2ren,nf,lambda2);
	return g4-g3-g2+g1;
}

