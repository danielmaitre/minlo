/*
 * ExpIntegral.h
 *
 *  Created on: 26 Jul 2017
 *      Author: daniel
 */

#ifndef EXPINTEGRAL_H_
#define EXPINTEGRAL_H_

double expIntegral(double x);
double HF111222(double x);
double numExpSudakovMode2(double q2,double Q2,int nf,double lambda2);
double numExpSudakovMode3(double q2,double Q2,int nf,double lambda2);
double numExpSudakovGluonMode2(double q2,double Q2,int nf,double lambda2);
double numExpSudakovGluonMode3(double q2,double Q2,int nf,double lambda2);

double analyticSudakovQuarkMode2(double q02,double Q2,double q2,double lambda2);
double analyticSudakovQuarkMode3(double q02,double Q2,double q2,double lambda2);
double analyticSudakovGluonMode2(double q02,double Q2,double q2,double lambda2);
double analyticSudakovGluonMode3(double q02,double Q2,double q2,double lambda2);

double analyticSudakovSubQuarkMode2(double q02,double Q2,double q2,double Q2ren,double lambda2);
double analyticSudakovSubQuarkMode3(double q02,double Q2,double q2,double Q2ren,double lambda2);
double analyticSudakovSubGluonMode2(double q02,double Q2,double q2,double Q2ren,double lambda2);
double analyticSudakovSubGluonMode3(double q02,double Q2,double q2,double Q2ren,double lambda2);


double numExpSudakovGluonSubMode3(double q2,double Q2,double Q2ren,int nf,double lambda2);
double numExpSudakovSubMode3(double q2,double Q2,double Q2ren,int nf,double lambda2);
double numExpSudakovGluonSubMode2(double q2,double Q2,double Q2ren,int nf,double lambda2);
double numExpSudakovSubMode2(double q2,double Q2,double Q2ren,int nf,double lambda2);

const double EulerGamma=0.57721566490153286061;

#endif /* EXPINTEGRAL_H_ */
