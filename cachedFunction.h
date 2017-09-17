/*
 * cachedFunction.h
 *
 *  Created on: 16 Sep 2017
 *      Author: daniel
 */

#ifndef CACHEDFUNCTION_H_
#define CACHEDFUNCTION_H_


//#include <tuple>
#include <map>
#include <deque>
#include "MinloInfo.h"

//typedef std::tuple<double,int,double,double,double,int> args_t;
//typedef std::tuple<double,int,double,double,double,double,int> arge_t;

struct args_t {
	double lambda;
	int mode;
	double q20;
	double q2h;
	double q2l;
	int flav;
	args_t(double la,int m,double z,double h,double l,int f ){
		lambda=la; mode=m; q20=z; q2h=h; q2l=l; flav=f;
	}
};

struct arge_t {
	double lambda;
	int mode;
	double q20;
	double q2h;
	double q2l;
	double q2ren;
	int flav;
	arge_t(double la,int m,double z,double h,double l,double r,int f ){
		lambda=la; mode=m; q20=z; q2h=h; q2l=l; q2ren=r; flav=f;
	}
};


#define COMP(mem) 	if (a1.mem <a2.mem ) return true; if (a2.mem < a1.mem ) return false;


bool operator<(const args_t& a1,const args_t& a2);

bool operator<(const arge_t& a1,const arge_t& a2);



template <class argT> class cache {
	std::map<argT,double> d_store;
	std::deque<argT> d_args;
	size_t d_maxSize;
public:
	cache(int size): d_maxSize(size){}
	bool get(const argT&,double& value);
	void set(const argT&,double value);
};


//nll_sudakov_withInfo(const MinloInfo& MI,const  double &q20,const double &q2h,const double &q2l,const int &flav)

#endif /* CACHEDFUNCTION_H_ */
