/*
 * cachedFunction.h
 *
 *  Created on: 16 Sep 2017
 *      Author: daniel
 */

#ifndef CACHEDFUNCTION_H_
#define CACHEDFUNCTION_H_


#include <tuple>
#include <map>
#include <deque>
#include "MinloInfo.h"

typedef std::tuple<double,int,double,double,double,int> args_t;
typedef std::tuple<double,int,double,double,double,double,int> arge_t;

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
