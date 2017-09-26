/*
 * cachedFunction.cpp
 *
 *  Created on: 16 Sep 2017
 *      Author: daniel
 */


#include "cachedFunction.h"
#include "debug.h"


template <class argT> bool cache<argT>::get(const argT& args,double& value){
	typename std::map<argT,double>::iterator it = d_store.find(args);
	if (it != d_store.end()){
		value=d_store[args];
		return true;
	}
	return false;

}


template <class argT> void cache<argT>::set(const argT& args,double value){
	NAMED_DEBUG("CACHE_STATS", std::cout << "cache args size!" << d_args.size() << std::endl;  )
	NAMED_DEBUG("CACHE_STATS", std::cout << "cache storage size!" << d_store.size() << std::endl;  )
	if ( d_args.size()==d_maxSize ){
		typename std::map<argT,double>::iterator it = d_store.find(d_args.back());
		d_store.erase (it);
		d_args.pop_back();
	}
	d_store[args]=value;
	d_args.push_front(args);
}

bool operator<(const args_t& a1,const args_t& a2){
	COMP(lambda);
	COMP(mode);
	COMP(q20);
	COMP(q2h);
	COMP(q2l);
	COMP(flav);
	return false;
}

bool operator<(const arge_t& a1,const arge_t& a2){
	COMP(lambda);
	COMP(mode);
	COMP(q20);
	COMP(q2h);
	COMP(q2l);
	COMP(q2ren);
	COMP(flav);
	return false;
}


template class cache<args_t>;
template class cache<arge_t>;
