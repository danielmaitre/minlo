#ifndef _H_MINLOINFO_H_
#define _H_MINLOINFO_H_

#include <iomanip>

struct MinloInfo {
	double d_energy;
	int d_njetsOrig;
	int d_njetsClus;
	enum wType { born,bornLO,nlo,real};
	wType d_type;

  void print(std::ostream& os) {
    os << "njetsOrig  : " << d_njetsOrig << std::endl;
    os << "njetsClus  : " << d_njetsClus << std::endl;
    os << "beam energy: " << d_energy << std::endl;
    if ( d_type == MinloInfo::born){
      os << "type  : born" << std::endl;
    }
    if ( d_type == MinloInfo::nlo){
      os << "type  : nlo" << std::endl;
    }
    if ( d_type == MinloInfo::bornLO){
      os << "type  : bornLO" << std::endl;
    }
    if ( d_type == MinloInfo::real){
      os << "type  : real" << std::endl;
    }
  }
};


#endif /* _H_MINLOINFO_H_ */

