#ifndef _H_PRINTHELPER_H
#define _H_PRINTHELPER_H

#include "MyFlavKtPlugin/MyFlavKtPlugin.hh"
#include "fastjet/ClusterSequence.hh"
#include "TLorentzVector.h"
#include <iomanip>


inline std::ostream & operator<<(std::ostream& ostr, const fastjet::PseudoJet & j) {
  ostr << std::setw(12) << j.perp()
       << std::setw(12) << j.rap()
	   << std::setw(12) << j.phi()
	   << " E: " << std::setw(12) << j.E()
	   << " x,y,z: " <<  std::setw(12) << j.px()
	   << std::setw(12) << j.py()
	   << std::setw(12) << j.pz() ;
  if (j.has_user_info<fastjet::FlavInfo>()) ostr << "  " << j.user_info<fastjet::FlavInfo>().description();
  return ostr;
}

inline std::ostream & operator<<(std::ostream& ostr, const TLorentzVector & m) {
  ostr << std::setw(12) << m.Perp()
       << std::setw(12) << m.Rapidity()
	   << std::setw(12) << m.Phi()
	   << " E: " << std::setw(12) << m.E()
	   << " x,y,z: " <<  std::setw(12) << m.X()
	   << std::setw(12) << m.Y()
	   << std::setw(12) << m.Z() ;
  return ostr;
}

#endif /* PRINT_HELPER_H */
