#ifndef BOOSTS_H_
#define BOOSTS_H_

#include "fastjet/ClusterSequence.hh"
#include "TLorentzVector.h"
#include <vector>

void ISRboost(fastjet::PseudoJet& bj1,fastjet::PseudoJet& bj2,
	      const fastjet::PseudoJet& j1,
	      const std::vector<int>& candidates,
	      fastjet::ClusterSequence & cs,
	      const std::vector<int>& beam_particles
	      );

void FSRboost(fastjet::PseudoJet& bj1,fastjet::PseudoJet& bj2,
	      const fastjet::PseudoJet& j1,
	      const fastjet::PseudoJet& j2,
	      fastjet::PseudoJet& j12,
	      const std::vector<int>& candidates,
	      fastjet::ClusterSequence & cs,
	      const std::vector<int>& beam_particles
	      );


#endif /* BOOSTS_H_ */
