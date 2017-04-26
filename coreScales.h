/*
 * coreScales.h
 *
 *  Created on: 24 Apr 2017
 *      Author: daniel
 */

#ifndef CORESCALES_H_
#define CORESCALES_H_

#include "MinloInfo.h"
#include "TLorentzVector.h"
#include "fastjet/JetDefinition.hh"

double coreScale(MinloInfo::coreScaleChoice cst,const TLorentzVector& basicProcess4Vector,std::vector<fastjet::PseudoJet>& jetsLeft,const std::vector<TLorentzVector>& nonPartons);



#endif /* CORESCALES_H_ */
