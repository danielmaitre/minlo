/*
 * MINLOfunctions.h
 *
 *  Created on: 28 Feb 2016
 *      Author: daniel
 */

#ifndef MINLOFUNCTIONS_H_
#define MINLOFUNCTIONS_H_

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include <iomanip>
#include "TLorentzVector.h"
#include "MinloInfo.h"


class minloImpl{
public:
	static double g_MinloScale;
	static double g_MinloFactor;
	static bool g_MinloScaleValid;
	static double g_nclusterings;  //needs to be a double for simplicity of the interface
};




struct sudakovCandidate {
	double highScale;
	double lowScale;
	int historyIndex;  // this is where the recombination would happen, but it can be that the clustering stops before that
	int flavor;
};

std::ostream & operator<<(std::ostream& ostr, const sudakovCandidate & sc) ;



std::ostream & operator<<(std::ostream& ostr, const fastjet::PseudoJet & j) ;





double MINLOcomputeSudakov(const MinloInfo& MI,const NtupleInfo<MAX_NBR_PARTICLES>& Ev,double &q0,double &scaleForNLO,int &status,bool useNewNtupleFormat,bool useDouble);
void fillJetVector(NtupleInfo<MAX_NBR_PARTICLES>& Ev,std::vector<fastjet::PseudoJet>& particles,const int* flavors,bool useFlavor,bool useDouble);
NtupleInfo<MAX_NBR_PARTICLES> boostedToCMF(const NtupleInfo<MAX_NBR_PARTICLES>& orig,bool useDouble);
double getBeamEnergy(const NtupleInfo<MAX_NBR_PARTICLES>& Ev,bool useDouble);
void printEventInclView(fastjet::ClusterSequence& cs,double ptCut);
void printEventExclView(int n_excl,fastjet::ClusterSequence& cs);
void printFlavourPart(const std::vector<fastjet::PseudoJet>& input_particles);
double MINLO_computeSudakovKeith(const NtupleInfo<MAX_NBR_PARTICLES>& Ev,int flg_bornonly,int imode,int isReal,double beamEnergy,int nlegborn,int st_bornorder,bool useMinloIDs,bool useDouble,int alltheway);
double getMINLOweight(MinloInfo& MI,const NtupleInfo<MAX_NBR_PARTICLES>& Ev,int weightType);


#endif /* MINLOFUNCTIONS_H_ */
