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
#include "MyFlavKtPlugin/MyFlavKtPlugin.hh"
#include <iomanip>
#include "TLorentzVector.h"



class minloImpl{
public:
	static double g_MinloScale;
	static double g_MinloFactor;
	static bool g_MinloScaleValid;
	static double g_nclusterings;  //needs to be a double for simplicity of the interface
};


struct MinloInfo {
	double d_energy;
	int d_njetsOrig;
	int d_njetsClus;
	enum wType { born,bornLO,nlo,real};
	wType d_type;
};


struct sudakovCandidate {
	double highScale;
	double lowScale;
	int historyIndex;  // this is where the recombination would happen, but it can be that the clustering stops before that
	int flavor;
};




std::ostream & operator<<(std::ostream& ostr, const fastjet::PseudoJet & j) ;


typedef fastjet::MyFlavKtPlugin THEPLUGIN;


double MINLOcomputeSudakov(MinloInfo& MI,NtupleInfo<MAX_NBR_PARTICLES>& Ev, int weightType,double &q0,double &scaleForNLO);
void fillJetVector(NtupleInfo<MAX_NBR_PARTICLES>& Ev,std::vector<fastjet::PseudoJet>& particles,int* flavors,bool useFlavor);
NtupleInfo<MAX_NBR_PARTICLES> boostedToCMF(NtupleInfo<MAX_NBR_PARTICLES>& orig);
double getBeamEnergy(NtupleInfo<MAX_NBR_PARTICLES>& Ev);
void printEventInclView(fastjet::ClusterSequence& cs,double ptCut);
void printEventExclView(int n_excl,fastjet::ClusterSequence& cs);
void printFlavourPart(const std::vector<fastjet::PseudoJet>& input_particles);
double MINLO_computeSudakovKeith(NtupleInfo<MAX_NBR_PARTICLES>& Ev,int flg_bornonly,int imode,int isReal);
double getMINLOweight(MinloInfo& MI,NtupleInfo<MAX_NBR_PARTICLES>& Ev,int weightType);

#endif /* MINLOFUNCTIONS_H_ */
