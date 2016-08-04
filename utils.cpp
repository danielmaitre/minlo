#include "fastjet/ClusterSequence.hh"
#include "MyFlavKtPlugin/MyFlavKtPlugin.hh"
#include "FlavKtPlugin/FlavKtPlugin.hh"
#include <vector>
#include "debug.h"
//const fastjet::FlavInfo & flavour_of(const fastjet::PseudoJet& j){
//	if (j.has_user_info<fastjet::FlavInfo>()) {
//		return j.user_info<fastjet::FlavInfo>();
//	} else {
//	   return fastjet::_no_flav;
// 	}
//}


// assumes the object is either a gluon or a quark, if multiflavoured returns 0
int pdgFromFlavor(const fastjet::FlavInfo& fi){
	if (fi.is_multiflavored()){
		return 0;
	}
	for (int i=1;i<=6;i++){
		if (fi[i]==1) return i;
		if (fi[i]==-1) return -i;
	}
	return 21;
}


bool hasQuarks(const fastjet::ClusterSequence& cs,int njets){

	std::vector<fastjet::PseudoJet> jetsLeft=cs.exclusive_jets(njets/*njetsStart-nbrClusteringsDone*/);
	NAMED_DEBUG("HASQUARK",std::cout << "Testing "<< njets << " particles" << std::endl; )
	int nQ=0;
	for (int ij=0;ij<jetsLeft.size();ij++){
		const fastjet::PseudoJet& j=jetsLeft[ij];
		const fastjet::FlavInfo & flavour = fastjet::MyFlavKtPlugin::flavour_of(j);
		NAMED_DEBUG("HASQUARK",std::cout << flavour.description() << std::endl; )
		int pdg=pdgFromFlavor(flavour);
		if ( pdg!=21 and pdg!=0 ){
			nQ+=1;
		}
	}
	const fastjet::MyFlavKtPlugin::Extras * extras =
			dynamic_cast<const fastjet::MyFlavKtPlugin::Extras *>(cs.extras());
	int b1=pdgFromFlavor(extras->beam_flav_forward(njets));
	int b2=pdgFromFlavor(extras->beam_flav_backward(njets));

	if ( b1!=21 and b1!=0 ){nQ+=1;}
	if ( b2!=21 and b2!=0 ){nQ+=1;}

	NAMED_DEBUG("HASQUARK",
		std::cout << "forward: " << b1 << "backward:" << b2 << std::endl;
	)
	;
	if (nQ>=2){
	NAMED_DEBUG("HASQUARK",
		std::cout << "number of quarks/antiquarks: " << nQ << " --> test passed"  << std::endl;
	)
		return true;
	} else {
	NAMED_DEBUG("HASQUARK",
		std::cout << "number of quarks/antiquarks: " << nQ << " --> test failed"  << std::endl;
	)
		return false;
	}
}
