#include <iostream>
#include "coreScales.h"
#include "debug.h"


double coreScale(MinloInfo::coreScaleChoice cst,const TLorentzVector& basicProcess4Vector,std::vector<fastjet::PseudoJet>& jetsLeft,const std::vector<TLorentzVector>& nonPartons){
	if (cst==MinloInfo::hthalf){
		//double ht=sqrt(basicProcess4Vector.Perp2()+80.385*80.385);
		//this is stefan's version using the invariant mass of the W instead of the pole mass
		double ht=sqrt(basicProcess4Vector.Perp2()+basicProcess4Vector.M2());
		for (int ij=0;ij<jetsLeft.size();ij++){
			fastjet::PseudoJet& j = jetsLeft[ij];
			NAMED_DEBUG("CORE_PROCESS_SCALE",
				std::cout << "now adding jet pt " << j.pt() << std::endl;
				std::cout << "ht so far: " << ht << std::endl;);
				ht+=j.pt();
		}
		NAMED_DEBUG("CORE_PROCESS_SCALE", std::cout << "final core process scale (ht/2): " << ht/2 << " Q^2: "<< ht*ht/4 <<std::endl;);
		return ht/2;
	}
	if (cst==MinloInfo::shat){
		TLorentzVector coreProcess(basicProcess4Vector);
		for (int ij=0;ij<jetsLeft.size();ij++){
			fastjet::PseudoJet& j = jetsLeft[ij];
			NAMED_DEBUG("CORE_PROCESS_SCALE",
					std::cout << "core process vector so far: " << coreProcess.E() <<" " << coreProcess.X() << " " << coreProcess.Y() << " " <<  coreProcess.Z() << std::endl;
				std::cout << "now adding jet " << j.E() <<" " << j.px() << " " << j.py() << " " <<  j.pz() << std::endl;
				std::cout << "core process scale so far: " << coreProcess.M() << std::endl;);
				coreProcess+=TLorentzVector(j.px(),j.py(),j.pz(),j.E());
		}
		NAMED_DEBUG("CORE_PROCESS_SCALE", std::cout << "final core process scale: " << coreProcess.M() << " Q^2: "<< coreProcess.M() *coreProcess.M() <<std::endl;);
		return coreProcess.M();
	}
	if (cst==MinloInfo::stefan){
		double phi_enu=nonPartons[0].DeltaPhi(nonPartons[1]);
		double cos_phi_enu =1 - cos(phi_enu);
	    //double mt = sqrt(2*nonPartons[0].Perp()*nonPartons[1].Perp()*cos_phi_enu);
	    double et=sqrt(basicProcess4Vector.Perp2()+basicProcess4Vector.M2());
	    double st=0;
	    for (int ij=0;ij<jetsLeft.size();ij++){
			fastjet::PseudoJet& j = jetsLeft[ij];
			st+=j.perp();
		}
	    double scale=st/2+et;
	    NAMED_DEBUG("CORE_PROCESS_SCALE", std::cout << "final core process scale: " << scale << " Q^2: "<< scale *scale <<std::endl;);
		return scale;
	}
}
