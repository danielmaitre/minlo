
#include "boosts.h"
#include "../debug.h"
#include "../printHelper.h"

using namespace std;

void ISRboost(fastjet::PseudoJet& bj1,fastjet::PseudoJet& bj2,
				const fastjet::PseudoJet& j1,
				const vector<int>& candidates,
				fastjet::ClusterSequence & cs,
				const vector<int>& beam_particles,
				double* boostzVec,double *boostxyVec
				){

			  double q0=bj1.E()+bj2.E()-j1.E();
			  double q2=q0*q0;
			  double keith_beta=j1.pz()/q0;
			  double krec2=j1.px()*j1.px()+j1.py()*j1.py()+j1.pz()*j1.pz();
			  double krec=sqrt(krec2);
			  double recMass2=q0*q0-krec2;
				boostzVec[0]=0;boostzVec[1]=0;boostzVec[2]=keith_beta;

			  double jperp=j1.perp();
			double norm=sqrt(j1.px()*j1.px()+j1.py()*j1.py()+j1.pz()*j1.pz());
			double Keith_norm=sqrt(recMass2+j1.perp2());

			double boostx=-j1.px()/Keith_norm;
			double boosty=-j1.py()/Keith_norm;
			boostxyVec[0]=-boostx;boostxyVec[1]=-boosty;boostxyVec[2]=0;



			TLorentzVector q(-j1.px(),-j1.py(),-j1.pz(),-j1.E());
			TLorentzVector m1(bj1.px(),bj1.py(),bj1.pz(),bj1.E());
			TLorentzVector m2(bj2.px(),bj2.py(),bj2.pz(),bj2.E());
			q=q+m1+m2;
			NAMED_DEBUG("BOOST",
				std::cout << "q0: " << q0 << " keith_beta: " << keith_beta <<  std::endl ;
				std::cout << "rebM2: " << recMass2 << std::endl ;
				std::cout << "keith norm: " << Keith_norm << std::endl ;
			)
			NAMED_DEBUG("BOOST",
			std::cout << setw(12) << "pt" << setw(12) << "rapidity" << setw(12) << "phi" << "  flav" << std::endl;
			)
			for (int ic=0;ic<candidates.size();ic++){
				int jetIndex=candidates[ic];
				const fastjet::PseudoJet& jc=cs.jets()[jetIndex];
				NAMED_DEBUG("BOOST",std::cout << "Jet before boost: " << jc  << std::endl; );
				fastjet::PseudoJet& j=const_cast<fastjet::PseudoJet&>(jc);
				TLorentzVector bm(j.px(),j.py(),j.pz(),j.E());
				NAMED_DEBUG("BOOST",std::cout << "Before boosts: " << bm  << std::endl; );
				bm.Boost(0,0,keith_beta);
				NAMED_DEBUG("BOOST",std::cout << "After z boost: " << bm  << std::endl; );
				bm.Boost(-boostx,-boosty,0);
				NAMED_DEBUG("BOOST",std::cout << "After xy boost: " << bm  << std::endl; );
				j.reset_momentum(bm.X(),bm.Y(),bm.Z(),bm.E());
				NAMED_DEBUG("BOOST",std::cout << "Jet after boost: " << cs.jets()[jetIndex]  << std::endl; );
			}
			NAMED_DEBUG("BOOST",std::cout << "Now boosting the beam particles" << std::endl; );

			for (int iib=0;iib<2;iib++){
				int jetIndex=beam_particles[iib];
				const fastjet::PseudoJet& jc=cs.jets()[jetIndex];
				NAMED_DEBUG("BOOST",std::cout << "Beam jet before boost: " << jc  << std::endl; );
				fastjet::PseudoJet& j=const_cast<fastjet::PseudoJet&>(jc);
				TLorentzVector bm(j.px(),j.py(),j.pz(),j.E());
				bm.Boost(0,0,keith_beta);
				NAMED_DEBUG("BOOST",std::cout << "Beam Jet after boost: " << bm  << std::endl; );
			}




			NAMED_DEBUG("BOOST",std::cout << "CMF before z boost: " << q  << std::endl; );
			q.Boost(0,0,-keith_beta);
			NAMED_DEBUG("BOOST",std::cout << "CMF after z boost: " << q  << std::endl; );
			q.Boost(-boostx,-boosty,0);
			NAMED_DEBUG("BOOST",std::cout << "CMF after xy boost: " << q  << std::endl; );

			double qz=q.E()/2.0;

			bj1.reset_momentum(0,0,qz,qz);
			bj2.reset_momentum(0,0,-qz,qz);
			NAMED_DEBUG("BOOST",std::cout << "beam Jet after boost: " << bj1  << std::endl; );
			NAMED_DEBUG("BOOST",std::cout << "beam Jet after boost: " << bj2  << std::endl; );

}


void FSRboost(fastjet::PseudoJet& bj1,fastjet::PseudoJet& bj2,
				const fastjet::PseudoJet& j1,const fastjet::PseudoJet& j2,
				fastjet::PseudoJet& j12,
				const vector<int>& candidates,
				fastjet::ClusterSequence & cs,
				const vector<int>& beam_particles,double* boostVec
				){

	double q0=2*bj1.E();
	double q2=q0*q0;
	double k0rec=q0-j1.E()-j2.E();
	double krecx=-j1.px()-j2.px();
	double krecy=-j1.py()-j2.py();
	double krecz=-j1.pz()-j2.pz();
	double krec=sqrt(krecx*krecx+krecy*krecy+krecz*krecz);
	double kk=(k0rec+krec);
	double beta=(q2-kk*kk)/(q2+kk*kk);
	double vecnx=krecx/krec;
	double vecny=krecy/krec;
	double vecnz=krecz/krec;
	double vecx=beta*vecnx;
	double vecy=beta*vecny;
	double vecz=beta*vecnz;

	boostVec[0]=vecx;
	boostVec[1]=vecy;
	boostVec[2]=vecz;

	TLorentzVector q(-j1.px()-j2.px(),-j1.py()-j2.py(),-j1.pz()-j2.pz(),-j1.E()-j2.E());
	TLorentzVector m1(bj1.px(),bj1.py(),bj1.pz(),bj1.E());
	TLorentzVector m2(bj2.px(),bj2.py(),bj2.pz(),bj2.E());
	q=q+m1+m2;
	NAMED_DEBUG("BOOST",
				std::cout << "k0rec: " << k0rec << " beta: " << beta <<  std::endl ;
			)
	NAMED_DEBUG("BOOST",
		std::cout << setw(12) << "pt" << setw(12) << "rapidity" << setw(12) << "phi" << "  flav" << std::endl;
	)
	for (int ic=0;ic<candidates.size();ic++){
		int jetIndex=candidates[ic];
		const fastjet::PseudoJet& jc=cs.jets()[jetIndex];
		NAMED_DEBUG("BOOST",std::cout << "Jet before boost: " << jc  << std::endl; );
		fastjet::PseudoJet& j=const_cast<fastjet::PseudoJet&>(jc);
		TLorentzVector bm(j.px(),j.py(),j.pz(),j.E());
		bm.Boost(vecx,vecy,vecz);
		j.reset_momentum(bm.X(),bm.Y(),bm.Z(),bm.E());
		NAMED_DEBUG("BOOST",std::cout << "Jet after boost: " << cs.jets()[jetIndex]  << std::endl; );
	}

	double k=(q2-(k0rec*k0rec-krec*krec))/(2*q0);
	j12.reset_momentum(
		-vecnx*k,
		-vecny*k,
		-vecnz*k,
		k
	);
}




