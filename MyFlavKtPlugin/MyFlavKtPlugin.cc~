#include "MyFlavKtPlugin.hh"
#include "fastjet/ClusterSequence.hh"
#include "ktflavf77.h"
#include <sstream>
#include <utility>
#include "../debug.h"
#include <iterator>
#include <cmath> //for log
#include "TLorentzVector.h"
#include <iomanip>


using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

inline std::ostream & operator<<(std::ostream& ostr, const PseudoJet & j) {
  ostr << std::setw(12) << j.perp()
       << std::setw(12) << j.rap()
	   << std::setw(12) << j.phi()
	   << std::setw(12) << j.px()
	   << " E: " << std::setw(12) << j.E()
	   << " x,y,z: " <<  std::setw(12) << j.px()
	   << std::setw(12) << j.py()
	   << std::setw(12) << j.pz() ;
  if (j.has_user_info<FlavInfo>()) ostr << "  " << j.user_info<FlavInfo>().description();
  return ostr;
}


const FlavInfo MyFlavKtPlugin::_no_flav;

//----------------------------------------------------------------------
std::string MyFlavKtPlugin::description() const {
  ostringstream ostr;
  ostr << "Flavour-enabled Kt family Algorithm, imode = ";
  ostr << _imode;
  ostr << ", R = " << _R;
  return ostr.str();
}


//----------------------------------------------------------------------
void MyFlavKtPlugin::run_clustering(ClusterSequence & cs) const {

  int    np = cs.jets().size();
  double pp[4*np], p[4];
  int    nflav = 6;
  int    flav[nflav*np], beamflav[nflav*2];
  
  // lookup table: lookup[i] tells you what position i in the
  // nodelist corresponds to in our list of jets
  vector<int> lookup(2*np+1);

  vector<int> beam_particles;
  vector<int> spectators;
  bool found_flavour = false;
  
  vector<double> pt(np);
  

  int j = 0; // index in fortran array
  for (int i = 0; i < np; i++) {
    const PseudoJet & particle = cs.jets()[i];
    const FlavInfo & flavour = flavour_of(particle);
    found_flavour |= (!flavour.is_flavourless());
    if (flavour.is_beam()) {
      beam_particles.push_back(i); 
      continue;
    }

    pp[0+4*j] = particle.px();
    pp[1+4*j] = particle.py();
    pp[2+4*j] = particle.pz();
    pp[3+4*j] = particle.E() ;
    // set the fortran flav array based on the C++ flavour array
    for (int iflav = 0; iflav < nflav; iflav++) flav[iflav+nflav*j] = flavour[iflav+1];
    lookup[np+1+j] = i;
    j++;
  }
  // the lookup table started from np; but actually needs to start
  // from the "fortran" np, i.e. the number of particles passed to the
  // fortran code; so now shift things appropriately
  for (int i = 0; i < j; i++) {
    lookup[j+1+i] = lookup[np+1+i];
  }
  np = j; // np is now the number of "fortran" particles
  
  // set up the extras structure
  MyFlavKtPlugin::Extras * extras = new MyFlavKtPlugin::Extras(cs);
  // beam indices will run in parallel with the history (ideally, they would
  // be stored in the history) [NB: do not use "np" here, because it may have changed]
  extras->_beam_indices.reserve(2*cs.history().size());
  extras->_beam_indices.resize(cs.history().size(),0); 
  extras->_stepIsFlavourFriendly.reserve(2*cs.history().size());
  extras->_stepIsFlavourFriendly.resize(cs.history().size(),true);

  // now get the beam flavour
  BeamFlavPair beam_flavour;
  double yCMF;
  double z1;
  double z2;
  if (beam_particles.size() == 0) {
    //-----------------------------------
    // no beams: invent flavourless beams
    beam_flavour.backward = _no_flav;
    beam_flavour.forward  = _no_flav;
    if (spectators.size() != 0) throw Error("There were spectators but no beams; this is not allowed");
  } else if (beam_particles.size() == 2){
    //-----------------------------------
    // we have two beams
    // first entry should have smaller rapidity (i.e. negative pz)
    int b1 = beam_particles[0];
    int b2 = beam_particles[1];
    z1=cs.jets()[b1].pz();
    z2=cs.jets()[b2].pz();
    yCMF=(cs.jets()[b1]+cs.jets()[b2]).rap();
    NAMED_DEBUG("CMS_RAPIDITY",cout << " CMS rapidity: " << yCMF  << endl;)
    yCMF=0.5*log(z1/(-z2));
    NAMED_DEBUG("CMS_RAPIDITY",cout << " CMS rapidity: " << yCMF  << endl;)

    if (cs.jets()[b1].rap() > cs.jets()[b2].rap()) std::swap(b1,b2);
    beam_flavour.backward = flavour_of(cs.jets()[b1]);
    beam_flavour.forward  = flavour_of(cs.jets()[b2]);
    // cluster the beam particles with the "beam", simply to get rid of them.
    // b2 goes first to help correspondence with tiled-strategies ordering of plain kt.
    cs.plugin_record_iB_recombination(b2, 0.0); // dij = 0
    extras->_beam_indices.push_back(+1);
    extras->_stepIsFlavourFriendly.push_back(true);
    cs.plugin_record_iB_recombination(b1, 0.0); // dij = 0
    extras->_beam_indices.push_back(-1);
    extras->_stepIsFlavourFriendly.push_back(true);

  } else {
    //-----------------------------------
    ostringstream str;
    str << "Found " << beam_particles.size() << " beams; only legitimate values are 0 or 2";
    throw Error(str.str());
  }

  // transfer beam flavour info to f77 array
  for (int iflav = 0; iflav < nflav; iflav++) {
    beamflav[iflav         ] = beam_flavour.backward[iflav+1];
    beamflav[iflav + nflav ] = beam_flavour.forward [iflav+1];
  }
  

  beam_flavour.label_as_beam(); // always do this to make sure we have correct labelling
  extras->_beam_flavs.push_back(beam_flavour);
  vector<int> candidates;
  for (int ii=0;ii<np;ii++){
	  candidates.push_back(ii+2);
  }

  //bool doItMyself=false;
  bool doItMyself=true;
  int beam_index;
  bool friendly;
  bool doBoost;
  if(_nBoosts==1){
	  doBoost=true;
  } else {
	  doBoost=false;
  }

  for (int nclus=0;nclus<np;nclus++){
	  NAMED_DEBUG("BEAM_CLUSTERING",cout << "clustering round: " << nclus << endl;)
	  NAMED_DEBUG("BEAM_CLUSTERING", cout << "candidates at this point: " ;copy(candidates.begin(),candidates.end(),ostream_iterator<int>(cout," ")); cout << endl)
	  friendly = true;
	  double minDistance=1e32;
	  int parent1=-3;
	  int parent2=-3;
	  //at this point np is the number of jets
	  for (int candidateIndex=0;candidateIndex<candidates.size();candidateIndex++){
		  int jetIndex1=candidates[candidateIndex];
		  double diB=cs.jets()[jetIndex1].pt();
		  double diB2=diB*diB;
		  const PseudoJet& j1=cs.jets()[jetIndex1];
		  NAMED_DEBUG("BEAM_CLUSTERING",cout << "jet flavor " << flavour_of(j1).description() << endl;)
		  FlavInfo fb=beam_flavour.backward - flavour_of(j1);
		  FlavInfo ff=beam_flavour.forward - flavour_of(j1);

		  double yjet=j1.rap();
		  double pt=j1.pt();
		  double pt2=pt*pt;
		  double penaltyForward=1;
		  double penaltyBackward=1;
		  if (yjet>yCMF){
			  penaltyBackward*=100000;
		  } else {
			  penaltyForward*=100000;
		  }
		  if (ff.is_multiflavored() ){  // || (ff.is_flavorless() && !beam_flavour.forward.is_flavorless() )){   // it is ok to merge gluons into gluon
			  penaltyForward*=10000000;
		  }

		  if (fb.is_multiflavored() ){ // || (fb.is_flavorless() && !beam_flavour.backward.is_flavorless())){
			  penaltyBackward*=10000000;
		  } ;



		  NAMED_DEBUG("BEAM_CLUSTERING",cout << "considering forward merge: diB("<<jetIndex1<<") = " << diB << " penalty: " << penaltyForward <<endl;)
		  NAMED_DEBUG("BEAM_CLUSTERING",cout << "forward beam flavor " << beam_flavour.forward.description()	 << endl;)
		  NAMED_DEBUG("BEAM_CLUSTERING",cout << "projected flavor forward: " << ff.description() << endl;)
		  if (diB2*penaltyForward<minDistance && penaltyForward < penaltyBackward ){ // at least one of the two has the penalty for the wrong y direction
			  minDistance=diB2*penaltyForward; // I can do this because I don't use this as the scale
			  parent1=candidateIndex;
			  parent2=-1;
			  beam_index = 1;
		  }
		  NAMED_DEBUG("BEAM_CLUSTERING",cout << "considering backward diB("<<jetIndex1<<") = " << diB << " penalty: " << penaltyBackward <<endl;)
		  NAMED_DEBUG("BEAM_CLUSTERING",cout << "backward beam flavor " << beam_flavour.backward.description()	 << endl;)
		  NAMED_DEBUG("BEAM_CLUSTERING",cout << "projected flavor backward: " << fb.description() << endl;)
		  if (diB2*penaltyBackward<minDistance && penaltyForward > penaltyBackward){
			  minDistance=diB2*penaltyBackward;  // I can do this because I don't use this as the scale
			  parent1=candidateIndex;
			  parent2=-1;
			  beam_index =-1;
		  }

		  for (int candidateIndex2=candidateIndex+1;candidateIndex2<candidates.size();candidateIndex2++){
			  int index2=candidates[candidateIndex2];
			  double d;
			  FlavInfo fi= flavour_of(cs.jets()[jetIndex1])+ flavour_of(cs.jets()[index2]);
			  if (fi.is_multiflavored()){
				  d=1e32;
			  } else {
				  d=cs.jets()[jetIndex1].kt_distance(cs.jets()[index2]);
			  }
			  NAMED_DEBUG("BEAM_CLUSTERING",cout << "dij("<<jetIndex1<<","<< index2 << ") = " << sqrt(d) << endl;)
			  if (d<minDistance){
				  minDistance=d;
				  parent1=candidateIndex;
				  parent2=candidateIndex2;
			  }
		  }
	  }



	  for (int ii=2;ii<cs.jets().size();ii++){
		pt[ii]=cs.jets()[ii].pt();
	  }

	  if (parent2==-1){
		  int jetIndex1=candidates[parent1];
		  const PseudoJet& j1=cs.jets()[jetIndex1];
		  double yjet=j1.rap();
		  double pt=j1.pt();
		  double pt2=pt*pt;
		  if (beam_index==1){
			  FlavInfo ff=beam_flavour.forward - flavour_of(j1);
			  NAMED_DEBUG("BEAM_CLUSTERING",cout << "Flavor forward after this clustering: " << ff.description() << endl;)

			  if (yjet>yCMF){
				  NAMED_DEBUG("BEAM_CLUSTERING",cout << "Good, it is preferred to merge with forward jet. yjet=" << yjet << ">" << yCMF << "=yCMF"  << endl;)
				} else {
				  NAMED_DEBUG("BEAM_CLUSTERING",cout << "would have preferred to merge with backward jet. yjet=" << yjet << "<" << yCMF << "=yCMF"  << endl;)
						friendly=false;
				};

				cs.plugin_record_iB_recombination(jetIndex1, pt2);
				beam_flavour.forward  = beam_flavour.forward  - flavour_of(j1);
				beam_index = 1;
				if (ff.is_multiflavored() ||  (ff.is_flavorless() && !beam_flavour.forward.is_flavorless() )){
					NAMED_DEBUG("BEAM_CLUSTERING",cout << "This clustering is bad from the flavor point of view!..." << endl;)
					friendly=false;
				} else {
					NAMED_DEBUG("BEAM_CLUSTERING",cout << "Flavor is fine, I am happy with this clustering..." << endl;)
				}
		  } else {
			  FlavInfo fb=beam_flavour.backward - flavour_of(j1);
			  NAMED_DEBUG("BEAM_CLUSTERING",cout << "Flavor backward after this clustering: " << fb.description() << endl;)
				if (yjet<yCMF){
				  NAMED_DEBUG("BEAM_CLUSTERING",cout << "Good, it is preferred to merge with backward jet. yjet=" << yjet << "<" << yCMF << "=yCMF"  << endl;)
				} else {
				  NAMED_DEBUG("BEAM_CLUSTERING",cout << "Would have preferred to merge with forward jet. yjet=" << yjet << ">" << yCMF << "=yCMF"  << endl;)
							friendly=false;
				};
				cs.plugin_record_iB_recombination(jetIndex1, pt2);
				beam_flavour.backward  = beam_flavour.backward  - flavour_of(j1);
				beam_index = -1;
				if (fb.is_multiflavored() || (fb.is_flavorless() && !beam_flavour.backward.is_flavorless() ) ){
					NAMED_DEBUG("BEAM_CLUSTERING",cout << "This clustering is bad from the flavor point of view!..." << endl;)
					friendly=false;
				} else {
					NAMED_DEBUG("BEAM_CLUSTERING",cout << "Flavor is fine, I am happy with this clustering..." << endl;)
				}
		  }

		  if (beam_index==1){
			  z1=z1-j1.pz();
			  //E1=z1; only need z1
			  if (z1<0){ throw; } // this should not happen!
		  } else {
			  z2=z2-j1.pz();
			  //E2=-z2; only need z2
			  if (z2>0){ throw; } // this should not happen!
		  }

		yCMF=0.5*log(z1/(-z2));
		NAMED_DEBUG("CMS_RAPIDITY",cout << " CMS rapidity: " << yCMF  << " (before boost)" << endl;)

		  if (doBoost){

			  double betaz;
			  double keith_k0;
			  if (beam_index==1){
				  betaz=-j1.pz()/z1;
				  //E1=z1; only need z1
				  if (z1<0){ throw; } // this should not happen!
			  } else {
				  betaz=-j1.pz()/z2;
				  //E2=-z2; only need z2
				  if (z2>0){ throw; } // this should not happen!
			  }

	//	      k0rec=q0-cmpin(0,j)
	//	      krecv=-cmpin(1:3,j)
	//	      beta=-krecv(3)/k0rec
			  keith_k0=2*z1-j1.E(); //if we are in CMF z1=-z2
			  double keith_beta=j1.pz()/keith_k0;
			  double krec2=j1.px()*j1.px()+j1.py()*j1.py()+j1.pz()*j1.pz();
			  double krec=sqrt(krec2);
			  double recMass2=keith_k0*keith_k0-krec2;


			  double jperp=j1.perp();
			double norm=sqrt(j1.px()*j1.px()+j1.py()*j1.py()+j1.pz()*j1.pz());
			double Keith_norm=sqrt(recMass2+j1.perp2());

			double boostx=-j1.px()/Keith_norm;
			double boosty=-j1.py()/Keith_norm;



			NAMED_DEBUG("BOOST",
				std::cout << "keith_k0: " << keith_k0 << " keith_beta: " << keith_beta <<  std::endl ;
				std::cout << "betaz: " << betaz << " boostx: " << boostx <<" boosty: " << boosty << std::endl ;
				std::cout << "rebM2: " << recMass2 << std::endl ;
				std::cout << "keith norm: " << Keith_norm << std::endl ;
			)
			NAMED_DEBUG("BOOST",
			std::cout << setw(12) << "pt" << setw(12) << "rapidity" << setw(12) << "phi" << "  flav" << std::endl;
			)
			for (int ic=0;ic<candidates.size();ic++){
				int jetIndex=candidates[ic];
				const PseudoJet& jc=cs.jets()[jetIndex];
				NAMED_DEBUG("BOOST",std::cout << "Jet before boost: " << jc  << std::endl; );
				PseudoJet& j=const_cast<PseudoJet&>(jc);
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
				const PseudoJet& jc=cs.jets()[jetIndex];
				NAMED_DEBUG("BOOST",std::cout << "Jet before boost: " << jc  << std::endl; );
				PseudoJet& j=const_cast<PseudoJet&>(jc);
				TLorentzVector bm(j.px(),j.py(),j.pz(),j.E());
				NAMED_DEBUG("BOOST",std::cout << "Before boosts: " << bm  << std::endl; );
				bm.Boost(0,0,keith_beta);
				NAMED_DEBUG("BOOST",std::cout << "Jet after boost: " << cs.jets()[jetIndex]  << std::endl; );
			}


			doBoost=false;  //only run the boost once!
		    z1=cs.jets()[beam_particles[0]].pz();
		    z2=cs.jets()[beam_particles[1]].pz();
			yCMF=0.5*log(z1/(-z2));
			NAMED_DEBUG("CMS_RAPIDITY",cout << " CMS rapidity: " << yCMF  << " (after boost)" <<endl;)

		  }

		candidates.erase(candidates.begin()+parent1);


	  } else {
		  const PseudoJet &j1=cs.jets()[candidates[parent1]];
		  const PseudoJet &j2=cs.jets()[candidates[parent2]];
		  PseudoJet newjet(j1+j2);
		  // track the jet flavour
		  FlavInfo newflav =   flavour_of(j1)
							 + flavour_of(j2);
		  newjet.set_user_info(new FlavInfo(newflav));
		  int k;
		  NAMED_DEBUG("BEAM_CLUSTERING",cout << "combining jets " << candidates[parent1] << " and " << candidates[parent2] << endl;)
		  cs.plugin_record_ij_recombination(candidates[parent1],candidates[parent2], minDistance,   // for jet clustering there is no penatly factor so it's ok to use minDistance
											newjet, k);
		  candidates.erase(candidates.begin()+parent1);
		  candidates.erase(candidates.begin()+(parent2-1)); // now this element is one position more to the left. We also know that parent1< parent2
		  candidates.push_back(k);
		  NAMED_DEBUG("BEAM_CLUSTERING", cout << "candidates at this point: " ;copy(candidates.begin(),candidates.end(),ostream_iterator<int>(cout," ")); cout << endl)

	  }

		beam_flavour.label_as_beam(); // always do this to make sure we have correct labelling
		extras->_beam_flavs.push_back(beam_flavour);
		extras->_beam_indices.push_back(beam_index);
		extras->_stepIsFlavourFriendly.push_back(friendly);

	}

  extras->_bf_offset = int(extras->_beam_flavs.size()) - int(cs.history().size());
  auto_ptr<ClusterSequence::Extras> extras_auto_ptr(extras);
  cs.plugin_associate_extras(extras_auto_ptr);

}

//----------------------------------------------------------------------
/// return the index of _beam_flavs that corresponds to the beam
/// flavour at the stage where we have exactly n_exclusive_jets
int MyFlavKtPlugin::Extras::bf_index(int n_exclusive_jets) const {
  int index = 2*int(_cs.n_particles()) - 1 - n_exclusive_jets + _bf_offset;
  if (index < 0) index = 0;
  if (index >= int(_beam_flavs.size())) throw Error("beam flavours do not exist for the requested number of exclusive jets (maybe clustering stopped earlier?)");
  return index;
}

//----------------------------------------------------------------------
/// returns
/// - +1 if the specified jet clusters with the forward-going beam
/// - -1 if it clusters with the backward-going beam
/// -  0 if it clusters with another PseudoJet or does not cluster at all
int MyFlavKtPlugin::Extras::beam_it_clusters_with(const PseudoJet & jet) const {
  assert(jet.associated_cs() == &_cs); // make sure jets are right ones
  int child = _cs.history()[jet.cluster_hist_index()].child;
  if (child == ClusterSequence::Invalid) return 0;
  else {
    // if child index is outside this range, we have a problem
    assert (child >= 0 && child < int(_beam_indices.size()));
    return _beam_indices[child];
  }
    
}

bool MyFlavKtPlugin::Extras::beam_clustering_flavour_friendly(const PseudoJet & jet) const {
  assert(jet.associated_cs() == &_cs); // make sure jets are right ones
  int child = _cs.history()[jet.cluster_hist_index()].child;
  if (child == ClusterSequence::Invalid) throw;
  else {
    // if child index is outside this range, we have a problem
    assert (child >= 0 && child < int(_stepIsFlavourFriendly.size()));
    return _stepIsFlavourFriendly[child];
  }

}

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
