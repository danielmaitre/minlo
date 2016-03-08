#include "FlavKtPlugin.hh"
#include "fastjet/ClusterSequence.hh"
#include "ktflavf77.h"
#include <sstream>
#include <utility>

using namespace std;

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh


//----------------------------------------------------------------------
FlavInfo::FlavInfo(int n_d, int n_u, int n_s, int n_c, int n_b, int n_t, int flags) : 
  _pdg_code(0) {
  _flav_content[0] = flags;
  _flav_content[1] = n_d;
  _flav_content[2] = n_u;
  _flav_content[3] = n_s;
  _flav_content[4] = n_c;
  _flav_content[5] = n_b;
  _flav_content[6] = n_t;
  update_flavourless_attribute();
}

//----------------------------------------------------------------------
FlavInfo::FlavInfo(int pdg_code, int flags) : _pdg_code(pdg_code){
  _flav_content[0] = flags;
  for(unsigned i = 1; i <= 6; i++) _flav_content[i] = 0;

  // for particles with illicit (zero) pdg_code, no work to be done
  if (_pdg_code == 0) return;

  int netsign = (pdg_code >= 0 ? +1 : -1);
  pdg_code = abs(pdg_code);

  // extract digits of the pdg_code, since these contain information 
  // on flavour of component quarks
  valarray<int> digit(4);
  int           ndigits = 0;
  for (int i = 0; i < 4; i++) {
    digit[i] = pdg_code % 10;
    if (digit[i] != 0) ndigits = i+1;
    pdg_code /= 10; // "shift" things along
  }
  
  // start this part with _flav_content already initialised to zero
  // in constructor
  if (ndigits == 1) { // a lone quark
    if (digit[0] > 6 || digit[0] == 0) {
      cerr << "FlavInfo failed to understand pdg_code = "<<_pdg_code<<endl; exit(-1);}
    _flav_content[digit[0]] = netsign;

  } else if (ndigits == 2) { // a lepton, photon or cluster [flav lost...]
    // do nothing...

  } else { // must be a meson, cluster or baryon
    // check sanity of codes
    for (int i=1; i < ndigits; i++) {
      if (digit[i] > 6) {cerr << "FlavInfo failed to understand pdg_code = "
			       <<_pdg_code<<endl; exit(-1);}}
    
    // now deal with different cases
    if (ndigits == 4) { // diquark [nm0x] or baryon [nmpx]
      for (int i=1; i < ndigits; i++) {
	if (digit[i] > 0) _flav_content[digit[i]] += netsign;}
    } else if (ndigits == 3) { // meson [nmx]
      // Beware of PDG convention that says that a K+ or B+ are a
      // particle and so have positive pdg_code (i.e. flavcodes > 1). So
      if (digit[2] == 3 || digit[2] == 5) netsign = -netsign;
      _flav_content[digit[2]] += netsign;
      _flav_content[digit[1]] -= netsign;
    } else {
      cerr << "FlavInfo failed to understand pdg_code = " <<_pdg_code<<endl; exit(-1);}
  }
  update_flavourless_attribute();
}


//----------------------------------------------------------------------
void FlavInfo::reset_all_but_flav(int iflv) {
  for (int i = 1; i <= 6; i++) {
    if (i != iflv) _flav_content[i] = 0;
  }
  update_flavourless_attribute();
}

//----------------------------------------------------------------------
FlavInfo FlavInfo::operator+(const FlavInfo & other) const {
  FlavInfo sum(operator[](1)+other[1],
               operator[](2)+other[2],
               operator[](3)+other[3],
               operator[](4)+other[4],
               operator[](5)+other[5],
               operator[](6)+other[6]
               );
  sum.update_flavourless_attribute();
  return sum;
}
//----------------------------------------------------------------------
FlavInfo FlavInfo::operator-(const FlavInfo & other) const {
  FlavInfo sum(operator[](1)-other[1],
               operator[](2)-other[2],
               operator[](3)-other[3],
               operator[](4)-other[4],
               operator[](5)-other[5],
               operator[](6)-other[6]
               );
  sum.update_flavourless_attribute();
  return sum;
}

//----------------------------------------------------------------------
void FlavInfo::update_flavourless_attribute() {
  for (unsigned i = 1; i <=6; i++) {
    if (_flav_content[i] != 0) {
      _flav_content[0] &= ~ _is_flavourless; // ~ is C++ bitwise not
      return;
    }
  }
  _flav_content[0] |= _is_flavourless;
}

//----------------------------------------------------------------------
bool FlavInfo::is_multiflavoured() const {
  int flavsum = 0;
  for (unsigned i = 1; i <=6; i++) flavsum += abs(_flav_content[i]);  
  return flavsum > 1;
}


//----------------------------------------------------------------------
// an object with no flavour, that we can conveniently point to
const FlavInfo FlavInfo::_no_flav;
const FlavInfo FlavKtPlugin::_no_flav;
LimitedWarning FlavKtPlugin::_spectators_unhandled;

//----------------------------------------------------------------------
string FlavInfo::description() const {
  const char * flavs = "duscbt";
  ostringstream result;
  
  if (is_flavourless()) {
    result << "g "; 
  } else {
    for (int iflav = 1; iflav <= 6; iflav++) {
      int n = operator[](iflav);
      for (unsigned i = 0; i < abs(n); i++) {
        result << flavs[iflav-1];
        if (n<0) result << "bar";
        result << " ";
      }
    }
  }
  if (is_beam()) result << "(beam) ";
  if (is_spectator()) result << "(spectator) ";
  return result.str();
}

//----------------------------------------------------------------------
std::string FlavKtPlugin::description() const {
  ostringstream ostr;
  ostr << "Flavour-enabled Kt family Algorithm, imode = ";
  ostr << _imode;
  ostr << ", R = " << _R;
  return ostr.str();
}


//----------------------------------------------------------------------
void FlavKtPlugin::run_clustering(ClusterSequence & cs) const {

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
  
  
  int j = 0; // index in fortran array
  for (int i = 0; i < np; i++) {
    const PseudoJet & particle = cs.jets()[i];
    const FlavInfo & flavour = flavour_of(particle);
    found_flavour |= (!flavour.is_flavourless());
    if (flavour.is_beam()) {
      beam_particles.push_back(i); 
      continue;
    }
    if (flavour.is_spectator()) {
      // do we check spectators are flavourless?
      if (!flavour.is_flavourless()) throw Error("FlavKtPlugin: spectators should be flavourless");
      spectators.push_back(i);
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
  FlavKtPlugin::Extras * extras = new FlavKtPlugin::Extras(cs);
  // beam indices will run in parallel with the history (ideally, they would
  // be stored in the history) [NB: do not use "np" here, because it may have changed]
  extras->_beam_indices.reserve(2*cs.history().size());
  extras->_beam_indices.resize(cs.history().size(),0); 

  // now get the beam flavour
  BeamFlavPair beam_flavour;
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
    if (cs.jets()[b1].rap() > cs.jets()[b2].rap()) std::swap(b1,b2);
    beam_flavour.backward = flavour_of(cs.jets()[b1]);
    beam_flavour.forward  = flavour_of(cs.jets()[b2]);
    if (spectators.size() != 0) {
      _spectators_unhandled.warn("Found spectators: but these will not be included in the beam distance (and are thus being treated wrongly)");
    }
    // cluster the beam particles with the "beam", simply to get rid of them.
    // b2 goes first to help correspondence with tiled-strategies ordering of plain kt.
    cs.plugin_record_iB_recombination(b2, 0.0); // dij = 0
    extras->_beam_indices.push_back(+1);
    cs.plugin_record_iB_recombination(b1, 0.0); // dij = 0
    extras->_beam_indices.push_back(-1);
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
  
  // the f90 code via its f77 interface
  ktf77MakeNodeListFlav(_imode, _R, np, pp, nflav, flav, beamflav);

  beam_flavour.label_as_beam(); // always do this to make sure we have correct labelling
  extras->_beam_flavs.push_back(beam_flavour);


  // decant the information from the f90 node list
  int    parent1, parent2;
  double dij;

  int beam_index;
  for (int i = np; i > 0; i--) {
    ktf77ReadNodeList(i, p, parent1, parent2, dij);
    
    if (parent1 > 0) {
      PseudoJet newjet(p);
      // track the jet flavour
      FlavInfo newflav =   flavour_of(cs.jets()[lookup[parent1]])
                         + flavour_of(cs.jets()[lookup[parent2]]);
      newjet.set_user_info(new FlavInfo(newflav));
      int k;
      //cout << i << " " << parent1 << " " << parent2 << " " << dij << endl;
      cs.plugin_record_ij_recombination(lookup[parent1], lookup[parent2], dij, 
                                        newjet, k);
      lookup[i] = k;
      beam_index = 0;
      //cout << "combining " << lookup[parent1] << " " << lookup[parent2] << " into " << k << " " << flavour_of(newjet).description() << endl;
    } else {
      cs.plugin_record_iB_recombination(lookup[parent2], dij);
      if        (parent1 == -2) {
        beam_flavour.backward = beam_flavour.backward - flavour_of(cs.jets()[lookup[parent2]]);
        beam_index = -1;
      } else if (parent1 == -1) {
        beam_flavour.forward  = beam_flavour.forward  - flavour_of(cs.jets()[lookup[parent2]]);
        beam_index = +1;
      } else {
        throw Error("inconsistent parent1 value");
      }
    }
    beam_flavour.label_as_beam(); // always do this to make sure we have correct labelling
    extras->_beam_flavs.push_back(beam_flavour);
    extras->_beam_indices.push_back(beam_index);
  }
  ktf77DeleteNodeList();
  extras->_bf_offset = int(extras->_beam_flavs.size()) - int(cs.history().size());
  auto_ptr<ClusterSequence::Extras> extras_auto_ptr(extras);
  cs.plugin_associate_extras(extras_auto_ptr);
}

//----------------------------------------------------------------------
/// return the index of _beam_flavs that corresponds to the beam
/// flavour at the stage where we have exactly n_exclusive_jets
int FlavKtPlugin::Extras::bf_index(int n_exclusive_jets) const {
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
int FlavKtPlugin::Extras::beam_it_clusters_with(const PseudoJet & jet) const {
  assert(jet.associated_cs() == &_cs); // make sure jets are right ones
  int child = _cs.history()[jet.cluster_hist_index()].child;
  if (child == ClusterSequence::Invalid) return 0;
  else {
    // if child index is outside this range, we have a problem
    assert (child >= 0 && child < int(_beam_indices.size()));
    return _beam_indices[child];
  }
    
}

FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
