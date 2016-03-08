#ifndef __FLAVKTPLUGIN_HH__
#define __FLAVKTPLUGIN_HH__
#include "fastjet/JetDefinition.hh"
#include "fastjet/LimitedWarning.hh"
#include "fastjet/ClusterSequence.hh"
#include "JetIModes.h"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh

//----------------------------------------------------------------------
/// class to allow representation of flavour, including concepts
/// such as the fact that a particle is an incoming "beam", or
/// should be a "spectator" during the clustering. 
///
/// The class also provides a facility to interpret PDG codes
/// and deduce the corresponding flavour content of the particle.
class FlavInfo : public PseudoJet::UserInfoBase {
public:
  /// constructs a flavour info from a pdg code. This can handle most
  /// standard model codes without difficulty. A zero code produces a flavourless object. 
  /// 
  /// The flags argument is optional, and can be set either to
  /// FlavInfo::beam or FlavInfo::spectator
  ///
  FlavInfo (int pdg_code = 0, int flags = 0);
  /// constructs a flavour info object from the individual flavours
  FlavInfo (int n_d, int n_u, int n_s, int n_c, int n_b, int n_t, int flags = 0);

  /// returns the net number of quarks of flavour iflv 
  /// (iflv runs from 1=d to 6=top).
  int operator[](int iflv) const {return _flav_content[iflv];}

  /// sets the number of objects of flavour iflv (1=d to 6=top) to be n
  void set_flav(int iflv, int n) {_flav_content[iflv] = n; update_flavourless_attribute();}

  /// resets all flavours to zero except iflv; this should be called,
  /// for example, when considering b-flavour at hadron level, so that
  /// the algorithm doesn't get confused by the many other quark
  /// flavours that are around (and also because, experimentally,
  /// those other flavours may not be well known).
  void reset_all_but_flav(int iflv);

  /// returns the pdg_code, or 0 if it's unknown (e.g. due to result
  /// of recombination)
  int pdg_code() const {return _pdg_code;}

  /// label this particle as being an incoming beam particle
  void label_as_beam() {_flav_content[0] |= beam;}
  /// returns true if this particle is a beam particle.
  bool is_beam() const {return (_flav_content[0] & beam);}

  /// label this object as being a "spectator", such as a W, which is 
  /// relevant for calculating the beam distance in flavour clusterings
  /// but does itself take part in the clustering
  void label_as_spectator() {_flav_content[0] |= spectator;}
  /// returns true if this particle is a spectator.
  bool is_spectator() const {return (_flav_content[0] & spectator);}

  /// returns true if the object has no net flavour
  bool is_flavourless() const {return (_flav_content[0] & _is_flavourless);}
  bool is_flavorless() const {return is_flavourless();}

  /// returns true if the object has more than one unit of flavour
  bool is_multiflavoured() const;
  bool is_multiflavored() const {return is_multiflavoured();}

  /// allows addition of flavour: note that beam, spectator and PDG status are lost
  FlavInfo operator+(const FlavInfo &) const;
  /// allows subtraction of flavour: note that beam, spectator and PDG status are lost
  FlavInfo operator-(const FlavInfo &) const;

  /// returns a string such as "u d", "cbar", etc.; "g" means gluon or anything
  /// else with no flavour
  std::string description() const;

  /// returns the flavour of a particle if that particle has flavour, otherwise
  /// just the default flavour
  static const FlavInfo & flavour_of(const PseudoJet & particle) {
    if (particle.has_user_info<FlavInfo>()) return particle.user_info<FlavInfo>();
    else                                    return _no_flav;
  }

  /// value of flag to indicate that the particle is an incoming beam particle
  static const int beam = 2;
  /// value of flag to indicate that the particle is a "spectator",
  /// such as a W, which is relevant for calculating the beam distance
  /// in flavour clusterings but does itself take part in the
  /// clustering
  static const int spectator = 4; 

private:
  int _pdg_code; 
  int _flav_content[7];
  static const int _is_flavourless = 1; // 

  static const FlavInfo _no_flav;
  void update_flavourless_attribute();
};


//----------------------------------------------------------------------
/// A class that carries out flavour recombination (and E-scheme
/// momentum recombination), to facilitate tracking of flavour
/// in the case of normal FastJet algorithms.
class FlavRecombiner : public JetDefinition::DefaultRecombiner {
public:
  FlavRecombiner() : DefaultRecombiner() {}
  std::string description() const {return DefaultRecombiner::description() + " and flavour recombination ";}
  void recombine(const PseudoJet & pa, const PseudoJet & pb, 
                 PseudoJet & pab) const  {
    DefaultRecombiner::recombine(pa, pb, pab);
    pab.set_user_info(new FlavInfo(FlavInfo::flavour_of(pa) + FlavInfo::flavour_of(pb)));
  }
};




//======================================================================
/// A plugin that provides a FastJet interface to the Kt family
/// flavour f90 clustering code by Banfi, Salam and Zanderighi
/// (Eur.Phys.J. C47 (2006) 113 [hep-ph/0601139]).
///
/// It provides both plain kt clustering (while tracking flavours) and
/// flavour-kt clustering. Access to beam flavour is possible via a
/// FlavKtPlugin::Extras class.
///
/// For flavour to be properly handled, the particles passed to
/// FastJet for clustering should each have a FlavInfo object set as
/// their user_info. Incoming beam particles should be included in the
/// list of input particles, and their flavour should also have the
/// "beam" attribute set.
///
/// If the user asks for inclusive jets, they will receive the
/// standard set of inclusive jets, as well as the two beam particles.
///
/// If the user asks for exclusive jets, they may also obtain the
/// flavour of the incoming beams at that exclusive stage of
/// clustering via the beam_flav_forward(n) etc members of the
/// extras. See the example for usage.
class FlavKtPlugin : public JetDefinition::Plugin {
public:
  /// Main constructor for the FlavKt Plugin class.  It takes an "imode", which 
  /// indicates the kind of algorithm to be used, and a radius value (which
  /// will be ignored for e+e- algorithms). 
  /// 
  /// IMode values likely to be of interest include
  ///
  ///   - ialg_pp_ktLI_E:  standard pp longitudinally invariant kt 
  ///   - ialg_pp_ktf1_E:  pp flavour kt with alpha=1
  ///   - ialg_pp_ktf2_E:  pp flavour kt with alpha=2
  ///   - ialg_pp_ktfs1_E: pp flavour kt with alpha=1 and only "safe" flavour recombinations (e.g. no u+d->ud)
  ///   - ialg_pp_ktfs2_E: pp flavour kt with alpha=2 and only "safe" flavour recombinations (e.g. no u+d->ud)
  ///
  /// Consult the JetIModes.h file for the full available list of
  /// imode values. 
  ///
  /// Note that with the "safe" flavour variants, the jet algorithm
  /// may reach a stage where the only possible clustering is an
  /// unsafe one. This is most often the case in the vicinity of the
  /// beam remnants, e.g. proton -> neutron + pi+ -> neutron + u +
  /// dbar, where the dbar goes into the hard reaction, while the
  /// neutron and u go down the beampipe. There is no way of combining
  /// the neutron and u quark to get something with just a single unit
  /// of flavour.
  /// 
  /// In such cases, the algorithm will then with these unsafe
  /// clusterings at the very end, assigning a dmerge value of
  /// 1e100. This can lead to unexpected results if you then ask (say)
  /// for 2 exclusive jets.
  FlavKtPlugin (int imode, double R) : _imode(imode), _R(R) {}

  /// standard plugin call to run the clustering
  virtual void run_clustering(ClusterSequence & cs) const;
  /// returns a description of the plugin
  virtual std::string description() const;

  /// returns the R value used.
  virtual double R() const {return _R;}

  /// returns true, which tells the FJ that the exclusive sequence can be used here.
  virtual bool exclusive_sequence_meaningful() const {return true;}

  class Extras;
  
private:
  int    _imode;
  double _R;
  const FlavInfo & flavour_of(const PseudoJet & particle) const {
    if (particle.has_user_info<FlavInfo>()) return particle.user_info<FlavInfo>();
    else                                    return _no_flav;
  }
  static const FlavInfo _no_flav;
  static LimitedWarning _spectators_unhandled;

  struct BeamFlavPair {
    FlavInfo forward, backward;
    void label_as_beam() {forward.label_as_beam(); backward.label_as_beam();}
  };
};




//----------------------------------------------------------------------
class FlavKtPlugin::Extras : public ClusterSequence::Extras {
public:

  /// returns the flavour of the forward-going (more positive rapidity) beam
  /// at the stage of the clustering with n_exclusive_jets
  const FlavInfo & beam_flav_forward (int n_exclusive_jets) const {
    return _beam_flavs[bf_index(n_exclusive_jets)].forward;
  }

  /// returns the flavour of the backward-going (more negative rapidity) beam
  /// at the stage of the clustering with n_exclusive_jets
  const FlavInfo & beam_flav_backward(int n_exclusive_jets) const {
    return _beam_flavs[bf_index(n_exclusive_jets)].backward;
  }

  /// returns
  /// - +1 if the specified jet clusters with the forward-going beam
  /// - -1 if it clusters with the backward-going beam
  /// -  0 if it clusters with another PseudoJet or does not cluster at all
  int beam_it_clusters_with(const PseudoJet & jet) const;

  virtual std::string description() const {return "the extras class for FlavKtPlugin";}  

protected:
  // constructor is not intended for public use...
  Extras(const ClusterSequence & cs) : _cs(cs) {}

  std::vector<FlavKtPlugin::BeamFlavPair> _beam_flavs;
  std::vector<int>                        _beam_indices;
  const ClusterSequence & _cs;
  int _bf_offset;
  friend class FlavKtPlugin;

  /// return the index of _beam_flavs that corresponds to the beam
  /// flavour at the stage where we have exactly n_exclusive_jets
  int bf_index(int n_exclusive_jets) const;
};


FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
#endif // __FLAVKTPLUGIN_HH__

