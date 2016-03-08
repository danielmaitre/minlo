#ifndef __MYFLAVKTPLUGIN_HH__
#define __MYFLAVKTPLUGIN_HH__
#include "fastjet/JetDefinition.hh"
#include "fastjet/LimitedWarning.hh"
#include "fastjet/ClusterSequence.hh"
#include "JetIModes.h"
#include "FlavKtPlugin.hh"

FASTJET_BEGIN_NAMESPACE      // defined in fastjet/internal/base.hh




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
class MyFlavKtPlugin : public JetDefinition::Plugin {
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
  MyFlavKtPlugin (int imode, double R,int nb) : _imode(imode), _R(R), _nBoosts(nb) {}

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
  int _nBoosts;
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
class MyFlavKtPlugin::Extras : public ClusterSequence::Extras {
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
  bool beam_clustering_flavour_friendly(const PseudoJet & jet) const;


  virtual std::string description() const {return "the extras class for FlavKtPlugin";}  

protected:
  // constructor is not intended for public use...
  Extras(const ClusterSequence & cs) : _cs(cs) {}

  std::vector<MyFlavKtPlugin::BeamFlavPair> _beam_flavs;
  std::vector<int>                        _beam_indices;
  std::vector<bool>                        _stepIsFlavourFriendly;
  const ClusterSequence & _cs;
  int _bf_offset;
  friend class MyFlavKtPlugin;

  /// return the index of _beam_flavs that corresponds to the beam
  /// flavour at the stage where we have exactly n_exclusive_jets
  int bf_index(int n_exclusive_jets) const;
};


FASTJET_END_NAMESPACE        // defined in fastjet/internal/base.hh
#endif // __MYFLAVKTPLUGIN_HH__

