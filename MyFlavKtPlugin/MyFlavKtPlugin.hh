#ifndef __MYFLAVKTPLUGIN_HH__
#define __MYFLAVKTPLUGIN_HH__
#include "fastjet/JetDefinition.hh"
#include "fastjet/LimitedWarning.hh"
#include "fastjet/ClusterSequence.hh"
#include "JetIModes.h"
#include "FlavKtPlugin.hh"

namespace fastjet {     // defined in fastjet/internal/base.hh



class MyFlavKtPlugin : public JetDefinition::Plugin {
public:
  /// based on the FlavKtPlugin
  MyFlavKtPlugin (double R,bool useModifiedR,int nb,bool useRapidity,bool doSimpleBoost) :  _R(R), _useModifiedR(useModifiedR), _nBoosts(nb), _useRapidity(useRapidity) , _doSimpleBoost(doSimpleBoost) {}

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
  double _R;
  bool _useModifiedR;
  bool _useRapidity;
  bool _doSimpleBoost;
  int _nBoosts;
  static const FlavInfo _no_flav;
  static LimitedWarning _spectators_unhandled;

  struct BeamFlavPair {
    FlavInfo forward, backward;
    void label_as_beam() {forward.label_as_beam(); backward.label_as_beam();}
  };

public:
  static const FlavInfo & flavour_of(const PseudoJet & particle) {
    if (particle.has_user_info<FlavInfo>()) return particle.user_info<FlavInfo>();
    else                                    return _no_flav;
  }


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

  void setFSRboost(double x, double y,double z){_hasFSRboost=true;  _FSRboost[0]=x,_FSRboost[1]=y;_FSRboost[2]=z; }
  void getFSRboost(double* vec) const { vec[0]=_FSRboost[0]; vec[1]=_FSRboost[1]; vec[2]=_FSRboost[2]; }

  void setISRxyboost(double x, double y,double z){_hasISRboost=true;  _ISRxyboost[0]=x,_ISRxyboost[1]=y;_ISRxyboost[2]=z; }
  void getISRxyboost(double* vec) const { vec[0]=_ISRxyboost[0]; vec[1]=_ISRxyboost[1]; vec[2]=_ISRxyboost[2]; }
  void setISRzboost(double x, double y,double z){_hasISRboost=true;  _ISRzboost[0]=x,_ISRzboost[1]=y;_ISRzboost[2]=z; }
  void getISRzboost(double* vec) const { vec[0]=_ISRzboost[0]; vec[1]=_ISRzboost[1]; vec[2]=_ISRzboost[2]; }

  bool hasFSRboost() const { return _hasFSRboost;}
  bool hasISRboost() const { return _hasISRboost;}


protected:
  // constructor is not intended for public use...
  Extras(const ClusterSequence & cs) : _cs(cs) {
  _hasISRboost=false;
  _hasFSRboost=false;
  }

  std::vector<MyFlavKtPlugin::BeamFlavPair> _beam_flavs;
  std::vector<int>                        _beam_indices;
  std::vector<bool>                        _stepIsFlavourFriendly;
  const ClusterSequence & _cs;
  int _bf_offset;
  friend class MyFlavKtPlugin;

  /// return the index of _beam_flavs that corresponds to the beam
  /// flavour at the stage where we have exactly n_exclusive_jets
  int bf_index(int n_exclusive_jets) const;

  bool _hasISRboost;
  bool _hasFSRboost;
  double _ISRzboost[3];
  double _ISRxyboost[3];
  double _FSRboost[3];
};


}       // defined in fastjet/internal/base.hh
#endif // __MYFLAVKTPLUGIN_HH__

