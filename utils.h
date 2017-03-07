
#include "fastjet/ClusterSequence.hh"
#include <vector>


int pdgFromFlavor(const fastjet::FlavInfo& fi);
bool hasQuarks(const fastjet::ClusterSequence& cs,int njets);
void displayClusterHistoryDot(fastjet::ClusterSequence& cs,std::ostream& os);

