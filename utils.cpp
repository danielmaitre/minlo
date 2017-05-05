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

typedef fastjet::MyFlavKtPlugin THEPLUGIN;


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

std::string getStyle(const std::string& desc){
				if (desc=="g (beam) "){
					return "style=dotted";
				} else {
					return "style=solid";
				};
}


std::string getStyle(int pdg){
	if (pdg==21){
		return "style=dotted";
	} else {
		return "style=solid";
	};
}

std::string getStyle(const fastjet::ClusterSequence& cs,int historyIndex){
	int pdg=pdgFromFlavor(cs.jets()[cs.history()[historyIndex].jetp_index].user_info<fastjet::FlavInfo>());
	return getStyle(pdg);
}
std::string getStyleParent(const fastjet::ClusterSequence& cs,int historyIndex){
	int parent=cs.history()[historyIndex].parent1;
	int jetp=cs.history()[parent].jetp_index;
	int pdg=pdgFromFlavor(cs.jets()[jetp].user_info<fastjet::FlavInfo>());
	return getStyle(pdg);
}

std::string getStyleFromIndex(const fastjet::ClusterSequence& cs,int index){
	int pdg=pdgFromFlavor(cs.jets()[index].user_info<fastjet::FlavInfo>());
	return getStyle(pdg);
}


void displayClusterHistoryDot(fastjet::ClusterSequence& cs,std::ostream& os){
	const std::vector<fastjet::ClusterSequence::history_element> & history = 			cs.history();
	int n = history.size();

	const THEPLUGIN::Extras * extras =
			dynamic_cast<const THEPLUGIN::Extras *>(cs.extras());

	int forward=0;
	int backward=1;

	int njetsCurrent=(n/2)-2 ; // change for something else than Z+4j!!!!

	std::vector<int> forwardNodes;
	std::vector<int> backwardNodes;
		os << "graph test {" << std::endl;

	os << 0 <<"[style=filled]"<<std::endl;
	os << n <<"[style=filled,color=red]"<<std::endl;
	os << 1 <<"[style=filled]"<<std::endl;

	for (int historyIndex=(n/2 + 2);historyIndex<n;historyIndex++){

	    int parent1=history[historyIndex].parent1;
	    int parent2=history[historyIndex].parent2;

		int maxHist=n;

	    if (parent2 == cs.BeamJet){
	    	fastjet::PseudoJet j=cs.jets()[history[parent1].jetp_index];

			if (extras->beam_it_clusters_with(j)==+1){
				os << "// cluster with forward beam " << pdgFromFlavor(extras->beam_flav_forward(njetsCurrent))<<std::endl;
		        os << "// forward beam before: '" << extras->beam_flav_forward(njetsCurrent+1).description() <<"'"<< std::endl;
			    os << "// forward beam now:" << extras->beam_flav_forward(njetsCurrent).description() << "(njetsCurrent=" << njetsCurrent << ")" << std::endl;
				std::string beamstyle=getStyle(extras->beam_flav_forward(njetsCurrent).description());
				std::string otherstyle=getStyleParent(cs,historyIndex);
				os << forward << " -- " << historyIndex << " [" << beamstyle <<",weight=0.5]" << std::endl;
				os << parent1 << " -- " << historyIndex << " ["<< otherstyle << "]"<< std::endl;
				forward=historyIndex;
				forwardNodes.push_back(historyIndex);
	        } else {
				os << "// cluster with backward beam " << pdgFromFlavor(extras->beam_flav_backward(njetsCurrent))<<std::endl;
		        os << "// backward beam before:" << extras->beam_flav_backward(njetsCurrent+1).description() << std::endl;
			    os << "// backward beam now:" << extras->beam_flav_backward(njetsCurrent).description() << std::endl;
				std::string beamstyle=getStyle(pdgFromFlavor(extras->beam_flav_backward(njetsCurrent)));
				std::string otherstyle=getStyleParent(cs,historyIndex);
				os << historyIndex  << " -- " << backward<< " [" << beamstyle  <<",weight=0.5]"<< std::endl;
				os << parent1 << " -- " << historyIndex << " ["<< otherstyle << "]" << std::endl;
				backward=historyIndex;
				backwardNodes.push_back(historyIndex);
			}
	    } else if ( parent1 >=0 and parent2 >= 0) {
			os << "// jet clustered together: " << parent1 << " " << parent2<< std::endl ;
			std::string style1=getStyle(pdgFromFlavor(cs.jets()[history[parent1].jetp_index].user_info<fastjet::FlavInfo>()));
			std::string style2=getStyle(pdgFromFlavor(cs.jets()[history[parent2].jetp_index].user_info<fastjet::FlavInfo>()));
    		os << parent1 << " -- " << historyIndex << " [" << style1 << "]" << std::endl;
    		os << parent2 << " -- " << historyIndex << " [" << style2 << "]" <<  std::endl;
    	}
		os << historyIndex << "[ label= \""<< historyIndex<<" \\n " << cs.history()[historyIndex].dij << "\" ]" << std::endl;
		njetsCurrent--;
	}
	os << "// connecting to the main process: " << std::endl ;

			std::string beamstyleF=getStyle(extras->beam_flav_forward(0+1).description());
			std::string beamstyleB=getStyle(extras->beam_flav_backward(0+1).description());
			os << forward << " -- " << n << " [" << beamstyleF << ",weight=0.5]" << std::endl;
			os << n << " -- " << backward << " ["<< beamstyleB << ",weight=0.5]" << std::endl;

    		os << "{ rank=same; 0;";
			for (int ii=0;ii<forwardNodes.size();ii++){
				os << forwardNodes[ii] << ";";
			}
			os << n << ";";
			for (int ii=backwardNodes.size()-1;ii>=0;ii--){
				os << backwardNodes[ii] << ";";
			}
			os << "1;}  \n}" << std::endl;

}
