/*
 * MINLOfunctions.cpp
 *
 *  Created on: 28 Feb 2016
 *      Author: daniel
 */

#include "ntuplereader/NtupleInfo.h"
#include "MINLOfunctions.h"
#include "debug.h"
#include <vector>
#include "fastjet/Selector.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"

#include "MyFlavKtPlugin/MyFlavKtPlugin.hh"
#include "MinloKeithDeclarations.h"
#include "LHAPDF/LHAPDF.h"

#include "InterpolatedFunction.h"
#include "pdf.h"
#include "utils.h"
#include "nll.h"

using namespace std;

#define myflavorPlugin

typedef fastjet::MyFlavKtPlugin THEPLUGIN;

double minloImpl::g_MinloScale;
double minloImpl::g_MinloFactor;
double minloImpl::g_nclusterings;
bool minloImpl::g_MinloScaleValid;

bool minloDot::s_MinloWriteDot=false;
std::ostream *minloDot::s_stream=0;


const double fixedScaleForNLO=80.419;


// this is to make sure scales are raising all the way
// not just at the last step
// false is the behaviour of indie-minlo
//bool raisingAllTheWay=true;
// if set to true then the local scale cannot be smaller than the largest clustering scale
bool raisingScales=true;

double alphasScaleMin=1;
double minRenScale=1;

#include "MinloKeithDeclarations.h"
// can't be const so they can be passed to fortran
double      KEITH_rad_charmthr2=KEITH_rad_charmthr2_CST;   //! Charm  (mass) threshold squared: (1.5 GeV)^2
double      KEITH_rad_bottomthr2=KEITH_rad_bottomthr2_CST;  //! Bottom (mass) threshold squared: (5.0 GeV)^2
int         KEITH_st_nlight=KEITH_st_nlight_CST;       //! Nominal number of light flavs in computation (5)
double      KEITH_st_lambda5MSB=KEITH_st_lambda5MSB_CST; /*same as in sudakov.f*/;   //! Lambda_5 in the MS bar scheme


//typedef FlavKtPlugin THEPLUGIN;
typedef fastjet::MyFlavKtPlugin THEPLUGIN;



// helpers functions to print pseudojets
ostream & operator<<(ostream& ostr, const fastjet::PseudoJet & j) {
  ostr << setw(12) << j.perp()
       << setw(12) << j.rap()
       << setw(12) << j.phi();
  if (j.has_user_info<fastjet::FlavInfo>()) ostr << "  " << j.user_info<fastjet::FlavInfo>().description();
  return ostr;
}
ostream& print_header(ostream& ostr) {
  ostr << setw(12) << "pt"
       << setw(12) << "rapidity"
       << setw(12) << "phi";
  ostr << "  flav";
  return ostr;
}

std::ostream & operator<<(std::ostream& ostr, const sudakovCandidate & sc){
	return ostr << "high: " << sc.highScale << " low: " << sc.lowScale << " falv: " << sc.flavor << " hist: " << sc.historyIndex ;
} ;


double getAlphasQ(double Q,const MinloInfo& MI){
	double res;
	if (Q<alphasScaleMin){
		res=getAlphasQ(alphasScaleMin,MI);
	} else {
		if (!MI.d_usePDFalphas){
			if (MI.d_useLOalphas){
				double Q2=Q*Q;
				double l2=MI.d_lambda*MI.d_lambda;
				int nf=5;
				res = I_PWHG_ALPHAS0(Q2,l2,nf);
			} else {
				double Q2=Q*Q;
				double l=MI.d_lambda;
				res = I_PWHG_ALPHAS(Q2,l,KEITH_st_nlight, KEITH_rad_charmthr2, KEITH_rad_bottomthr2);
			}
		} else {
		#ifdef LHAPDF_NEW_VERSION
			res = currentPDF::s_PDF->alphasQ(Q);
		#else
			res = LHAPDF::alphasPDF(Q);
		#endif
		}
	}
	return res;
}

double getAlphasQ2(double Q2,bool useKeithAlphas){
	if (useKeithAlphas){
		return I_PWHG_ALPHAS(Q2,KEITH_st_lambda5MSB,KEITH_st_nlight, KEITH_rad_charmthr2, KEITH_rad_bottomthr2);
	} else {
	#ifdef LHAPDF_NEW_VERSION
		return currentPDF::s_PDF->alphasQ2(Q2);
	#else
		double q=sqrt(Q2);
		return LHAPDF::alphasPDF(q);
	#endif

	}
}




//this assumes that the first two particles are leptons and not partons to be clustered.
// it should be made more flexible/smart

void fillJetVector(NtupleInfo<MAX_NBR_PARTICLES>& Ev,std::vector<fastjet::PseudoJet>& particles,const int* flavors,bool useFlavor,bool useDouble){
	for (int p=2;p<Ev.nparticle;p++){
		if (useDouble){
			particles.push_back(fastjet::PseudoJet( Ev.pxD[p],Ev.pyD[p],Ev.pzD[p],Ev.ED[p] ));
		} else {
			particles.push_back(fastjet::PseudoJet( Ev.px[p],Ev.py[p],Ev.pz[p],Ev.E[p] ));
		}
		if (!useFlavor){
			particles.back().set_user_index(flavors[p]);
		} else {
			#if defined(flavorPlugin) or defined(myflavorPlugin)
					particles.back().set_user_info(new fastjet::FlavInfo(flavors[p]));
			#endif
		}
	}
}



void displayClusterHistory(fastjet::ClusterSequence& cs){
	const std::vector<fastjet::ClusterSequence::history_element> & history = 			cs.history();
	int n = history.size();
	for (int i=0;i<n;i++){
		std::cout << i << ": parents: " << history[i].parent1 << "," << history[i].parent2 << std::endl;
		std::cout << " child: " << history[i].child << "  dij: " << history[i].dij << std::endl;
	}
}


void computeSudakov(const MinloInfo& MI,double highScale,double lowScale,double q02,int flavor,double& sudakov, double& bornSub) {
	if (highScale > lowScale) {
		//sudakov=lo_sudakov(q02,highScale,lowScale,f % 21 );
		sudakov = nll_sudakov_withInfo(MI,q02, highScale, lowScale,
				flavor % 21);
		bornSub = EXPSUDAKOV(q02, highScale, lowScale, flavor % 21);
	} else {
		sudakov = 1;
		bornSub = 0;
	}
	NAMED_DEBUG("CLUSTERING_STEPS",
			cout << "This gives a sudakov factor of " << sudakov << endl; cout << "born LL sudakov subtraction " << bornSub << std::endl;)
}


double getSudakovFactor(
		fastjet::ClusterSequence& cs,
		const MinloInfo& MI,
		vector<double>& clusteringScales,
		double& Qlocal,
		double& q0,
		double scaleIfFailed,
		double &bornSubtraction,
		int &status,
		const TLorentzVector& basicProcess4Vector,
		bool isReal=false
	){

	int njetsStart=MI.d_njetsOrig;
	int nclusterings=MI.d_njetsOrig - MI.d_njetsClus;
	NAMED_DEBUG("DISPLAY_CS",displayClusterHistory(cs));
    NAMED_DEBUG("DISPLAY_CS_DOT",displayClusterHistoryDot(cs,std::cout));
	if (minloDot::s_MinloWriteDot){
		displayClusterHistoryDot(cs,*minloDot::s_stream);
	}
    //double MEscale2 = Q * Q;

	double factor = 1;
	bornSubtraction =0;

	//cout << "===========" << endl;
	const std::vector<fastjet::ClusterSequence::history_element> & history =
			cs.history();

	vector<double> scales_beamForward;
	vector<double> scales_beamBackward;
	vector<int> pdg_beamForward;
	vector<int> pdg_beamBackward;

	int n = history.size();
	int nparticles = n / 2;
	int njetsCurrent = njetsStart - 1;
	int start = n / 2 + 2; //the two first clusterings are "beam with beam"
	// if this is a real event we skip the first clustering
	if (isReal) {
		start += 1;
	}
	int maxHist = start + nclusterings - 1;

	NAMED_DEBUG("START_AND_END",
		cout << "history size: " << n << endl;
		cout << "njetsStart: " << njetsStart << endl;
		cout << "njetsCurrent: " << njetsCurrent << endl;
		cout << "maxHist: " << maxHist << endl;
		cout << "start: " << start << endl;
	);

	if (isReal) {
	      NAMED_DEBUG("CLUSTERING_STEPS",
	    		  cout << " skipped clustering scale: " << sqrt(history[start-1].dij) << endl;
	      )
	}

	double lastScale2=0; // this is a place where a minimum scale for the recombination scales could be set. If differernt from 0, we don't go lower than that value.
	int additionalExternal=-1;
	int externalToIgnore1=-1;
	int externalToIgnore2=-1;
	const THEPLUGIN::Extras * extras =
			dynamic_cast<const THEPLUGIN::Extras *>(cs.extras());

	if (isReal) {
		int i = start-1;
		int njetsCurrent=njetsStart;

		NAMED_DEBUG("CLUSTERING_STEPS (SKIPPED)",
	    		  cout << " ---- Step " << i << " ----" << endl;
	      	  	  cout <<"dij: "<< cs.history()[i].dij<< endl;
	      	)
    	      int parent1=history[i].parent1;
    	      int parent2=history[i].parent2;
    	      NAMED_DEBUG("CLUSTERING_STEPS (SKIPPED)",
				  cout << " parents: " << parent1 << " " << parent2 << endl;
    	      )

    	      if (parent2 == cs.BeamJet){
    	    	  fastjet::PseudoJet j=cs.jets()[history[parent1].jetp_index];
    	    	  NAMED_DEBUG("BEAM_HISTORY",
    	        	cout <<"jet clustered with beam:" << cs.jets()[history[parent1].jetp_index]<< endl;
    	    	    if (extras->beam_it_clusters_with(j)==+1){
    	          	  cout <<"  forward beam before:" << extras->beam_flav_forward(njetsCurrent+1).description() << endl;
    	          	  cout <<"  forward beam now:" << extras->beam_flav_forward(njetsCurrent).description() << endl;
    	    	    } else {
    	    	    	typedef fastjet::MyFlavKtPlugin THEPLUGIN;  cout <<"  backward beam before:" << extras->beam_flav_backward(njetsCurrent+1).description() << endl;
    	    	      cout <<"  backward beam now:" << extras->beam_flav_backward(njetsCurrent).description() << endl;
    	    	    }
    	    	  )
    	    	  externalToIgnore1=parent1;
    	      } else if ( parent1 >=0 and parent2 >= 0) {
    	    	  double lowScale=history[i].dij;
    	    	  NAMED_DEBUG("CLUSTERING_STEPS",
    	          cout <<"jet clustered together: " << endl
    	        		  << cs.jets()[history[parent1].jetp_index]<< endl
    	        		  << cs.jets()[history[parent2].jetp_index]<< endl;
    	          cout <<"jet scale: " << sqrt(lowScale) << endl;
    	          )
							additionalExternal=i;
							externalToIgnore1=parent1;
							externalToIgnore2=parent2;
    	      }
	}


	double q02 = history[start].dij;
	// this can happen if the second clustering has a smaller dij than the first.
	// it doesn't happen for born type configurations, only for real configurations
	if (q02<lastScale2){
		q02=lastScale2;
	}
	q0 = sqrt(q02);

	scales_beamForward.push_back(q02);
	scales_beamBackward.push_back(q02);


	int fbf = pdgFromFlavor(extras->beam_flav_forward(njetsCurrent + 1));
	int fbb = pdgFromFlavor(extras->beam_flav_backward(njetsCurrent + 1));
	pdg_beamForward.push_back(fbf);
	pdg_beamBackward.push_back(fbb);


	int historyIndex = start;
	//this will store the scales and
	std::vector<sudakovCandidate> scalesConnectedToCoreProcess;
	//this will store the scales and
	std::vector<sudakovCandidate> sudakovCandidates;
	double maxScale2=0;
	bool setQ0ToQlocal=false;
	while (historyIndex <= maxHist) {
      NAMED_DEBUG("CLUSTERING_STEPS",
    		cout << " ---- Step " << historyIndex << " ----" << endl;
  	  	cout <<"dij: "<< cs.history()[historyIndex].dij<< endl;
      )
		double nextScale2=cs.history()[historyIndex].dij;
    if (nextScale2<lastScale2){
			NAMED_DEBUG("CLUSTERING_STEPS",	cout << " -!!! Step " << historyIndex << " has a lower scale "<< sqrt(nextScale2) << " than before "<< sqrt(lastScale2)<< endl;)
			if (MI.d_stopAfterFirstDrop){
				break;
			} else if (MI.d_alltheway){
				nextScale2=lastScale2;
				const fastjet::ClusterSequence::history_element& che=cs.history()[historyIndex];
				fastjet::ClusterSequence::history_element& he=const_cast<fastjet::ClusterSequence::history_element&>(che);
				he.dij=lastScale2;
				NAMED_DEBUG("CLUSTERING_STEPS",cout << " keep last scale (alltheway=true)" << endl;)
			} else {
				lastScale2=nextScale2;
				NAMED_DEBUG("CLUSTERING_STEPS",cout << " didn't keep last scale (alltheway=false)" << endl;)
			}
	} else {
		lastScale2=nextScale2;
	}
    clusteringScales.push_back(sqrt(nextScale2));
    int parent1=history[historyIndex].parent1;
    int parent2=history[historyIndex].parent2;

    if (parent2 == cs.BeamJet){
    	fastjet::PseudoJet j=cs.jets()[history[parent1].jetp_index];
    	NAMED_DEBUG("BEAM_HISTORY",
				cout <<"jet clustered with beam:" << cs.jets()[history[parent1].jetp_index];
				if (extras->beam_clustering_flavour_friendly(j)){
				  cout <<" flavour friendly!" << endl;
				} else {
					cout <<" flavour not friendly!" << endl;
				}
    	)
			// Testing whether we lost the quark to which the vector boson couples
    	bool hasQs=hasQuarks(cs,maxHist-historyIndex);
    	NAMED_DEBUG("CLUSTERING_STEPS",
				if (!hasQs){
				 	cout << " stopping recombination here because there is no quarks to couple the vector boson to." << std::endl;
				}
			)
		  if (!extras->beam_clustering_flavour_friendly(j) || !hasQs){
		    //need to remove the scale we just introduced...
			  clusteringScales.pop_back();
			  NAMED_DEBUG("CLUSTERING_STEPS",
			  	if (!extras->beam_clustering_flavour_friendly(j) ){
						cout << " stopping recombination here because of flavour incompatibility." << std::endl;
					}
			  	if (!hasQs){
						cout << " stopping recombination here because of bad born." << std::endl;
					}
			  )
			  if (historyIndex==start){ // if the first recombination fails then q0 should be the ME scale
				  scales_beamForward.pop_back();
				  scales_beamBackward.pop_back();
				  setQ0ToQlocal=true;
			  }
			  break;
		  }
			// keep track of the highest scale
  		if (nextScale2>maxScale2){
  			maxScale2=nextScale2;
  		}

			if (extras->beam_it_clusters_with(j)==+1){
				scales_beamForward.push_back(nextScale2);
        pdg_beamForward.push_back(pdgFromFlavor(extras->beam_flav_forward(njetsCurrent)));
        NAMED_DEBUG("BEAM_HISTORY",
          cout <<"  forward beam before:" << extras->beam_flav_forward(njetsCurrent+1).description() << endl;
			    cout <<"  forward beam now:" << extras->beam_flav_forward(njetsCurrent).description() << endl;
        )
      } else {
	      scales_beamBackward.push_back(nextScale2);
        pdg_beamBackward.push_back(pdgFromFlavor(extras->beam_flav_backward(njetsCurrent)));
        NAMED_DEBUG("BEAM_HISTORY",
        	cout <<"  backward beam before:" << extras->beam_flav_backward(njetsCurrent+1).description() << endl;
       		cout <<"  backward beam now:" << extras->beam_flav_backward(njetsCurrent).description() << endl;
        )
			}

			if ( history[parent1].parent1 == cs.InexistentParent && history[parent1].parent2 == cs.InexistentParent ){
        	// this is the case if the jet is an initial particle
      }
    } else if ( parent1 >=0 and parent2 >= 0) {
    	double lowScale=nextScale2;
    	//keep track of highest scale
    	if (nextScale2>maxScale2){
  			maxScale2=nextScale2;
  		}
    	NAMED_DEBUG("CLUSTERING_STEPS",
				cout <<"jet clustered together: " << endl
             << cs.jets()[history[parent1].jetp_index]<< endl
				     << cs.jets()[history[parent2].jetp_index]<< endl;
      	cout <<"jet scale: " << sqrt(lowScale) << endl;
			)
			int child=history[historyIndex].child;
			double highScale;
			if (child>maxHist or child == cs.Invalid){
				NAMED_DEBUG("CLUSTERING_STEPS",
					cout <<"jet doesn't get recombined, or too late: " << child << endl;
					cout <<"jet scale of ME: " <<  sqrt(highScale) << endl;
				)
			  sudakovCandidate sc;
			  sc.lowScale=lowScale;
			  sc.flavor=pdgFromFlavor(cs.jets()[history[historyIndex].jetp_index].user_info<fastjet::FlavInfo>());
			  scalesConnectedToCoreProcess.push_back(sc);
			} else {
				highScale=history[child].dij;
        NAMED_DEBUG("CLUSTERING_STEPS",
					cout <<"jet gets recombined at history point: " << child << endl;
					cout <<"jet scale of child: " <<  sqrt(highScale) << endl;
				)
        		sudakovCandidate sc;
				sc.highScale=highScale;
				sc.lowScale=lowScale;
				sc.historyIndex=child;
				sc.flavor=pdgFromFlavor(cs.jets()[history[historyIndex].jetp_index].user_info<fastjet::FlavInfo>());
        sudakovCandidates.push_back(sc);
      }
	    int f=pdgFromFlavor(cs.jets()[history[historyIndex].jetp_index].user_info<fastjet::FlavInfo>());
    }
    historyIndex++;
    njetsCurrent--;
  }



	int nbrClusteringsDone=historyIndex-start;
	minloImpl::g_nclusterings=nbrClusteringsDone;
	NAMED_DEBUG("CLUSTERING_STEPS",
		cout <<"number of clusterings done: " << nbrClusteringsDone << endl;
	)

	std::vector<fastjet::PseudoJet> jetsLeft=cs.exclusive_jets(njetsStart-nbrClusteringsDone);

	if (MI.d_useHT2){
		double ht=sqrt(basicProcess4Vector.Perp2()+80.419*80.419);
		for (int ij=0;ij<jetsLeft.size();ij++){
			fastjet::PseudoJet& j = jetsLeft[ij];
			NAMED_DEBUG("CORE_PROCESS_SCALE",
				std::cout << "now adding jet pt " << j.pt() << std::endl;
				std::cout << "ht so far: " << ht << std::endl;);
				ht+=j.pt();
		}
		NAMED_DEBUG("CORE_PROCESS_SCALE", std::cout << "final core process scale (ht/2): " << ht/2 << " Q^2: "<< ht*ht/4 <<std::endl;);
		Qlocal=ht/2;
	} else {
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
		Qlocal=coreProcess.M();
	}



	double Qlocal2=Qlocal*Qlocal;

	if (raisingScales){
	if (Qlocal2<maxScale2){
		Qlocal2=maxScale2;
		Qlocal=sqrt(Qlocal2);
	}
	}

	if ( setQ0ToQlocal){
		  q0=Qlocal;
		  q02=Qlocal2;
	}


	NAMED_DEBUG("SUDAKOV_CANDIDATES",
		cout << "Sudakovs for scales connected to core process ("<<scalesConnectedToCoreProcess.size()  << "):"<<endl;
		for (int isc=0;isc<scalesConnectedToCoreProcess.size();isc++){
			cout <<  scalesConnectedToCoreProcess[isc] << endl;
		}
		cout << "intermediate Sudakovs ("<< sudakovCandidates.size() << "):" <<endl;
		for (int isc=0;isc<sudakovCandidates.size();isc++){
			cout <<  sudakovCandidates[isc] << endl;
		}
	)

	// need to do this here because I need to know the core scale
	for (int isc=0;isc<scalesConnectedToCoreProcess.size();isc++){
    sudakovCandidate& sc=scalesConnectedToCoreProcess[isc];
		double sudakov,bornSub;//,bornSubNLL;
		computeSudakov(MI,Qlocal2,sc.lowScale,q02,sc.flavor,sudakov,bornSub);
    factor*=sudakov;
    // Keith says we should use the LL one
    bornSubtraction+=bornSub;
	}

	for (int isc=0;isc<sudakovCandidates.size();isc++){
		sudakovCandidate& sc=sudakovCandidates[isc];
		double sudakov,bornSub,bornSubNLL;
		double highScale;
		if (sc.historyIndex>=historyIndex){  // that is the sudakov connects to a point in the history further that we stopped!
			highScale=Qlocal2;
		} else {
			highScale=history[sc.historyIndex].dij ;
			//highScale=sc.highScale;
		}
		computeSudakov(MI,highScale,sc.lowScale,q02,sc.flavor,sudakov,bornSub);
    factor*=sudakov;
    // Keith says we should use the LL one
    bornSubtraction+=bornSub;
	}


	double lastScale;
	if (clusteringScales.size()>0){
		lastScale=clusteringScales.back();
	} else {
		lastScale=Qlocal;
	}
	if (lastScale>Qlocal){
		scales_beamForward.push_back(lastScale*lastScale);
		scales_beamBackward.push_back(lastScale*lastScale);
	} else {
		scales_beamForward.push_back(Qlocal2);
		scales_beamBackward.push_back(Qlocal2);
	}

	NAMED_DEBUG("BEAM_SCALES",
    cout << "Beam forward scales:" ; for (int ii=0;ii<scales_beamForward.size();ii++){ cout << sqrt(scales_beamForward[ii]) << " ";} cout << endl;
    cout << "Beam forward pdg codes:" ; for (int ii=0;ii<pdg_beamForward.size();ii++){ cout << pdg_beamForward[ii] << " ";} cout << endl;
    cout << "Beam backward scales:" ; for (int ii=0;ii<scales_beamBackward.size();ii++){ cout << sqrt(scales_beamBackward[ii]) << " ";} cout << endl;
    cout << "Beam backward pdg codes:" ; for (int ii=0;ii<pdg_beamBackward.size();ii++){ cout << pdg_beamBackward[ii] << " ";} cout << endl;
	  )
	//needs to be adapted!
	if ( nbrClusteringsDone==njetsStart && (pdg_beamForward.back()==21 || pdg_beamBackward.back()==21)){
	  	  NAMED_DEBUG("BORN_CONSISTENCY",
	  	    cout << "born not consistent! " << endl;
	  	  )

		for (int isc=0;isc<clusteringScales.size();isc++){
			clusteringScales[isc]=scaleIfFailed;
		}
	  	status=-1;
	  	bornSubtraction=0;
		minloImpl::g_nclusterings=-1;
		return 1;  // the sudakov
	}

	for (int ii=0;ii<scales_beamForward.size()-1;ii++){
  	double sudakov=nll_sudakov_withInfo(MI,q02,scales_beamForward[ii+1],scales_beamForward[ii],pdg_beamForward[ii] % 21 );
    factor*=sudakov;
    double bornSub=EXPSUDAKOV(q02,scales_beamForward[ii+1],scales_beamForward[ii],pdg_beamForward[ii] % 21 );
    bornSubtraction+=bornSub;
  	NAMED_DEBUG("BEAM_SCALES",
  		cout << "need sudakov " << sudakov << endl;
    	cout << "born sudakov subtraction " << bornSub << endl;
  	)
  }
  NAMED_DEBUG("BEAM_SCALES",    cout << "~~~~~"<< endl;)
  for (int ii=0;ii<scales_beamBackward.size()-1;ii++){
   	double sudakov=nll_sudakov_withInfo(MI,q02,scales_beamBackward[ii+1],scales_beamBackward[ii],pdg_beamBackward[ii] % 21);
		factor*=sudakov;
    double bornSub=EXPSUDAKOV(q02,scales_beamBackward[ii+1],scales_beamBackward[ii],pdg_beamBackward[ii] % 21 );
		bornSubtraction+=bornSub;
		NAMED_DEBUG("BEAM_SCALES",
			cout << "need sudakov from high scale " << scales_beamBackward[ii+1] << " to scale " << scales_beamBackward[ii] << sudakov << endl;
   		cout << "born sudakov subtraction " << bornSub << endl;
 	  )
	}
 	NAMED_DEBUG("EXTERNAL_FACTORS",
    cout << " --- external factors ---- " << endl;
  	)
    // don't take initial state particles (they were already taken into account in the beam sudakovs)
		std::vector<int> externals;
		for (int j=2; j<nparticles;j++){if (j!=externalToIgnore1 && j!=externalToIgnore2){ externals.push_back(j);}}
		// if the first clustering is jet clustering I need to include the result of
		// that clustering in the external particles. Technically I should remove the parents of
		// this node, but they will all lead to sudakovs of 1.
		if (additionalExternal>0){externals.push_back(additionalExternal);}
		NAMED_DEBUG("EXTERNAL_FACTORS",cout <<"external particles: "; for (int iext=0; iext<externals.size();iext++){ cout << " " << externals[iext];};cout << endl;)
    for (int iext=0; iext<externals.size();iext++){
        int j=externals[iext];
        int child=history[j].child;
        double highScale2;
        // the last clustering happened at history index historyIndex-1 if the clustering got interrupted
        // otherwise is maxHist (note if no interruption we have historyIndex-1==maxHist)
        if (child>min(maxHist,historyIndex-1) or child == cs.Invalid){
          	NAMED_DEBUG("EXTERNAL_FACTORS",cout <<"external particle doesn't get recombined, or too late: " << child << endl;)
            highScale2=Qlocal2;
        	NAMED_DEBUG("EXTERNAL_FACTORS",cout <<"  use ME scale for this recombination: " << highScale2 << endl;)
        } else {
        	NAMED_DEBUG("EXTERNAL_FACTORS",cout <<"external particle gets recombined in a jet at history point: " << child << endl;)
            highScale2=history[child].dij;
        	NAMED_DEBUG("EXTERNAL_FACTORS",cout <<"  dij for this recombination: " << highScale2 << endl;)
        }
        int f=pdgFromFlavor(cs.jets()[history[j].jetp_index].user_info<fastjet::FlavInfo>());
        double sudakov,bornSub;
        if (highScale2<q02){ // this happens for the external particles attached to the first clustering in the real radiation
            sudakov=1;
            bornSub=0;
        } else {
        	if (highScale2>q02){
            	sudakov=nll_sudakov_withInfo(MI,q02,highScale2,q02,f % 21 );
                bornSub=EXPSUDAKOV(q02,highScale2,q02,f % 21 );
        	} else {
                sudakov=1;
                bornSub=0;
        	}
        }
        NAMED_DEBUG("EXTERNAL_FACTORS",cout <<"sudakov for external particle "<< j<< ": " <<  sudakov <<" scales (" << sqrt(highScale2)<< ","<< sqrt(q02) << ")" << endl;)
        factor*=sudakov;
        NAMED_DEBUG("EXTERNAL_FACTORS",cout << "born sudakov subtraction " << bornSub << endl;)
        bornSubtraction+=bornSub;
    }
  	NAMED_DEBUG("SUBTRACTION",
  		cout << "total born sudakov subtraction " << bornSubtraction << endl;
    	cout <<"==========="<< endl;
  	)
  	status=nbrClusteringsDone;
  	return factor;

}


double MINLOcomputeSudakov(const MinloInfo& MI,const NtupleInfo<MAX_NBR_PARTICLES>& Ev,double &q0,double &scaleForNLO,int &status,bool useNewNtupleFormat,bool useDouble) {

	std::vector<fastjet::PseudoJet> input_particles;

	NtupleInfo<MAX_NBR_PARTICLES> evb=boostedToCMF(Ev,useDouble);
	NAMED_DEBUG("BOOST", std::cout << "Reconstructed beam energy:" << getBeamEnergy(evb,useDouble) << std::endl; )

	fastjet::PseudoJet beam1(0, 0, +MI.d_energy*evb.x1, MI.d_energy*evb.x1);
	fastjet::PseudoJet beam2(0, 0, -MI.d_energy*evb.x2, MI.d_energy*evb.x2);

	bool isReal = false;
	if (Ev.nparticle-MI.d_njetsOrig-2/*nbr non-partons*/>0 ){ isReal=true;}
	NAMED_DEBUG("EVENT_TYPE",if (isReal ){cout << "This is a real event" << endl;} else {cout << "This is a born-type event" << endl;})


	//this is to be replaced by getFlavourChange to use with old nTuples with missing
	// flavour info, although it didn't reall work...
	if (useNewNtupleFormat){
		beam1.set_user_info(new fastjet::FlavInfo(evb.id1p, fastjet::FlavInfo::beam));
		beam2.set_user_info(new fastjet::FlavInfo(evb.id2p, fastjet::FlavInfo::beam));
	} else {
		beam1.set_user_info(new fastjet::FlavInfo(evb.id1, fastjet::FlavInfo::beam));
		beam2.set_user_info(new fastjet::FlavInfo(evb.id2, fastjet::FlavInfo::beam));
	}
	input_particles.push_back(beam1);
	input_particles.push_back(beam2);


	fillJetVector(evb, input_particles,&Ev.kf[0], true,useDouble);


	NAMED_DEBUG("FLAVOUR_PART", printFlavourPart(input_particles);)
	// specify the jet algorithm
	// flavour kt with alpha=1 (see FlavKtPlugin.hh for other options)

	//int imode = ialg_pp_ktfs1_E;  // try ktf2 and ktLI to see how event changes
	int imode = itype_pp+iangle_DeltaRs+irecom_E;
	double R = MI.d_R;
	int nb;
	if (isReal){
		nb=1;
	} else {
		nb=0;
	}
	fastjet::JetDefinition jet_def = new fastjet::MyFlavKtPlugin(R,MI.d_useModifiedR,nb);
	fastjet::ClusterSequence cs(input_particles, jet_def);
	const THEPLUGIN::Extras * extras = dynamic_cast<const THEPLUGIN::Extras *>(cs.extras());
	NAMED_DEBUG("EVENT_INCL_VIEW", printEventInclView(cs,0.0);)
	// next look at the exclusive clustering, viewing this as a 2-jet event
	// (beware when doing this with "ktfs" type algorithms, which leave
	// clusterings that lead to non-partonic flavours, e.g. ud, until the very end)
	NAMED_DEBUG("EVENT_EXCLUSIV_VIEW", printEventExclView(MI.d_njetsClus,cs);)
	// this is the number of clustering we have to perform
	int nclusterings = MI.d_njetsOrig - MI.d_njetsClus;
	// compute the basic process scale
	TLorentzVector basicProcess;
	// need to use the jets defined by the clustering...
	// need to think about whether or not to apply jet cuts?
	std::vector<fastjet::PseudoJet> fjJets = sorted_by_pt(cs.inclusive_jets());
	for (int j=0;j<MI.d_njetsClus;j++){
		NAMED_DEBUG("CORE_PROCESS_SCALE",
			std::cout << "core process vector: " << basicProcess.E() <<" " << basicProcess.X() << " " << basicProcess.Y() << " " <<  basicProcess.Z() << std::endl;
		)
		basicProcess+=TLorentzVector(
				fjJets[j].px(),
				fjJets[j].py(),
				fjJets[j].pz(),
				fjJets[j].E());
	}
	// needs to be more flexible this will only work for V+jets
	// need to use the boosted momenta!
	if (useDouble){
		basicProcess+=TLorentzVector(
			evb.pxD[0],
			evb.pyD[0],
			evb.pzD[0],
			evb.ED[0]
			);
			} else {
		basicProcess+=TLorentzVector(
			evb.px[0],
			evb.py[0],
			evb.pz[0],
			evb.E[0]);
			}
	NAMED_DEBUG("CORE_PROCESS_SCALE",
		std::cout << "core process vector after first lepton: " << basicProcess.E() <<" " << basicProcess.X() << " " << basicProcess.Y() << " " <<  basicProcess.Z() << std::endl;
	)
	if (useDouble){
		basicProcess+=TLorentzVector(
			evb.pxD[1],
			evb.pyD[1],
			evb.pzD[1],
			evb.ED[1]
			);
	} else {
		basicProcess+=TLorentzVector(
			evb.px[1],
			evb.py[1],
			evb.pz[1],
			evb.E[1]
			);
	}
	NAMED_DEBUG("CORE_PROCESS_SCALE",
		std::cout << "core process vector after second lepton: " << basicProcess.E() <<" " << basicProcess.X() << " " << basicProcess.Y() << " " <<  basicProcess.Z() << std::endl;
	)
	if (extras->hasFSRboost()){
		double boostVec[3];
		extras->getFSRboost(&boostVec[0]);
		basicProcess.Boost(boostVec[0],boostVec[1],boostVec[2]);
		NAMED_DEBUG("CORE_PROCESS_SCALE",
			std::cout << "core process vector after boost: " << basicProcess.E() <<" " << basicProcess.X() << " " << basicProcess.Y() << " " <<  basicProcess.Z() << std::endl;
		)
	}
	if (extras->hasISRboost()){
		double boostVec[3];
		extras->getISRzboost(&boostVec[0]);
		basicProcess.Boost(boostVec[0],boostVec[1],boostVec[2]);
		extras->getISRxyboost(&boostVec[0]);
		basicProcess.Boost(boostVec[0],boostVec[1],boostVec[2]);
		NAMED_DEBUG("CORE_PROCESS_SCALE",
			std::cout << "core process vector after boost: " << basicProcess.E() <<" " << basicProcess.X() << " " << basicProcess.Y() << " " <<  basicProcess.Z() << std::endl;
		)
	}
	double Q=basicProcess.M();
	//take all jets, not only those who passed
	TLorentzVector fullProcess=basicProcess;

	for (int j=MI.d_njetsClus;j<fjJets.size();j++){
		fullProcess+=TLorentzVector(
				fjJets[j].px(),
				fjJets[j].py(),
				fjJets[j].pz(),
				fjJets[j].E());
	}
	double shat=fullProcess.M();

	vector<double> scales;
	double subtraction;
	double Qlocal; // the scale of the core process after doing the clusterings that are allowed
	double sfactor=getSudakovFactor(cs,MI,scales,Qlocal,q0,shat,subtraction,status,basicProcess,isReal);
	NAMED_DEBUG("SUDAKOV_FACTOR",cout << "Factor from sudakovs: " << sfactor << endl;)
	NAMED_DEBUG("SUDAKOV_FACTOR",cout << "born subtraction: " << subtraction << endl;)
	NAMED_DEBUG("ALPHAS_SCALES",cout << "number of clustering scales: " << scales.size();
		for (int si=0; si<scales.size();si++){ cout << " " << scales[si]; }; cout << endl
	)
	double oldAlpha=Ev.alphas;
	double alphaFactor=1;
	double alphaFactorNum=1;
	double alphaFactorDen=1;
	double sumAlpha=0;
	double newAlpha;
	double scaleProduct=1;
	for (int i =0; i<scales.size();i++){
		newAlpha=getAlphasQ(scales[i],MI);
		NAMED_DEBUG("ALPHAS_SCALES",cout << "alphas from clustering scale "<< scales[i] << ": "<< newAlpha << endl;)
		alphaFactor/=oldAlpha;
		alphaFactorDen*=oldAlpha;
		alphaFactor*=newAlpha;
		alphaFactorNum*=newAlpha;
		sumAlpha+=newAlpha;
		scaleProduct*=max(scales[i],minRenScale);
	}
	// this is to set the alphas in the primary process to the ME scale
	for (int i =1;i<=MI.d_njetsClus+MI.d_njetsOrig-scales.size();i++){
		newAlpha=getAlphasQ(Qlocal,MI);
		NAMED_DEBUG("ALPHAS_SCALES",cout << "alphas from primary process scale "<< Qlocal << ": "<< newAlpha << endl;)
		alphaFactor/=oldAlpha;
		alphaFactorDen*=oldAlpha;
		alphaFactor*=newAlpha;
		alphaFactorNum*=newAlpha;
		sumAlpha+=newAlpha;
		scaleProduct*=Qlocal;
	}
	double alphaForNLO=sumAlpha/double(MI.d_njetsOrig);
	if (MI.d_type==MinloInfo::nlo || MI.d_type==MinloInfo::real ){
		alphaFactor/=oldAlpha;
		alphaFactorDen*=oldAlpha;
		alphaFactor*=alphaForNLO;
		alphaFactorNum*=alphaForNLO;
	}
	if (MI.d_scaleMode==MinloInfo::geometric){
		scaleForNLO=exp(log(scaleProduct)/double(MI.d_njetsOrig));
	}
	static QofAlphasInterpolated* QoA=new QofAlphasInterpolated();

	if (MI.d_scaleMode==MinloInfo::inverseAlpha){
		double effectiveAlphasNoSudakov=exp(log(alphaFactor)/double(MI.d_njetsOrig))*oldAlpha;
		scaleForNLO=(*QoA)(effectiveAlphasNoSudakov);
	}

	NAMED_DEBUG("ALPHAS_SCALES",cout <<"scale for NLO: " << scaleForNLO <<  endl;)
	NAMED_DEBUG("ALPHAS_SCALES",cout <<"ME scale: " << Q <<  " original: "<< Ev.muR << endl;)
	NAMED_DEBUG("ALPHAS_SCALES",cout <<"alphas factor: " << alphaFactor <<  endl;)
	NAMED_DEBUG("ALPHAS_SCALES",cout <<"alphas factor (numerator): " << alphaFactorNum <<  endl;)
	NAMED_DEBUG("ALPHAS_SCALES",cout <<"alphas factor (denominator): " << alphaFactorDen <<  endl;)

	if (MI.d_type==MinloInfo::born){
		double b0=(33-2*5)/12.0/3.14159265358979323846264;
		//double fullSubtraction= alphaForNLO*(subtraction+d_njetsOrig*b0*2*log(scaleForNLO/fixedScaleForNLO));
		double fullSubtraction=
			alphaForNLO*(
				subtraction+MI.d_njetsOrig*b0*2*log(scaleForNLO/Ev.muR)
				);
		double bornFactor=(1+fullSubtraction);
		NAMED_DEBUG("ALPHAS_SCALES",cout <<"full subtraction: " << fullSubtraction <<  endl;)
		NAMED_DEBUG("ALPHAS_SCALES",cout <<"full born factor: " << bornFactor <<  endl;)
		NAMED_DEBUG("KEITH",cout << "To compare: basicfac: " << sfactor*alphaFactor << endl;)
		NAMED_DEBUG("KEITH",cout << "To compare: bornfac: " << bornFactor << endl;)
		NAMED_DEBUG("KEITH",cout << "To compare: nlofac: " << alphaForNLO/oldAlpha << endl;)
		NAMED_DEBUG("KEITH",cout << "To compare: st_mufac2: " << q0*q0 << endl;)
		double factor=sfactor*alphaFactor*bornFactor;
		NAMED_DEBUG("KEITH",cout << "To compare: weight: " << factor << endl;)
		double effectiveAlphas=exp(log(sfactor*alphaFactor)/double(MI.d_njetsOrig))*oldAlpha;  // this spreads the factor to all alphas, not including the subtraction
		minloImpl::g_MinloFactor=sfactor*alphaFactor;
		minloImpl::g_MinloScale=(*QoA)(effectiveAlphas);
		minloImpl::g_MinloScaleValid=true;
		return factor;

		//return sfactor*alphaFactor*(1-alphaForNLO*subtraction);
	} else {
		double b0=(33-2*5)/12.0/3.14159265358979323846264;
		double fullSubtraction= alphaForNLO*(subtraction+MI.d_njetsOrig*b0*2*log(scaleForNLO/Ev.muR));
		double bornFactor=(1+fullSubtraction);
		double nlofac=alphaForNLO/oldAlpha;

		NAMED_DEBUG("ALPHAS_SCALES",
			cout <<"full subtraction: " << fullSubtraction <<  endl;
			cout <<"full born factor: " << bornFactor <<  endl;
		)
		NAMED_DEBUG("KEITH",cout << "To compare: basicfac: " << sfactor*(alphaFactor/nlofac) << endl;) // in this case my alphafactor include nlofac already, need to divide out to compare
		NAMED_DEBUG("KEITH",cout << "To compare: bornfac: " << bornFactor << endl;)
		NAMED_DEBUG("KEITH",cout << "To compare: nlofac: " << nlofac << endl;)
		NAMED_DEBUG("KEITH",cout << "To compare: st_mufac2: " << q0*q0 << endl;)
		NAMED_DEBUG("KEITH",cout << "To compare: weight: " << sfactor*alphaFactor << endl;)
		return sfactor*alphaFactor;
	}

}

NtupleInfo<MAX_NBR_PARTICLES> boostedToCMF(const NtupleInfo<MAX_NBR_PARTICLES>& orig,bool useDouble){
	NtupleInfo<MAX_NBR_PARTICLES> evb;
	evb.nparticle=orig.nparticle;

	double LabEx=(orig.x1+orig.x2);
	double Labzx=(orig.x1-orig.x2);
	double beta=-Labzx/LabEx;
	double ECMFx=sqrt(orig.x1*orig.x2);
	std::vector<TLorentzVector> boostedVectors;
	evb.x1=ECMFx;
	evb.x2=ECMFx;
	evb.id1=orig.id1;
	evb.id2=orig.id2;
	evb.id1p=orig.id1p;
	evb.id2p=orig.id2p;
	if (useDouble){
		for (int ii=0;ii<orig.nparticle;ii++){
			TLorentzVector boosted(orig.pxD[ii],orig.pyD[ii],orig.pzD[ii],orig.ED[ii]);
			boosted.Boost(0,0,beta);
			evb.pxD[ii]=boosted.X();
			evb.pyD[ii]=boosted.Y();
			evb.pzD[ii]=boosted.Z();
			evb.ED[ii]=boosted.E();
			evb.kf[ii]=orig.kf[ii];
		}
	} else {
		for (int ii=0;ii<orig.nparticle;ii++){
			TLorentzVector boosted(orig.px[ii],orig.py[ii],orig.pz[ii],orig.E[ii]);
			boosted.Boost(0,0,beta);
			evb.px[ii]=boosted.X();
			evb.py[ii]=boosted.Y();
			evb.pz[ii]=boosted.Z();
			evb.E[ii]=boosted.E();
			evb.kf[ii]=orig.kf[ii];
		}
	}
	return evb;
};


double getBeamEnergy(const NtupleInfo<MAX_NBR_PARTICLES>& Ev,bool useDouble){
	TLorentzVector sum(0,0,0,0);
	for (int ii=0;ii<Ev.nparticle;ii++){
		if (useDouble){
			sum+=TLorentzVector(Ev.pxD[ii],Ev.pyD[ii],Ev.pzD[ii],Ev.ED[ii]);
		} else {
			sum+=TLorentzVector(Ev.px[ii],Ev.py[ii],Ev.pz[ii],Ev.E[ii]);
		}
	}
	double E=sum.E();
	double Z=sum.Z();
	if (sum.X()/E> 1e-6){
		std::cerr << "suspicious momentum conservation! sum.X()/E=" << sum.X()/E << std::endl;
	}
	if (sum.Y()/E> 1e-6){
		std::cerr << "suspicious momentum conservation! sum.Y()/E" << sum.Y()/E << std::endl;
	}

	double ECMF=sqrt((E*E-Z*Z)/Ev.x1/Ev.x2);
	return ECMF/2;
};

void printEventInclView(fastjet::ClusterSequence& cs,double ptCut){
	vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets(ptCut));
	cout << endl;
	cout << "Inclusive jets  " << endl;
	cout << "------------------------------" << endl;
	print_header(cout) << endl;
	for (unsigned i = 0; i < jets.size(); i++) cout << jets[i] << endl;
}

void printEventExclView(int n_excl,fastjet::ClusterSequence& cs){

	  vector<fastjet::PseudoJet> jets = cs.exclusive_jets(n_excl);
	  cout << endl;
	  cout << "Event viewed as " << n_excl << " exclusive jets" << endl;
	  cout << "--------------------------------" << endl;
	  print_header(cout) << endl;
	  for (unsigned i = 0; i < jets.size(); i++) cout << jets[i] << endl;
	  // to get the beam flavours in the exclusive view, we need "extras" from the cs;
	  const THEPLUGIN::Extras * extras = dynamic_cast<const THEPLUGIN::Extras *>(cs.extras());
	  if (extras != 0) {
		  cout << "   forward  beam flavour: " << extras->beam_flav_forward(n_excl).description()  << endl;
		  cout << "   backward beam flavour: " << extras->beam_flav_backward(n_excl).description() << endl;
		  cout << endl;
	  }

	  cout << endl;

	  for (unsigned i = 0; i < jets.size(); i++){
		  //printHistory(jets[i],0);
	  }
}

void printFlavourPart(const std::vector<fastjet::PseudoJet>& input_particles){
	  // output the flavoured part of the event\
	  cout << "All flavoured input particles " << endl;
	  cout << "------------------------------" << endl;
	  cout << "index ";
	  print_header(cout) << endl;
	  double Etot = 0;
	  for (unsigned i = 0; i < input_particles.size(); i++) {
	    Etot += input_particles[i].E();
//	    if (input_particles[i].user_info<FlavInfo>().is_flavourless()) continue;
	    std::cout << std::setw(5) << i << " " << input_particles[i] << endl;
	  }

}

const int nMomenta=20;


double MINLO_computeSudakovKeith(
	const NtupleInfo<MAX_NBR_PARTICLES>& Ev,
	int flg_bornonly,
	int imode,
	int isReal,
	double BeamEnergy,
	int nlegborn,
	int st_bornorder,
	bool useMinloIDs,
	bool useDouble,
	int alltheway
	){

	double kn_pborn[4*nMomenta];
	double kn_cmpborn[4*nMomenta];
	int flav[nMomenta];
	int index;

	double E=BeamEnergy;


	kn_pborn[0]=E*Ev.x1;
	kn_pborn[1]=0;
	kn_pborn[2]=0;
	kn_pborn[3]=E*Ev.x1;
	kn_pborn[4]=E*Ev.x2;
	kn_pborn[5]=0;
	kn_pborn[6]=0;
	kn_pborn[7]=-E*Ev.x2;
	if (useMinloIDs){
		flav[0]=Ev.id1p%21;
		flav[1]=Ev.id2p%21;
	} else {
		flav[0]=Ev.id1%21;
		flav[1]=Ev.id2%21;
	}
	double LabE=E*(Ev.x1+Ev.x2);
	double Labz=E*(Ev.x1-Ev.x2);
	double beta=-Labz/LabE;

	double ECMF=E*sqrt(Ev.x1*Ev.x2);

	kn_cmpborn[0]=ECMF;
	kn_cmpborn[1]=0;
	kn_cmpborn[2]=0;
	kn_cmpborn[3]=ECMF;
	kn_cmpborn[4]=ECMF;
	kn_cmpborn[5]=0;
	kn_cmpborn[6]=0;
	kn_cmpborn[7]=-ECMF;

	for (int i=0;i<Ev.nparticle;i++){
//		index=(nMomenta)*0+(i+2);
//		kn_pborn[index]=Ev.E[i];
//		index=(nMomenta)*1+(i+2);
//		kn_pborn[index]=Ev.px[i];
//		index=(nMomenta)*2+(i+2);
//		kn_pborn[index]=Ev.py[i];
//		index=(nMomenta)*3+(i+2);
//		kn_pborn[index]=Ev.pz[i];
		if (useDouble){
			kn_pborn[i*4+8]=Ev.ED[i];
			kn_pborn[i*4+9]=Ev.pxD[i];
			kn_pborn[i*4+10]=Ev.pyD[i];
			kn_pborn[i*4+11]=Ev.pzD[i];
		} else {
			kn_pborn[i*4+8]=Ev.E[i];
			kn_pborn[i*4+9]=Ev.px[i];
			kn_pborn[i*4+10]=Ev.py[i];
			kn_pborn[i*4+11]=Ev.pz[i];
		}


		flav[i+2]=Ev.kf[i]%21;

		if (useDouble){
			TLorentzVector boosted(Ev.pxD[i],Ev.pyD[i],Ev.pzD[i],Ev.ED[i]);
			boosted.Boost(0,0,beta);
			kn_cmpborn[i*4+8]=boosted.E();
			kn_cmpborn[i*4+9]=boosted.X();
			kn_cmpborn[i*4+10]=boosted.Y();
			kn_cmpborn[i*4+11]=boosted.Z();
		} else {
			TLorentzVector boosted(Ev.px[i],Ev.py[i],Ev.pz[i],Ev.E[i]);
			boosted.Boost(0,0,beta);
			kn_cmpborn[i*4+8]=boosted.E();
			kn_cmpborn[i*4+9]=boosted.X();
			kn_cmpborn[i*4+10]=boosted.Y();
			kn_cmpborn[i*4+11]=boosted.Z();
		}
	}
    //int      nlegborn=6;      //  ! Number of final state particle in Born plus 2.
    //int      st_bornorder=2;    //! Number of powers of aS in the Born (ZJ=1, ZJJ=2)
//    double      st_muren2=fixedScaleForNLO*fixedScaleForNLO;       //! The mu_R^2 value that was used to evaluate virt
    double      st_muren2=Ev.muR*Ev.muR;       //! The mu_R^2 value that was used to evaluate virt
//    double      st_muren2=Ev.muR*Ev.muR;       //! The mu_R^2 value that was used to evaluate virt
    double      rad_charmthr2=2.25;   //! Charm  (mass) threshold squared: (1.5 GeV)^2
    double      rad_bottomthr2=25.0;  //! Bottom (mass) threshold squared: (5.0 GeV)^2
    int      st_nlight=5;       //! Nominal number of light flavs in computation (5)
    double      st_lambda5MSB=0.226 /*same as in sudakov.f*/;   //! Lambda_5 in the MS bar scheme
    double      st_renfact=1.0;      //! Scaling factor for mu_R (not mu_R^2): 0.5,1,2.0
    double      st_facfact=1.0;      //! Scaling factor for mu_F (not mu_F^2): 0.5,1,2.0
    int      raisingscales=1;   //! if set to 1: Freeze clustering scale if it attempts to
                           //! go down - no more sudakovs/subtractions, only
                           //! coupling reweighting (for all remaining nodes).
//			C - Momenta tha tare only used for debugging
//    double*      kn_pborn,        //! Born particles momenta in LAB frame
//			C - Old inputs and output
//    double*      kn_cmpborn,      //! Born particles momenta in C.O.M frame
//    int*      flav,            //! Born particle flavours (1,2) is incoming
    double      rescfac;         //! Comb. of suds, subtractions & aS weights to
                           //! rescale the weight of the contrib assoc. to this
                           //! mom. and flav. configuration.
    double      basicfac;        //! Comb. of suds, subtractions & aS weights to
    double      bornfac;         //! Comb. of suds, subtractions & aS weights to
    double      nlofac;          //! Comb. of suds, subtractions & aS weights to
//			C - Extra factorisation scale output as it normally lives in pwhg common block
    double      st_mufact2;      //! The mu_F^2 value that the PDFs should be
                            //! re-evaluated at now

    double q2=Ev.muR*Ev.muR;
	#ifdef LHAPDF_NEW_VERSION
    	double alphasMZ=currentPDF::s_PDF->alphasQ(fixedScaleForNLO);
	#else
		double alphasMZ=LHAPDF::alphasPDF(fixedScaleForNLO);
	#endif
	double alphaspwg=I_PWHG_ALPHAS(q2,st_lambda5MSB,st_nlight, rad_charmthr2, rad_bottomthr2);
    double mz2=fixedScaleForNLO*fixedScaleForNLO;
	double alphaspwgMZ=I_PWHG_ALPHAS(mz2,st_lambda5MSB,st_nlight, rad_charmthr2, rad_bottomthr2);

	    double      st_alpha=Ev.alphas;        //! Fixed aS used in MEs before they got here
//	    double      st_alpha=alphaspwg;        //! Fixed aS used in MEs before they got here
	//    double      st_alpha=alphaspwgMZ;        //! Fixed aS used in MEs before they got here

    NAMED_DEBUG("KEITH",cout << "Original scale : " << st_muren2 << " Ev.muR" << Ev.muR << endl;)
	NAMED_DEBUG("KEITH",cout << "nlegborn:    " << nlegborn << endl;)

	double res=I_SETLOCALSCALES(
//			C - New extra inputs for independence from includes
	           &flg_bornonly,    //! Are we feeding through only Born stuff, no NLO?
	           imode,           //! imode=1 for Born, imode=2 for all NLO contribs
			   isReal,         // ! Set isReal=1 for real kinematics, 0 otherwise.
			   nlegborn,      //  ! Number of final state particle in Born plus 2.
	           kn_cmpborn,      //! Born particles momenta in C.O.M frame
//	           kn_pborn,      //! Born particles momenta in C.O.M frame
	           flav,            //! Born particle flavours (1,2) is incoming
	           raisingscales,   //! Freeze clustering scale if it attempts to
	                            //! go down - no more sudakovs/subtractions, only
	                            //! coupling reweighting (for all remaining nodes).
		           st_alpha,        //! Fixed aS used in MEs before they got here
		           st_muren2,       //! The mu_R^2 value that was used to evaluate virt
		           st_bornorder,    //! Number of powers of aS in the Born (ZJ=1, ZJJ=2)
		           rad_charmthr2,   //! Charm  (mass) threshold squared: (1.5 GeV)^2
		           rad_bottomthr2,  //! Bottom (mass) threshold squared: (5.0 GeV)^2
		           st_nlight,       //! Nominal number of light flavs in computation (5)
		           st_lambda5MSB,   //! Lambda_5 in the MS bar scheme
		           st_renfact,      //! Scaling factor for mu_R (not mu_R^2): 0.5,1,2.0
		           st_facfact,      //! Scaling factor for mu_F (not mu_F^2): 0.5,1,2.0
		           rescfac,         //! Comb. of suds, subtractions & aS weights to
		                            //! rescale the weight of the contrib assoc. to this
		                            //! mom. and flav. configuration.
		           st_mufact2,      //! The mu_F^2 value that the PDFs should be
		                             //! re-evaluated at now

				   //			C - Momenta tha tare only used for debugging
//			C - Old inputs and output
		           basicfac,        //! Comb. of suds, subtractions & aS weights to
		           bornfac,         //! Comb. of suds, subtractions & aS weights to
		           nlofac,          //! Comb. of suds, subtractions & aS weights to
//			C - Extra factorisation scale output as it normally lives in pwhg common block
//		           kn_pborn,        //! Born particles momenta in LAB frame
                   alltheway

	);
	NAMED_DEBUG("KEITH",cout << "Keith basicfac:  " << basicfac << endl;)
	NAMED_DEBUG("KEITH",cout << "Keith bornfac:   " << bornfac << endl;)
	NAMED_DEBUG("KEITH",cout << "Keith nlofac:    " << nlofac << endl;)
	NAMED_DEBUG("KEITH",cout << "Keith st_mufact2:" << st_mufact2 << endl;)

	return res;
}


extern "C" {

	double PDFFROMLHAPDF (double& Q){
#ifdef LHAPDF_NEW_VERSION
		return currentPDF::s_PDF->alphasQ(Q);
#else
		return LHAPDF::alphasPDF(Q);
#endif
	}
}

