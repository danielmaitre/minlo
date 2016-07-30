#include "MinloReader.h"
#include "MINLOfunctions.h"
#include "pdf.h"
using namespace std;




int main(){

  currentPDF::init("CT10nlo",0);

	std::vector<std::string> fs;
	fs.push_back("data/Wm2j_born_small.root");

	MINLOreader r;

	r.addFiles(fs);

	int weightType=1;

	MinloInfo MI;
	// number of jets at the beginning
	MI.d_njetsOrig=2;
	// number of jets to cluster to
	MI.d_njetsClus=0;
	MI.d_energy=3500;
	MI.d_type=MinloInfo::born;
	
	while(r.nextEntry()){
		int id=r.d_NI.id;

		double q0,scaleForNLO;
		double alphaFactor=MINLOcomputeSudakov(MI,r.d_NI,weightType,q0,scaleForNLO);


		
		int flg_bornonly=0;    //! Are we feeding through only Born stuff (1), or NLO (0)?
		int imode=1;           //! imode=1 for Born, imode=2 for all NLO contribs
		int isReal=0;          //! Set isReal=1 for real kinematics, 0 otherwise.
		double keith = MINLO_computeSudakovKeith(r.d_NI,flg_bornonly,imode,isReal);
		cout << "Keith weight: " << keith << endl;
		cout << "(" << r.d_NI.id<<") Keith weight/my weight: " << keith/alphaFactor << endl;

		
	}

	return 0;

}
