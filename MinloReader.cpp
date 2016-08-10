#include "MinloReader.h"
#include "MINLOfunctions.h"
#include <string>
#include <vector>
#include "pdf.h"

double MINLOreader::computeSudakov(const MinloInfo& MI,double &q0,double &scaleForNLO){
	return MINLOcomputeSudakov(MI,d_NI,q0,scaleForNLO,d_hasMinlo);
};
void MINLOreader::addFiles(const std::vector<std::string>& fs){
    RootFileReaderBase::addFiles(d_NI,fs);
    d_hasMinlo=nTupleHasMinlo(getVersion());

  };

MINLOreader::MINLOreader(){
//  RootFileReaderBase::init(d_NI);
}


double MINLOreader::computeSudakovKeith(const MinloInfo& MI,const KeithInfo& KI){
		int isReal;
		if ( d_NI.nparticle==(KI.nlegborn-2) ) {
			isReal=0; //! Set isReal=1 for real kinematics, 0 otherwise.
		} else {
			isReal=1;          //! Set isReal=1 for real kinematics, 0 otherwise.
		}
		return MINLO_computeSudakovKeith(d_NI,KI.flg_bornonly,KI.imode,isReal,MI.d_energy,KI.nlegborn,KI.st_bornorder,d_hasMinlo);
};

std::vector<double> MINLOreader::momentum(int i) {
    std::vector<double> p(4);
	  p[0]=d_NI.E[i];
    p[1]=d_NI.px[i];
    p[2]=d_NI.py[i];
    p[3]=d_NI.pz[i];
    return p;
}

void MINLOreader::initPDF(const std::string& pdfName){
	currentPDF::init(pdfName,0);
}
