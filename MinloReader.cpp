#include "MinloReader.h"
#include "MINLOfunctions.h"
#include <string>
#include <vector>

double MINLOreader::computeSudakov(MinloInfo& MI, int weightType,double &q0,double &scaleForNLO){
	return MINLOcomputeSudakov(MI,d_NI,weightType,q0,scaleForNLO,d_hasMinlo);
};
void MINLOreader::addFiles(const std::vector<std::string>& fs){
    RootFileReaderBase::addFiles(d_NI,fs);
    d_hasMinlo=nTupleHasMinlo(getVersion());

  };

MINLOreader::MINLOreader(){
  RootFileReaderBase::init(d_NI,"t3");
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
