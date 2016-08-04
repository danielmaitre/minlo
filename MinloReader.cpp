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
