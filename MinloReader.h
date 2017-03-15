
#ifndef MINLOREADER_H_
#define MINLOREADER_H_

#include "ntuplereader/nTupleReader.h"
#include "ntuplereader/NtupleInfo.h"
#include "ntuplereader/EventReaderBase.h"
#include "ntuplereader/version.h"
#include "MinloInfo.h"


class MINLOreader : public nTupleReader {
  bool d_hasMinlo;
  bool d_useDouble;
  const NtupleInfo<MAX_NBR_PARTICLES>& d_NI;
  public:
  MINLOreader() ;
  double computeSudakov(const MinloInfo& MI,double &q0,double &scaleForNLO,int &status);
  double computeSudakovKeith(const MinloInfo& MI,const KeithInfo& KI);
  std::vector<double> momentum(int i);
  void initPDF(const std::string& pdfName);
  virtual void addFiles(const std::vector<std::string>& fs);
  virtual ~MINLOreader(){};
};

#endif /* MINLOREADER_H_*/

