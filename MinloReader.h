
#ifndef MINLOREADER_H_
#define MINLOREADER_H_

#include "ntuplereader/nTupleReader.h"
#include "ntuplereader/NtupleInfo.h"
#include "ntuplereader/EventReaderBase.h"
#include "ntuplereader/version.h"
#include "MinloInfo.h"


class MINLOreader : public RootFileReaderBase {
  bool d_hasMinlo;
  public:
  NtupleInfo<MAX_NBR_PARTICLES> d_NI;
  MINLOreader() ;
  bool nextEntry(){return readNextEntry(d_NI);}
  void addFiles(const std::vector<std::string>& fs);
  double computeSudakov(const MinloInfo& MI,double &q0,double &scaleForNLO);
  double computeSudakovKeith(const MinloInfo& MI,const KeithInfo& KI);
  std::vector<double> momentum(int i);
  void initPDF(const std::string& pdfName);
  virtual ~MINLOreader(){};
};


#endif /* MINLOREADER_H_*/

