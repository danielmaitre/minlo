
#ifndef MINLOREADER_H_
#define MINLOREADER_H_

#include "ntuplereader/nTupleReader.h"
#include "ntuplereader/NtupleInfo.h"
#include "ntuplereader/EventReaderBase.h"

class MINLOreader : public RootFileReaderBase {
public:
  NtupleInfo<MAX_NBR_PARTICLES> d_NI;
  MINLOreader() {
    RootFileReaderBase::init(d_NI,"t3");
  }
  bool nextEntry(){return readNextEntry(d_NI);}
  void addFiles(const std::vector<std::string>& fs){
  		RootFileReaderBase::addFiles(d_NI,fs);
  };
  virtual ~MINLOreader(){};
};


#endif /* MINLOREADER_H_*/

