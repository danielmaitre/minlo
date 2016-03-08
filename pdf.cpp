#include "pdf.h"
#include "LHAPDF/LHAPDF.h"
#include <iostream>

using namespace std;


#ifdef LHAPDF_NEW_VERSION

LHAPDF::PDF* currentPDF::s_PDF=0;
LHAPDF::PDFSet* currentPDF::s_PDFSet=0;

void currentPDF::init(const std::string& name,int set){
  s_PDFSet=new LHAPDF::PDFSet("CT10nlo");
  s_PDF=s_PDFSet->mkPDF(set);
}

#endif



